import numpy as np
import xarray as xr 
import pandas as pd
from os.path import isfile
import sys
from dateutil.parser import parse
import searoute as sr
from scgraph.geographs.marnet import marnet_geograph
from scipy.stats import pearsonr

#Purpose of the script:
#    Data from observation like TideGauges are possibly sampled 
#    in time not evenly or there are gaps in the entries.
#
#    This script takes a tide gauge input from CMEMS catalog, identifies
#    the main sampling rate, resamples it and fills the gaps with masked
#    values

#ds["TIME"].values[:]  = np.array(pd.to_datetime(ds["TIME"]).round("1s"))

def resampleCmemsTG(dsIN=xr.Dataset,
                    pathOUT="out.nc",
                    rule=None,
                    verbose=False,
                    clipDates=None,
                    makePlot=False):

    #
    # this routine resamples the TG timeseries from CMEMS catalog 
    # using the most common sampling frequency within the series
    # (the same historical timeseries might contain more than 1 sampling
    # frequency)
    #

    # if you don't want to do a time-slice just set clipYear like = (1800,2100)

    if type(dsIN) is not xr.Dataset:
        raise Exception("Error: the input dsIN must be an xr.Dataset")

    ds  = dsIN.copy()

    if clipDates != None:
        if verbose:
            print("Slicing dataset for dates (%s,%s)" % (clipDates[0],clipDates[1]))
        date0 = parse(clipDates[0])
        date1 = parse(clipDates[1])
        ds = ds.sel(TIME=slice(date0,date1))

    if verbose:
        print("sizes of TG Dataset: ", ds.sizes)

    try:
        if ds.sizes["DEPTH"] == 2:
            print("WARNING! Found dataset with DEPTH = 2, selecting DEPTH=0")
            ds = ds.isel(DEPTH = 0)
    except:
        pass

    values      = ds["SLEV"]
    timestamps  = ds["TIME"]
    qc          = ds["SLEV_QC"]

    # I take only qc=1 values
    goodValues  = qc == 1

    ds.close()

    # Create a pandas Series with the timestamps as the index
    try:
        time_series = pd.Series(values, index=pd.to_datetime(timestamps))
    except:    
        time_series = pd.Series(values.squeeze(), index=pd.to_datetime(timestamps))     # sometimes is SLEV[TIME,DEPH]

    # set not-good values to Nan
    try:
        time_series[goodValues.data == False] = np.nan
    except:    
        time_series[goodValues.data.squeeze() == False] = np.nan

    # Calculate time differences between consecutive samples
    time_diffs = time_series.index.to_series().diff().dropna()

    # Round the time differences to a certain precision (e.g., seconds)
    rounded_time_diffs = time_diffs.dt.round('1s')

    # replace index of timeseries with rounded time index (it is necessary
    # otherwise some good values are set to NaN)
    time_series.index = time_series.index.round("1s")

    if rule == None:
        # Determine the most common time difference
        new_sampling_freq = rounded_time_diffs.mode().iloc[0]
        if verbose:
            print("Resampling using the most common time difference:", new_sampling_freq) 
    else:
        # check if rule is str of Timedelta
        if verbose:
            print("Resampling using ", rule)    
        new_sampling_freq = rule

    offSet  = 0.5 * new_sampling_freq

    # Resample the time series using the main sampling frequency 
    resampled_series = time_series.resample(new_sampling_freq).asfreq()

    # Fill missing values with masked values
    #masked_series = resampled_series.where(~resampled_series.isna(), np.nan)

    daOut = xr.DataArray(resampled_series.values, 
                        dims = ["TIME"], 
                        coords = {"TIME" : resampled_series.index, 
                                  "LONGITUDE": ds["LONGITUDE"].values.item(), 
                                  "LATITUDE" : ds["LATITUDE"].values.item() 
                                  } 
                        )
         
    daOut.name = "SLEV"

    # return DataArray
    return daOut

def DailyAverageCmemsTG(daIN=xr.DataArray,
                        pathOUT="out.nc",
                        verbose=False,
                        clipDates=None,         # list / tuple of (Year0,Year1) to clip dataset
                        anomaly=True,           # if True, the total mean of selected period is subtracted 
                        nMin=80,                # minimum number [%] of values per bin to calculate daily average
                        reindex=True,           # round the legth of timseries to whole (start/end) month if necessary  
                        makePlot=False
                        ):

    print(type(daIN))

    if type(daIN) is xr.Dataset:
        try:
            daIN    = daIN["SLEV"]
        except:            
            raise Exception("Error: input dataset does not contain SLEV")

    if type(daIN) is not xr.DataArray:
        raise Exception("Error: input is not a xr.DataArray nor xr.Dataset")

    da = daIN.copy()

    if clipDates != None:
        if verbose: 
            print("Slicing DataArray for Dates (%s,%s)" % (clipDates[0],clipDates[1]))
        date0 = parse(clipDates[0])
        date1 = parse(clipDates[1])
        da = da.sel(TIME=slice(date0,date1))

    values      = da.values
    timeStamps  = da["TIME"]

    # Create a pandas Series with the timestamps as the index
    timeSeries = pd.Series(values, index=pd.to_datetime(timeStamps))

    # create a MultiIndex to properly groupby
    multiIndex = pd.MultiIndex.from_arrays([timeSeries.index.year,timeSeries.index.month,timeSeries.index.day])

    # calculate mean only if number of counts per bin is > (nMin / 100) * nExpected 
    sampligFreq = timeSeries.index.to_series().diff().median()
    nExpected   = int(pd.Timedelta("1D") / sampligFreq)
    ratio       = float(nMin) / 100
    minNValues  = int(ratio * nExpected) 

    if verbose:
        print("Calculating daily mean if bins have [>=%d%%] of values" % nMin)

    # do temporary group mean     
    groupedMeanTmp = timeSeries.groupby(by=multiIndex).agg(["mean","count"])

    # filter group mean
    groupedMeanTmp["mean"] = groupedMeanTmp["mean"].where(groupedMeanTmp["count"]>=minNValues)

    groupedMean = groupedMeanTmp["mean"]

    # change the Index (made of YYYY,mm,dd) to pd.DatetimeIndex
    groupedMean.index = pd.to_datetime(groupedMean.index.map(lambda x: '-'.join(map(str, x)))) 

    # add 12h offset
    groupedMean.index   += pd.Timedelta("12h")

    """
    # create dates to re-index (full dates) 
    fullDateRange   = pd.date_range(start=pd.Timestamp(groupedMean.index.year.min(),1,1),
                                    end=pd.Timestamp(groupedMean.index.year.max(), 12, 31), 
                                    freq='D')
    
    fullDateRange   += pd.Timedelta("12h")

    # re-index the series
    groupedMean     = groupedMean.reindex(fullDateRange)
    """
    
    if reindex:
        if clipDates == None:
            raise Exception("Error: you must set clipDates to reindex!")

        fullDateStart   = parse(clipDates[0])                                   
        fullDateEnd     = parse(clipDates[1])     

        fullDateRange   = pd.date_range(start=fullDateStart,end=fullDateEnd,freq="D") 
        fullDateRange   += pd.Timedelta("12h")

        """
        print("     reindex dbg: ")
        print(groupedMean.index[0])
        print(groupedMean.index[-1])
        print(fullDateRange[0])
        print(fullDateRange[-1])
        """
        
        if verbose:
            print("   Reindexing to : ")
            print("           start : %s" % fullDateRange[0])
            print("             end : %s" % fullDateRange[-1])
            groupedMean     = groupedMean.reindex(fullDateRange)

    # remove anomaly
    if anomaly:
        groupedMean.values[:] = groupedMean.values[:] - groupedMean.mean()

    # create output DataArray
    daOut = xr.DataArray(groupedMean.values, 
                        dims = ["TIME"], 
                        coords = {"TIME" : groupedMean.index, 
                                  "LONGITUDE": da["LONGITUDE"].values.item(), 
                                  "LATITUDE" : da["LATITUDE"].values.item() 
                                  } 
                        )
         
    daOut.name = "SLEV"

    # return DataArray
    return daOut

def DailyAverageNemoSSH(daIN=xr.DataArray,
                        pathOUT="out.nc",
                        verbose=False,
                        clipDates=None,             # list / tuple of (Year0,Year1) to clip dataset
                        anomaly=False,              # if True, the total mean of selected period is subtracted
                        reindex=True,               # round the length of timeseries to the whole (start/end) month if necessary
                        makePlot=False):
    
    #
    # compute Daily average of ssh from native NEMO output
    #

    if type(daIN) is xr.Dataset:
        try:
            daIN    = daIN["ssh"]
        except:            
            raise Exception("Error: input dataset does not contain ssh")

    if type(daIN) is not xr.DataArray:
        raise Exception("Error: input is not a xr.DataArray nor xr.Dataset")

    da = daIN.copy()

    if clipDates != None:
        if verbose: 
            print("Slicing DataArray for Dates (%s,%s)" % (clipDates[0],clipDates[1]))
        date0 = parse(clipDates[0])
        date1 = parse(clipDates[1])
        da = da.sel(time_counter=slice(date0,date1))

    values      = da.values
    timeStamps  = da["time_counter"]

    # Create a pandas Series with the timestamps as the index
    timeSeries = pd.Series(values.squeeze(), index=pd.to_datetime(timeStamps))

    # create a MultiIndex to properly groupby
    multiIndex = pd.MultiIndex.from_arrays([timeSeries.index.year,timeSeries.index.month,timeSeries.index.day])

    # do the group mean    
    groupedMean = timeSeries.groupby(by=multiIndex).mean()

    # change the Index (made of YYYY,mm,dd) to pd.DatetimeIndex
    groupedMean.index = pd.to_datetime(groupedMean.index.map(lambda x: '-'.join(map(str, x)))) 

    # add 12h offset
    groupedMean.index   += pd.Timedelta("12h")

    """
        # create dates to re-index (full dates) 
        fullDateRange   = pd.date_range(start=pd.Timestamp(groupedMean.index.year.min(),1,1),
                                    end=pd.Timestamp(groupedMean.index.year.max(), 12, 31), 
                                    freq='D')
    
        fullDateRange   += pd.Timedelta("12h")
    """

    """
    if reindex:
        if clipDates == None:
            raise Exception("Error: you must set clipDates to reindex!")

        fullDateStart   = parse(clipDates[0])                                   
        fullDateEnd     = parse(clipDates[1])     

        print("     reindex dbg: ")
        print(groupedMean.index[0])
        print(groupedMean.index[-1])
        print(fullDateStart)
        print(fullDateEnd)

        fullDateRange   = pd.date_range(start=fullDateStart,end=fullDateEnd,freq="D") 
        fullDateRange   += pd.Timedelta("12h")

        if verbose:
            print("   Reindexing to : ")
            print("           start : %s" % fullDateStart)
            print("             end : %s" % fullDateEnd)
            groupedMean     = groupedMean.reindex(fullDateRange)
    """
            
    # remove anomaly
    if anomaly:
        groupedMean.values[:] = groupedMean.values[:] - groupedMean.mean()

    # create output DataArray
    daOut = xr.DataArray(groupedMean.values, 
                        dims = ["time_counter"], 
                        coords = {"time_counter" : groupedMean.index, 
                                  "nav_lon": da["nav_lon"].values.item(), 
                                  "nav_lat" : da["nav_lat"].values.item() 
                                  } 
                        )

    daOut.name = "ssh"     

    # return a xr.DataArray
    return daOut

def pearsonCorr(da1=xr.DataArray,da2=xr.DataArray,verbose=False):
    # calculate Pearson correlation coefficient between X and Y 
    # DataArrays

    # r = A / B

    # A =  sum_i ( (X_i - <X>) (Y_i - <Y>) ) = cov(X,Y) 
    # B1 = sum_i ( (X_i - <X>)**2 ) = variance(x)
    # B2 = sum_i ( (Y_i - <Y>)**2 ) = variance(Y)
    # B  = sqrt  ( B1 * B2)

    # check on size
    if da1.size != da2.size:
        raise Exception("DataArray sizes do not match! size1: %s size2: %s" % (da1.size,da2.size))

    if verbose:
        print("shape da1:", da1.shape)
        print("shape da2:", da2.shape)

    # calculate averages of DataArrays
    values1 = da1.values
    values2 = da2.values

    deviations1 = values1 - np.nanmean(values1)
    deviations2 = values2 - np.nanmean(values2)

    # calc numerator ( Covariance)
    A = np.nansum(np.multiply ( deviations1, deviations2 ) )

    # calc denominator factors ( Variances)
    B1 = np.nansum ( deviations1 ** 2)
    B2 = np.nansum ( deviations2 ** 2)

    B = np.sqrt( B1 * B2 )

    rCoeff = A / B

    return rCoeff,None

def pearsonCorr2(da1=xr.DataArray,da2=xr.DataArray,verbose=False):
    # wrapper to scipy.stats.pearsonr

    # check on size
    if da1.size != da2.size:
        raise Exception("DataArray sizes do not match! size1: %s size2: %s" % (da1.size,da2.size))

    if verbose:
        print("    pearsonr. shape da1:", da1.shape)
        print("     pearsonr.shape da2:", da2.shape)

    # check for Nan
    nas = np.logical_or(np.isnan(da1.values), np.isnan(da2.values))
    if verbose:
        print("nas")
        print(nas)
    result  = pearsonr(da1[~nas], da2[~nas])
    rCoeff  = result.statistic
    pValue  = result.pvalue

    return rCoeff,pValue

def calcNanRatio(daIN=xr.DataArray,verbose=False):
    # calculate the percentage of Nan with respect to the
    # length of the timeseries
    if verbose:
        print("calcNanRatio: size of input dA: %d" % daIN.size)

    nNan = np.sum(np.isnan(daIN))
    if verbose:
        print("calcNanRatio: nNan: %d " % nNan)
    nTot = daIN.size 
    if verbose:
        print("calcNanRatio: nTot: %d " % nTot)

    nanRatio = nNan / nTot
    if verbose:
        print("calcNanRatio: nanRatio: %.2f" % nanRatio.values.item())
    nanRatio = nanRatio.values.item()

    return nanRatio
    
def calcShortestPath(origin = [0,0],destination = [0,0],unit="km",method="scgraph",verbose=False):
    # use one of searoute / scgraph packages to 
    # compute an approximate shortest route between 2 points
    # and return the distance associated in km

    # note!!
    # origin        = [lon,lat]
    # destination   = [lon,lat]

    if method not in ["searoute","scgraph"]:
        raise Exception("Error (calcShortestPath): method not recognized")

    if verbose:
        print("calcShortestPath: calculating distance between:")
        print("         origin(lon,lat): (%.2f,%.2f)" % (origin[0],origin[1]))
        print("    destination(lon,lat): (%.2f,%.2f)" % (destination[0],destination[1]))

    dist = 0
    if method == "scgraph":
        outPut = marnet_geograph.get_shortest_path(
            origin_node={"latitude": origin[1], "longitude":origin[0] },
            destination_node = {"latitude": destination[1] ,"longitude": destination[0] } 
                                        )
        dist = outPut["length"]
        if verbose:
            print("   Routepath Coordinates:")
            print(outPut["coordinate_path"])
    if method == "searoute":
        outPut = sr.searoute(origin,
                             destination,
                             units=unit,
                             append_orig_dest=True)
        dist = outPut["properties"]["length"]
        if verbose:
            print("   Routepath Coordinates:")
            print(outPut["geometry"]["coordinates"])    

    if verbose:
        print("     total distance[km]: %.2f" % dist)

    return dist