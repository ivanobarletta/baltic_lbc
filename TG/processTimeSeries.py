import numpy as np
import xarray as xr 
import pandas as pd
from os.path import isfile
import sys
from dateutil.parser import parse

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
            print("Slicing dataset for years (%s,%s)" % (clipDates[0],clipDates[1]))
        date0 = parse(clipDates[0])
        date1 = parse(clipDates[1])
        ds = ds.sel(TIME=slice(date0,date1))

    values      = ds["SLEV"]
    timestamps  = ds["TIME"]
    qc          = ds["SLEV_QC"]

    # I take only qc=1 values
    goodValues  = qc == 1

    ds.close()

    # Create a pandas Series with the timestamps as the index
    time_series = pd.Series(values, index=pd.to_datetime(timestamps))

    # set not-good values to Nan
    time_series[goodValues.data == False] = np.nan

    # Calculate time differences between consecutive samples
    time_diffs = time_series.index.to_series().diff().dropna()

    # Round the time differences to a certain precision (e.g., seconds)
    rounded_time_diffs = time_diffs.dt.round('1s')

    # replace index of timeseries with rounded time index (it is necessary
    # otherwise some good values are set to NaN)
    time_series.index = time_series.index.round("1s")

    if rule == None:
        # Determine the most common time difference
        if verbose:
            print("Resampling using the most common time difference")    
        new_sampling_freq = rounded_time_diffs.mode().iloc[0]
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
        da = da.sel(time_counter=slice(date0,date1))

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

    # create dates to re-index (full dates) 
    fullDateRange   = pd.date_range(start=pd.Timestamp(groupedMean.index.year.min(),1,1),
                                    end=pd.Timestamp(groupedMean.index.year.max(), 12, 31), 
                                    freq='D')
    
    fullDateRange   += pd.Timedelta("12h")

    # re-index the series
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

    # create dates to re-index (full dates) 
    fullDateRange   = pd.date_range(start=pd.Timestamp(groupedMean.index.year.min(),1,1),
                                    end=pd.Timestamp(groupedMean.index.year.max(), 12, 31), 
                                    freq='D')
    
    fullDateRange   += pd.Timedelta("12h")

    # re-index the series
    groupedMean     = groupedMean.reindex(fullDateRange)

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
    # check on size
    if da1.size != da2.size:
        raise Exception("DataArray sizes do not match!")

    if verbose:
        print(da1.shape)
        print(da2.shape)

    # calculate averages of DataArrays
    values1 = da1.values
    values2 = da2.values

    deviations1 = values1 - np.nanmean(values1)
    deviations2 = values2 - np.nanmean(values2)

    # calc numerator
    A = np.nansum(np.multiply ( deviations1, deviations2 ) )

    # calc denominator factors
    B1 = np.nansum ( deviations1 ** 2)
    B2 = np.nansum ( deviations2 ** 2)

    B = np.sqrt( B1 * B2 )

    rCoeff = A / B

    return rCoeff

