import xarray as xr
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
from os.path import isfile,join
import sys
from glob import glob
from datetime import datetime
from processTimeSeriesTools import *
import cartopy.crs as ccrs

# I process the TG and corresponding Model outputs
# I create 2 DataArrays that gather the data to be 
# processes with EOF analysis (not in this script) 

logFileName = "log_%s" % datetime.now().strftime("%Y%m%d%H%M")

obc_coords  = [13.6,55]     # calculate distance of TGs from this point 

# temporal limit for the analysis
date1       = "2022-1-1"
date2       = "2022-4-1"
date2       = "2023-7-1"

# the inventory contains the list of TG to process
inventoryPath   = "inventoryTG_process.csv"
#inventoryPath   = "ciao.csv"

minDataRatio    = 0.9       # lower limit for the ratio (nNan) / (nTot) in a timeseries
pValueMax       = 0.01

if isfile(inventoryPath) == False:
    raise Exception("inventory file not existing")

# open inventory with tidegauges filenames
df      = pd.read_csv(inventoryPath,header=0,sep=",")

print("Processing TG in inventory: %s" % inventoryPath)
print()
print(df)
print()

# root folder for tidegauges files
rootTG          = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TG/in_situ_historic"
tgFilesList     = []

# root folder of model timeseries
rootModelSSH    = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE"

experimentNames = ["CONTROL","EXP_BT1","EXP_BT2"]
experimentColors    = ["b","c","r"]

experimentNames = ["CONTROL","EXP_BT1","EXP_BT2","TESTRUN"]
experimentColors    = ["b","c","r","g"]

experimentNames = ["CONTROL","TESTRUN"]
experimentColors    = ["b","g"]


if len(experimentNames) != len(experimentColors):
    raise Exception("Error: Experiment names and Colors must have same length")

nExperiments        = len(experimentNames)
nTG                 = df.shape[0]
nTGOK               = 0 

tgOK                = []            # list  
tgDistList          = []            # list containing distances of TGs from obc         
nTGOK               = 0 

# gather coords and names of TGOK
lonsOK              = []
latsOK              = []
tgNamesOK           = []

# estimate length of timeseries
nT                  = int((parse(date2)-parse(date1)) / pd.Timedelta("24h")) + 1

# create temporary arrays to store data
tempArrayTG         = np.zeros((nT,nTG))
tempArrayModel      = np.zeros((nExperiments,nT,nTG))

# loop through the tidegauges files
for itg,(tgName,tgLon,tgLat) in enumerate(zip(df["#ID"],df["lon"],df["lat"])):
    print("----------------------")
    print("TG Name: %s" % tgName)
    tgFilePath = glob(join(rootTG,"*"+tgName+"*.nc"))[0]
    if isfile(tgFilePath) == False:
        raise Exception("TG file not existing")

    tgFilesList.append(tgFilePath)

    tgDist = calcShortestPath( origin = obc_coords, destination=[tgLon,tgLat],verbose=True )

    # open TG Dataset      
    dsTG    = xr.open_dataset(tgFilePath)
    # resample and return DataArray
    daTG    = resampleCmemsTG(dsIN=dsTG,clipDates=(date1,date2),verbose=True)
    # calculate Daily average
    daTG    = DailyAverageCmemsTG(daIN=daTG,nMin=80,clipDates=(date1,date2),verbose=True,reindex=True)
    
    nanRatio = calcNanRatio(daIN=daTG,verbose=True)
    print("nanRatio: %.2f" % nanRatio )
    goodDataRatio = 1. - nanRatio

    tgDistList.append(tgDist)

    if goodDataRatio < minDataRatio:
        print("This TG has too few data [less than %d %%] " % (100*minDataRatio) )
        print("  ...skipping...  ")
        tgOK.append("False")
        continue

    # calc some stats about the timeseries        
    stats = calcTGStats(daIN=daTG)    

    if stats["lConsecutiveNans"] == True:
        print("The timeseries has consecutive Nans. Nan will not be filled")
        print("skipping this TG...")
        tgOK.append("False")
        continue

    tgOK.append("True")    
    lonsOK.append(tgLon)
    latsOK.append(tgLat)
    tgNamesOK.append(tgName)

    nTGOK += 1

    # fill Nan by interpolating
    daTG = interpolateTG(daIN=daTG,verbose=True)

    # fill temporary array for TG
    tempArrayTG[:,nTGOK-1] = daTG.values 

    # search for corresponding model timeseries
    for iExp,(experiment,color) in enumerate(zip(experimentNames,experimentColors)):
        print("   Experiment: %s" % experiment)
        modelPath = glob(join(rootModelSSH,experiment,"timeseries","*%s*" % tgName))[0]
        if isfile(modelPath) == False:
            raise Exception("Model timeseries %s not existing" % modelPath)
        print("   Model Path: %s" % modelPath)

        # open model ssh DA
        daModelSSH  = xr.open_dataset(modelPath)["ssh"]
        # calc daily Average
        daModelSSH  = DailyAverageNemoSSH(daIN=daModelSSH,clipDates=(date1,date2),verbose=True,reindex=True)
        # plot
        #daModelSSH.plot(ax=ax,color=color,linewidth=0.6,label=experiment)     

        tempArrayModel[iExp,:,nTGOK-1] = daModelSSH.values

    """
    ax.set_xlim(datetime(2022,1,1),datetime(2022,4,1))
    ax.set_ylim(-1.5,1.5)
    ax.set_title(tgName)
    ax.legend()
    plt.savefig("%s/tg_plots_%s_" % (outFolder,tgName),dpi=600,bbox_inches="tight")
    plt.close()
    """

print("TG remaining after process: %d" % nTGOK )

times = pd.date_range(start=date1,end=date2,freq="1d") + pd.Timedelta("12h")

# build output DataArray for TG
daOutTG = xr.DataArray(
                        data = tempArrayTG[:,:nTGOK],               # take only the part of array filled with TGOK
                        dims = ["time","nTG"],
                        coords = dict(
                                    time = times,
                                    lon = (["nTG"], lonsOK),
                                    lat = (["nTG"], latsOK),
                                    ),
                        attrs=dict(
                            description="ssh from daily averaged TG",
                            units="m"
                                    ),
                        )

daOutTG.name = "SLEV"

# build output DataArray for Models
daOutModel = xr.DataArray(
                        data = tempArrayModel[:,:,:nTGOK],          # take only the part of array filled with TGOK
                        dims = ["nExp","time","nTG"],
                        coords = dict(
                                    time = times,
                                    lon = (["nTG"], lonsOK),
                                    lat = (["nTG"], latsOK),
                                    ),
                        attrs=dict(
                            description="ssh from daily averaged model outputs",
                            units="m",
                                    ),
                        )

daOutModel.name = "ssh"

# create a Datarray with Names
daOutNames  = xr.DataArray(
                        data = np.asarray(tgNamesOK),
                        dims = ["nTG"], 
                        attrs=dict(
                            description="names of TideGauges Locations",    
                                )
                            )

daOutNames.name = "TGNames"

# create a Datarray with Experiment names
daOutExpNames  = xr.DataArray(
                        data = np.asarray(experimentNames),
                        dims = ["nExp"], 
                        attrs=dict(
                            description="names of experiments",            
                                )
                            )

daOutExpNames.name = "ExpNames"

# create a Datarray with Experiment colors
daOutExpColors  = xr.DataArray(
                        data = np.asarray(experimentColors),
                        dims = ["nExp"], 
                            )
daOutExpColors.name = "ExpNames"

dsOut       = daOutTG.to_dataset()
dsOut["ssh"]        = daOutModel
# integrate names into datasets
dsOut["TGNames"]    = daOutNames
dsOut["ExpNames"]   = daOutExpNames
dsOut["ExpColors"]  = daOutExpColors

str1    = parse(date1).strftime("%Y%m%d")
str2    = parse(date2).strftime("%Y%m%d")
# save to netcdf
dsOut.to_netcdf("dsOut_4EOF_%s_%s.nc" % (str1,str2), encoding={"time" : {"dtype": "i4"} } )
