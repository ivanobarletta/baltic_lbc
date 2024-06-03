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

logFileName = "log_%s" % datetime.now().strftime("%Y%m%d%H%M")

obc_coords  = [13.6,55]     # calculate distance of TGs from this point 

# the inventory contains the list of TG to process
inventoryPath   = "inventoryTG_plot_all.csv"
#inventoryPath   = "ciao.csv"

outFolder       = "plots_daily_ave"

dfOutName       = "tg_rPearson_all_4expetiments.csv"

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

if len(experimentNames) != len(experimentColors):
    raise Exception("Error: Experiment names and Colors must have same length")

nExperiments        = len(experimentNames)
nTG                 = df.shape[0]
nTGOK               = 0 

rPearsonArray       = np.zeros((nTG,nExperiments))
pValueArray         = np.zeros((nTG,nExperiments))
tgOK                = []            # list  
tgDistList          = []            # list containing distances of TGs from obc         

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
    daTG    = resampleCmemsTG(dsIN=dsTG,clipDates=("2022-1-1","2022-4-1"),verbose=True)
    # calculate Daily average
    daTG    = DailyAverageCmemsTG(daIN=daTG,nMin=80,clipDates=("2022-1-1","2022-4-1"),verbose=True,reindex=True)
    
    nanRatio = calcNanRatio(daIN=daTG,verbose=True)
    print("nanRatio: %.2f" % nanRatio )
    goodDataRatio = 1. - nanRatio

    tgDistList.append(tgDist)

    if goodDataRatio < minDataRatio:
        print("This TG has too few data [less than %d %%] " % (100*minDataRatio) )
        print("  ...skipping...  ")
        tgOK.append("False")
        continue

    tgOK.append("True")    

    # start a figure for each TG
    fig,ax = plt.subplots(figsize=(10,4))

    daTG.plot(ax=ax,color="k",linewidth=0.5)

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
        daModelSSH  = DailyAverageNemoSSH(daIN=daModelSSH,clipDates=("2022-1-1","2022-4-1"),verbose=True,reindex=True)
        # plot
        daModelSSH.plot(ax=ax,color=color,linewidth=0.6,label=experiment)     

        # calculate Pearson Correlation coefficient
        #rPearson,pVal = pearsonCorr(da1=daTG,da2=daModelSSH,verbose=True)
        rPearson,pVal = pearsonCorr2(da1=daTG,da2=daModelSSH,verbose=True)

        print("     rPearson: %.2f" % rPearson)
        print("       pValue: %.5e" % pVal)

        # fill rPearson array 
        rPearsonArray[itg,iExp] = rPearson
        pValueArray[itg,iExp] = pVal    

    ax.set_xlim(datetime(2022,1,1),datetime(2022,4,1))
    ax.set_ylim(-1.5,1.5)
    ax.set_title(tgName)
    ax.legend()
    plt.savefig("%s/tg_plots_%s_" % (outFolder,tgName),dpi=600,bbox_inches="tight")
    plt.close()


print(rPearsonArray)
# save a new DataFrame with rPearson coefficients

# initialize output DataFrame
dfOut = pd.DataFrame()
dfColumns = ["region","#ID","lon","lat"]

# populate new df with selected columns
for col in dfColumns:
    dfOut[col] = df[col]

for i,expName in enumerate(experimentNames):
    newColName  = "rP_%s" % expName
    newColName2 = "pV_%s" % expName
    #print(i,expName,newColName)
    dfOut[newColName] = rPearsonArray[:,i]
    dfOut[newColName2] = pValueArray[:,i]

# mark ruled out TG with "false" 
dfOut["tgOK"]   = np.asarray(tgOK)
dfOut["tgDist"] = np.asarray(tgDistList)

# this sets all float less than 0.001 to 0
# the pValues that I find are generally really small (1e-10 in the 
# biggest case )
# I retain the columns of pValues, but they are not really necesary
# correlation is always significant

dfOut.to_csv(dfOutName,float_format="%.3f")


#plt.ion()
#plt.show()




