import xarray as xr
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
from os.path import isfile,join
import sys
from glob import glob
from datetime import datetime
from processTimeSeries import *

# the inventory contains the list of TG to process
inventoryPath   = "inventoryTG_plot2.csv"

outFolder       = "plots_daily_ave"

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

# loop through the tidegauges files
for tgName in df["#ID"]:
    tgFilePath = glob(join(rootTG,"*"+tgName+"*.nc"))[0]
    if isfile(tgFilePath) == False:
        raise Exception("TG file not existing")

    tgFilesList.append(tgFilePath)

    # start a figure for each TG
    fig,ax = plt.subplots(figsize=(10,4))

    # open TG Dataset      
    dsTG    = xr.open_dataset(tgFilePath)
    # resample and return DataArray
    daTG    = resampleCmemsTG(dsIN=dsTG,clipDates=("2022-1-1","2022-4-1"),verbose=True)
    # calculate Daily average
    daTG    = DailyAverageCmemsTG(daIN=daTG,nMin=80,verbose=True)

    daTG.plot(ax=ax,color="k",linewidth=0.5)

    print("TG: %s:" % tgName)

    # search for corresponding model timeseries
    for experiment,color in zip(experimentNames,experimentColors):
        print("   Experiment: %s" % experiment)
        modelPath = glob(join(rootModelSSH,experiment,"timeseries","*%s*" % tgName))[0]
        if isfile(modelPath) == False:
            raise Exception("Model timeseries %s not existing" % modelPath)
        print("   Model Path: %s" % modelPath)

        # open model ssh DA
        daModelSSH  = xr.open_dataset(modelPath)["ssh"]
        # calc daily Average
        daModelSSH  = DailyAverageNemoSSH(daIN=daModelSSH,clipDates=("2022-1-1","2022-4-1"),verbose=True)
        # plot
        daModelSSH.plot(ax=ax,color=color,linewidth=0.6,label=experiment)     

        # calculate Pearson Correlation coefficient
        rPearson = pearsonCorr(da1=daTG,da2=daModelSSH,verbose=True)

        print("     rPearson: %.2f" % rPearson)


    ax.set_xlim(datetime(2022,1,1),datetime(2022,4,1))
    ax.set_ylim(-1.5,1.5)
    ax.set_title(tgName)
    ax.legend()
    plt.savefig("%s/tg_plots_%s_" % (outFolder,tgName),dpi=600,bbox_inches="tight")

plt.show()




