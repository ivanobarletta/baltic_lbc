import pandas as pd 
import numpy as np 
import xarray as xr
from utide import solve,reconstruct 
import sys
from glob import glob
import numpy as np

if len(sys.argv) < 4:
    print("Error:")
    print("1) Provide inventory of Timeseries (.csv)")
    print("2) Provide folder where timeseries .nc are stored")  
    print("3) Provide output inventory (.csv)")    
    sys.exit(1)

inventoryPath   = sys.argv[1]
timeseriesPath  = sys.argv[2]
outFile         = sys.argv[3]

constituents    = ["M2","S2","N2"]

inventory = pd.read_csv(inventoryPath,header=0,sep=",")

# Drop columns not needed
try:
    inventory = inventory.drop(labels=["ini_year","end_year","diff_years"],axis=1)
except:
    print("No need to drop columns. I go on")


# number of constituents to analyze 
Nc      = len(constituents)
Nfiles  = inventory.shape[0]

Amplitudes  = np.zeros((Nfiles,Nc))
Phases      = np.zeros((Nfiles,Nc))


Names   = inventory["#ID"]

# build list of files ( I do in 2 steps)
# 1) 
#    I create a list of filenames with wildcards (*)...
#    
#    list of /timeseriesPath/*name*nc
#
# 2)
#    I use the glob function to find the real name 
#   
#    realName = glob( /timeseriesPath/*name*nc )

fileList2   = [timeseriesPath+"/*"+file+"*.nc" for file in Names]
fileList    = []

for file in fileList2:
    fileList.append( glob(file)[0] )

print(fileList)

# loop through the files
for i,file in enumerate(fileList):
    ds = xr.open_dataset(file)
    # make analysis
    coef = solve(ds["time_counter"].data,ds["ssh"].data.squeeze(),lat=ds["nav_lat"].data[0][0],nodal=False,trend=False,method="ols",Rayleigh_min=0.95)
    # extract harmonics (order is not always the same, they are ordered by Amplitude in descending order )
    # I have to re-calculate the indexes
    idx = []
    for const in constituents:
        idx.append(np.where(coef["name"]==const)[0])   
    # extract Amplitudes
    Ampli = coef["A"][idx].squeeze() 
    Phase = coef["g"][idx].squeeze()

    # fill arrays
    Amplitudes[i,:] = Ampli
    Phases[i,:] = Phase

# update Inventory
inventory2 = inventory.copy()

for i,const in enumerate(constituents):
    print(i,const)
    inventory2.insert(inventory2.shape[1],const+"_A",Amplitudes[:,i],True)
    inventory2.insert(inventory2.shape[1],const+"_g",Phases[:,i],True)

# save new csv
inventory2.to_csv(outFile,float_format="%.5e")
