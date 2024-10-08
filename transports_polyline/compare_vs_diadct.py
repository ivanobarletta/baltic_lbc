import matplotlib.pyplot as plt
import xarray as xr
from glob import glob
import pandas as pd
from os.path import isfile

def resample2Day(ds):
    ds2 = ds.resample(index="D").mean()
    indexx = ds2["index"]
    indexx = indexx + pd.Timedelta("12h")
    ds2["index"] = indexx
    return ds2

transectsList   = ["EAST2","ARKONA","DARSS","GREATBELT","SOUND","KATTEGAT","SKAGERRAK"]

pathDiadct      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/TESTRUN2/transports/diadct"
pathPolyline    = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/transports/polyline"

transportType   = "heat"

fileDiadctList  = []
for transect in transectsList:
    print(transect)
    print(glob(pathDiadct + "/" + "*" + transportType + "*" + transect +  "*.nc" ))
    filePath = glob(pathDiadct + "/" + "*" + transportType + "*" + transect +  "*.nc" )[0]
    fileDiadctList.append(filePath)

filePolylineList    = []
for transect in transectsList:
    print(transect)
    print(glob(pathPolyline + "/" + "*" +"test4" + "*" + transect +  "*.nc" ))
    filePath = glob(pathPolyline + "/" + "*" +"test4" + "*" + transect +  "*.nc" )[0]
    filePolylineList.append(filePath)

# check existence of files
for filePath in fileDiadctList:
    print("existence of: %s" % filePath)
    print(isfile(filePath))

for filePath in filePolylineList:
    print("existence of: %s" % filePath)
    print(isfile(filePath))


# check length of list (must have same length)
nFilesDiadct    = len(fileDiadctList)
nFilesPolyline  = len(filePolylineList)

if nFilesDiadct != nFilesPolyline:
    raise Exception("Error: the lists have not the same length! Check the paths")

for iTransect,transect in enumerate(transectsList):
    dsDiadct    = xr.open_dataset(fileDiadctList[iTransect])
    resampled   = resample2Day(dsDiadct)
    dsPolyline  = xr.open_dataset(filePolylineList[iTransect])

    plt.figure()

    resampled["transport_total"].plot(label = "diadct resampled",linewidth =0.5, color="gray")
    varName = "%s_transport_total" % transportType
    dsPolyline[varName].plot(label="outputs polyline",linewidth=0.5,color="k")
    plt.legend()
    plt.title("%s Transport at %s " % (transportType,transect))
    plt.xlabel("")

    dsDiadct.close()
    dsPolyline.close()
