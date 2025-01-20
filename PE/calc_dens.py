import xarray as xr 
import numpy as np
import sys 
from glob import glob 
from os.path import isfile, join, basename 
import pandas as pd 
import seawater 

sys.path.insert(1,"/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/xGCM")

from xgcm_utils import create_nemo_file_list


pathMesh3D  = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh/static/meshmask_SIREN/mesh_zgr.nc"
if isfile(pathMesh3D) == False:
    raise Exception("Error: pathMesh3D not existing. check path")    

date1       = "20211229"
date2       = "20221230"

varNameS    = "so"
varNameT    = "thetao"
fileNameRoot = "NEATL36_TESTRUN_1d25h-m"
fileType    = "3DT"

outFolder   = "density"

pathRoot0   = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/"

pathFiles   = join(pathRoot0,"TESTRUN2_run_outs/ZNB_native_mesh")

fileListS   = create_nemo_file_list(rootPath=join(pathFiles,""),fileNameRoot=fileNameRoot,
                      fileType="3DT",varName="so",date1=date1,date2=date2,verbose=True)

fileListT   = create_nemo_file_list(rootPath=join(pathFiles,""),fileNameRoot=fileNameRoot,
                      fileType="3DT",varName="thetao",date1=date1,date2=date2,verbose=True)


dateList    = pd.date_range(date1,date2,freq="d").strftime("%Y%m%d")


# access to gdept_0 (Z,Y,X) 
try:
    gdept_0     = xr.open_dataset(pathMesh3D)["gdept_0"]
except:
    print("Error! gdept_0 not found in pathMesh3D")

# loop on files. I compute density for each date

count = 0 
countMax = 800
for fileS,fileT,date in zip(fileListS,fileListT,dateList):
    count += 1
    print("date: %s" % date)
    print(basename(fileS))
    print(basename(fileT))     
    daS = xr.open_dataset(fileS)[varNameS]
    daT = xr.open_dataset(fileT)[varNameT]
    print(daS.sizes)
    print(daT.sizes)    
    # calculate sigma ( rho - 1000)
    density = seawater.dens(s=daS,t=daT,p=gdept_0) - 1000  
    # convert to single precision
    density = density.astype(np.float32)
    # create dataArray of density
    densityDA = xr.DataArray(data = density,dims=daS.dims,coords=daS.coords)
    densityDA.name = "sigma"
    densityDA = densityDA.assign_attrs({"units":"Kg/m^3","long_name":"rho-1000"})
    # build filename
    fileNameOut = "%s_3DT-density_%s-%s.nc_ZNB" % (fileNameRoot,date,date)
    densityDA.to_dataset().to_netcdf(join(outFolder,fileNameOut))
    if count > countMax:
        break









