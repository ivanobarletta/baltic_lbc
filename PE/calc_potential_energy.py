import xarray as xr 
import sys
from os.path import basename
import numpy as np
import pandas as pd 
from load_mesh_file import load_nemo_mesh_file

sys.path.insert(1,"/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/xGCM")

from xgcm_utils import create_nemo_file_list

### Still under development!!!

"""
    I compute the potential energy of a series of density NEMO files

    the principle of computation is simple

    Epot = -\int g*rho*(z-z_ref)*dV [J]
                [m/s**2] *[kg/m**3] * [m] * [m**3] = [J]   

    the minus sign is because the depth is positive downwards
    
    a) rho is taken from the files in the list
    b) z is gdept_0 from the mesh file. The dimensions (z,y,x)
        of gdept_0 must be the same of the input density files
    c) z_ref is a value contained in a dataset provided in pathReferenceDepth. The value
        must refer on the same geometry of the input density files.        
    d) dV(Z,Y,X) are the volumes of the cells, computed from the mesh
        file in pathMesh

"""

date1 = "20211229"
date2 = "20231230"
#date2 = "20220105"
"""

# name of experiment
experimentName = "CTR-Run"

# path of mesh static file
pathMesh = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/"+\
            "CONTROL_run_outs/ZNB_native_mesh/static/meshmask_SIREN/"

# path of density files
rootPath = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/"+\
            "CONTROL_run_outs/ZNB_native_mesh/density"

fileNameRoot = "NEATL36_1d25h-m"            

"""
# name of experiment
experimentName = "BAL-Run"

# path of mesh static file
pathMesh = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/"+\
            "TESTRUN2_run_outs/ZNB_native_mesh/static/meshmask_SIREN/"

# path of density files
rootPath = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/"+\
            "TESTRUN2_run_outs/ZNB_native_mesh/3DT-density"
            
fileNameRoot = "NEATL36_TESTRUN_1d25h-m"

pathReferenceDepth = ""

# get mesh information
print("processing NEMO mesh file")
dsMesh = load_nemo_mesh_file(pathMesh=pathMesh,in_msh=3)

###############################################################
print("computing volumes")
# calculate volumes
e1t = dsMesh["e1t"]     #(Y,X)
e2t = dsMesh["e2t"]     #(Y,X)
e3t = dsMesh["e3t_0"]   #(Z,Y,X)
tmask = dsMesh["tmask"] #(Z,Y,X)
gdept_0 = dsMesh["gdept_0"] #(Z,Y,X)

# assign coordinates to all arrays
xcoord = np.arange(e1t.sizes["X"]).astype(np.int32)+1
ycoord = np.arange(e1t.sizes["Y"]).astype(np.int32)+1
zcoord = np.arange(e3t.sizes["Z"]).astype(np.int32)+1

e1t = e1t.assign_coords({"X":xcoord,"Y":ycoord,
                    "glamt":dsMesh["glamt"],"gphit":dsMesh["gphit"]})

e2t = e2t.assign_coords({"X":xcoord,"Y":ycoord,
                    "glamt":dsMesh["glamt"],"gphit":dsMesh["gphit"]})

e3t = e3t.assign_coords({"X":xcoord,"Y":ycoord,"Z":zcoord,
                    "glamt":dsMesh["glamt"],"gphit":dsMesh["gphit"],"gdept_1d":dsMesh["gdept_1d"]})

gdept_0 = gdept_0.assign_coords({"X":xcoord,"Y":ycoord,"Z":zcoord,
                    "glamt":dsMesh["glamt"],"gphit":dsMesh["gphit"],"gdept_1d":dsMesh["gdept_1d"]})

tmask   = tmask.assign_coords({"X":xcoord,"Y":ycoord,"Z":zcoord,
                    "glamt":dsMesh["glamt"],"gphit":dsMesh["gphit"],"gdept_1d":dsMesh["gdept_1d"]})

# compute volume
volT = e3t*e1t*e2t  # I compute in this order to get (Z,Y,X)
# mask land values
volT = volT.where(tmask==1)

volT.name = "volumeT_0"
volT = volT.assign_attrs({"units":"m^3"})
################################################################################


# create list of files
print("creating list of density files")
densityFileList = create_nemo_file_list(rootPath=rootPath
                                        ,fileNameRoot=fileNameRoot,
                                        fileType="3DT",
                                        varName="density",
                                        date1=date1,
                                        date2=date2,verbose=True)

zref = 15.0 # arbitrary zref
grav = 9.81 # m/s**2

dateList    = pd.date_range(date1,date2,freq="d").strftime("%Y%m%d")
valueList   = []

for file,date in zip(densityFileList,dateList):
    print("file: %s" % basename(file))
    rho = xr.open_dataset(file)["sigma"]+1000 
    integral = rho.data * (gdept_0 - zref).data * volT.data 
    valueList.append(-grav*np.nansum(integral))
    

daOut = xr.DataArray(data = np.array(valueList),dims=["time_counter"],
                     coords={"time_counter":pd.to_datetime(dateList)})

daOut.name = "PE"
daOut = daOut.assign_attrs({"units":"J","Experiment":experimentName,
                                "long_name":"Pontential Energy","zref":zref})

fileNameOut = "NEATL36_PE_%s_%s_%s.nc_ZNB" % (experimentName,date1,date2)

daOut.to_dataset().to_netcdf(fileNameOut,encoding={"time_counter" : {"dtype": "i4"} })
