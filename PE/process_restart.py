import xarray as xr 
from load_mesh_file import load_nemo_mesh_file
import numpy as np 
import seawater 



"""
    I take 1 restart file and I create a density file 
    with NEMO output standard 
        rhop(time_counter,deptht,y,x)
        * time_counter(time_counter)
        * deptht(deptht)
        nav_lon(y,x)
        nav_lat(y,x)
"""

dateRestart     = "2021-12-28T12:00:00"
pathRestart     = "restart_oce_20211228.nc_ZNB"


# take coordinates from an example density file
exampleDsPath   = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/"+ \
            "NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/CONTROL_run_outs/ZNB_native_mesh/density/"+\
            "NEATL36_1d25h-m_3DT-density_20211229-20211229.nc_ZNB"

dsExample   = xr.open_dataset(exampleDsPath) 

dsExample   = dsExample.transpose("time_counter","deptht","y","x")

# open dataset
dsRestart = xr.open_dataset(pathRestart)

# take data from rhop dataArray
rhop = dsRestart["rhop"].data - 1000

print("rhop.shape")
print(rhop.shape)

# I do this because changing the time dimension with other values is a nightmare in xarray
# note that [] is necessary for time_counter
newCoords = {"time_counter":[np.datetime64(dateRestart)],
             "deptht":dsExample["deptht"].data,
             "nav_lat":dsExample["nav_lat"],
             "nav_lon":dsExample["nav_lon"]
             }

rhopDa = xr.DataArray(data = rhop.astype(np.float32), 
                        dims=["time_counter","deptht","y","x"],
                        coords=newCoords)

print("rhopDa.shape")
print(rhopDa.shape)
print(rhopDa.sizes)

maskT = np.isnan(dsExample["sigma"])

# mask land values (strangely da.where works all way around! crazy!)
# Note! I have to use maskT.values, because if I use maskT, the difference
# in the time_counter coordinate between the 2 dataArrays raises an exeption.
rhopDa  = rhopDa.where(maskT.values==0)

rhopDa.name = "sigma"

# save to netcdf
rhopDa.to_netcdf("NEATL36_restart_density_20211228.nc_ZNB")



