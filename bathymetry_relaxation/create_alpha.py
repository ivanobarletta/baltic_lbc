import xarray as xr 
import matplotlib.pyplot as plt 
import numpy as np 

neatlBathyPath  = "bathy_meter2.nc"
outDataset      = "alpha.nc"

# open dataset and get sizes
dsNeatl = xr.open_dataset(neatlBathyPath)
print(dsNeatl)
nx = dsNeatl.sizes["x"]
ny = dsNeatl.sizes["y"]
nav_lon = dsNeatl["nav_lon"]
nav_lat = dsNeatl["nav_lat"]

xidx0,xidx1 = 1078-1,1092-1
yidx0,yidx1 = 1476-1,1548-1

# build integer coords
ii = np.arange(0,nx,1)
jj = np.arange(0,ny,1)

ii,jj = np.meshgrid(ii,jj)

m = -1.0 / (xidx1-xidx0)

q = -m * xidx1

# clip function values
f = m * ii + q
f = np.clip(f,a_min=0,a_max=1)
f[jj<yidx0] = 1
f[jj>yidx1] = 1

"""
plt.figure()
cf = plt.pcolormesh(ii,jj,f)
plt.scatter(ii,jj,s=0.5,color="k")
plt.colorbar(cf)
plt.show()
"""


da = xr.zeros_like(dsNeatl["Bathymetry"])
da = da.rename("alpha1")
da.data = f
da.to_netcdf(outDataset)





