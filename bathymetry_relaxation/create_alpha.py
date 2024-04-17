import xarray as xr 
import matplotlib.pyplot as plt 
import numpy as np 

#
#   nx = 1093
#   
#   i-indexes in python notation
#   
#       
#        |                           EAST2 LBC BUFFER ZONE                         | Last Value 
#        |                                                                         | (masked)    
#   
#   1076 1077 1078 1079 1080 1081 1082 1083 1084 1085 1086 1087 1088 1089 1090 1091 1092 
#   
#          o    o    o    o    o    o    o    o    o    o    o    o    o    o    o    o   
#   
#          
#   xidx0 = 1077
#   xidx1 = 1091
#          
#   search for a linear function f[]:
#       f = m * x + q 
#       with:   
#           f[xidx0] = f[1077] = 1
#           f[xidx1] = f[1091] = 0    
#   
#   
#   slope of the linear function:
#       m = 1 / (xidx1 - xidx0 ) = 1. / 14
#   
#   y at origin:
#       f must be 0 in i=xidx1
#   
#       0 = m * xidx1 + q
#   
#       q = -m * xidx1
#   

neatlBathyPath  = "bathy_meter2.nc" # bathy_meter2.nc data is exactly as in bathy_meter.nc 
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





