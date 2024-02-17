import xarray as xr 
import matplotlib.pyplot as plt 
import numpy as np 

path = "../cdo_create_LBC/test_out_cons2.nc"

try:
    ds = xr.open_dataset(path)
except:
    ds = xr.open_dataset(path,decode_times=False)    

ds = ds.isel(T=0,Z=0)

NxDomain = ds.dims['X']
NyDomain = ds.dims['Y']

# T - points
nav_lon = ds["nav_lon"]
nav_lat = ds["nav_lat"]

print(nav_lon)

plt.figure(figsize=(9,14))
cf = plt.pcolormesh(nav_lon,nav_lat,nav_lon,edgecolors="k")
# add text

for j in range(10):
    for i in range(10):
        x = nav_lon.isel(X=i,Y=j).data
        y = nav_lat.isel(X=i,Y=j).data
        plt.text(x,y,"%d,%d" % (i+1,j+1),fontsize="xx-small",color="r")


cb = plt.colorbar(cf)
plt.xlim(13.15,13.6)
plt.ylim(54.35,54.85)
cb.set_label("nav_lon")
figname = "coordinates_ordering"
plt.savefig(figname,bbox_inches="tight",dpi=600)

plt.show()




