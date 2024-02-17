import xarray as xr
import numpy as np
import matplotlib.pyplot as plt 
import sys
from lbc_utils import get_nemo_transect_coords

varname = "thetao"
varname = "so"
varname = "uo"

nemo_path = "cmems_baltic/BAL-NEMO_PHY-DailyMeans-20220115.nc"

transect_path = "neatl/NEATL36_obcdta_east_2_20220115P01_R20220123.nc"

ds = xr.open_dataset(nemo_path)

time = ds.time
lons = ds.lon
lats = ds.lat
depth = ds.depth

#build transect DataArrays
"""
x0,x1 = 13.5,14
y0,y1 = 54.4,55.4
n_points = 50
x_int = xr.DataArray(np.linspace(x0,x1,n_points),dims="new_coord")
y_int = xr.DataArray(np.linspace(y0,y1,n_points),dims="new_coord")
"""

x_int, y_int = get_nemo_transect_coords(transect_path)

print(type(time))

date_str = np.datetime_as_string(ds.time.isel(time=0),unit="h")

transect = ds[varname].isel(time=0).interp(lon=x_int,lat=y_int)

print(transect.shape)

#lats2,depth2 = np.meshgrid(lats,np.asarray(depth_list))

def contour_levels(variable):
    clevs = None
    if variable == "thetao":
        clevs = np.linspace(5,12,29)
        ccmap = "RdYlBu_r"
        clabel = "T"
    elif variable == "so":     
        clevs = np.linspace(7,15,33)
        ccmap = "RdYlBu_r"
        clabel = "PSU"    
    elif variable in ["uo","vo"]:
        clevs = np.linspace(-0.35,0.35,21)
        ccmap = "RdBu"
        clabel = "m/s"

    return clevs,ccmap,clabel  

clevs,ccmap,clabel = contour_levels(varname)

vvmin = clevs[0]
vvmax = clevs[-1]

#plt.figure()
#cf = plt.contourf(lats2,depth2,transect,levels=clevs,cmap=ccmap)
#cf = plt.pcolormesh(lats2,depth2,transect,cmap=ccmap,vmin=vvmin,vmax=vvmax)
#cb = plt.colorbar(cf)
#cb.set_label(clabel)

plt.figure()
transect.plot(cmap=ccmap,vmin=vvmin,vmax=vvmax)

plt.gca().invert_yaxis()
plt.ylim((50,0))

#plt.xlabel("Lat Along LBC")
#plt.ylabel("Depth [m]")

title = "%s - Time: %s" % (varname,date_str)

plt.title(title)
plt.savefig("cmems_baltic/%s_transect_%s" % (varname,date_str))

plt.show()

#sys.exit()

