import xarray as xr 
import matplotlib.pyplot as plt 
import numpy as np 


neatlBathyPath  = "bathy_meter2.nc"
balticBathyPath = "bathy_BALTIC_ger_swe.nc"
neatlEast2LBC   = "../neatl_lbc_files/NEATL36_obcdta_east_2_20220115P01_R20220123.nc"

dsNeatl     = xr.open_dataset(neatlBathyPath)
dsBaltic    = xr.open_dataset(balticBathyPath)

dsNeatl     = dsNeatl.isel(x=-2,y=slice(1475,1548))

print(dsNeatl)

bathyNeatl  = dsNeatl["Bathymetry"]
bathyBaltic = dsBaltic["deptho"]

nav_lon     = dsNeatl["nav_lon"]
nav_lat     = dsNeatl["nav_lat"]

x_int = xr.DataArray(nav_lon.data,dims="new_coord")
y_int = xr.DataArray(nav_lat.data,dims="new_coord") 

bathyBalticTransect = bathyBaltic.interp(longitude=x_int,latitude=y_int)

fig = plt.figure(figsize=(8,5))
plt.scatter(y_int,bathyNeatl.data,s=1,color="b")
plt.plot(y_int,bathyNeatl.data,label="neatl36",color="b")
plt.scatter(y_int,bathyBalticTransect.data,s=1,color="k")
plt.plot(y_int,bathyBalticTransect.data,label="cmems-baltic",color="k")
plt.gca().invert_yaxis()
plt.ylim(50,0)
plt.legend()
plt.ylabel("Depth [m]")
plt.xlabel("Latitude along Transect")

figname = "bathy_transect"
plt.savefig(figname,bbox_inches="tight",dpi=600)

plt.show()



