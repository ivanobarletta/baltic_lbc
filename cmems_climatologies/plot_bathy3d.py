import xarray as xr
import matplotlib.pyplot as plt
import cmocean
import os

path = "../domain/bathy_meter.nc"
outFolder = "bathymetry"

dsbathy = xr.open_dataset(path)

lons    = dsbathy["nav_lon"]
lats    = dsbathy["nav_lat"]
bathy   = dsbathy["Bathymetry"]

ax = plt.figure().add_subplot(projection="3d"); ax.plot_surface(lons,lats,bathy,cmap=cmocean.cm.deep);ax.invert_zaxis();ax.grid(visible=None)

figname = "neatl36_bathymetry_3d" 

plt.savefig(os.path.join(outFolder,figname),bbox_inches="tight",dpi=300)