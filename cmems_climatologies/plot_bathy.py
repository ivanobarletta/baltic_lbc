import xarray as xr
import matplotlib.pyplot as plt
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os

path = "../domain/bathy_meter.nc"
outFolder = "bathymetry"

dsbathy = xr.open_dataset(path)

lons    = dsbathy["nav_lon"]
lats    = dsbathy["nav_lat"]
bathy   = dsbathy["Bathymetry"]

# projection of map
#projj = ccrs.EuroPP()
projj = ccrs.NearsidePerspective(central_longitude=0,central_latitude=45,satellite_height=24e6)
# projection of data
data_crs = ccrs.PlateCarree()

fig, ax = plt.subplots(1,1, subplot_kw={"projection":projj},figsize=(7,12))

ax.pcolormesh(lons,lats,bathy,cmap=cmocean.cm.deep,transform=data_crs)

ax.coastlines(resolution='10m')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl.xlabels_top  = False
gl.ylabels_left = True
gl.ylabels_right = False

figname = "neatl36_bathymetry_2d" 

plt.savefig(os.path.join(outFolder,figname),bbox_inches="tight",dpi=300)

#plt.show()