import xarray as xr 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

nemo_lbc_path = "NEATL36_obcdta_east_2_20220515P01_R20220524.nc"

ds = xr.open_dataset(nemo_lbc_path)

lons = ds.nav_lon
lats = ds.nav_lat

proj = ccrs.PlateCarree()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([12, 15, 52, 56], crs=proj)

ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE)
#ax.add_feature(cfeature.BORDERS, linestyle=':')
#ax.add_feature(cfeature.LAKES, alpha=0.5)
#ax.add_feature(cfeature.RIVERS)

ax.scatter(lons,lats,s=0.1)
ax.scatter(lons.isel(X=0),lats.isel(X=0),s=0.1,color="r")

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl.xlabels_top = False
gl.ylabels_left = False
#ax.coastlines(resolution="10m")

plt.title("East2 Boundary Points")

plt.savefig("nemo_east2_lbc_points",bbox_inches="tight",dpi=600)

#plt.show()


