import matplotlib.pyplot as plt
import numpy as np 
import xarray as xr 
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

dsAlpha    = xr.open_dataset("alpha.nc")

cmapp       = "YlGn_r"

east2Box    = [1078,1092,1476,1548]

# load coordinates
nav_lon,nav_lat = dsAlpha["nav_lon"],dsAlpha["nav_lat"]
# load alpha field
alpha           = dsAlpha["alpha1"]


# step to subsample the coordinates
subsample = 10

# smaller coordinates for smaller zoom axes (takes a lot to make the plot)
nav_lon2    = nav_lon.isel(x=slice(-100,None),y=slice(1400,1600))
nav_lat2    = nav_lat.isel(x=slice(-100,None),y=slice(1400,1600))
alpha2      = alpha.isel(x=slice(-100,None),y=slice(1400,1600))


lonsT = nav_lon.isel(x=slice(0,None,subsample),y=slice(0,None,subsample))
latsT = nav_lat.isel(x=slice(0,None,subsample),y=slice(0,None,subsample))

lons_east2  = nav_lon.isel(x=slice(east2Box[0],east2Box[1]),y=slice(east2Box[2],east2Box[3]))
lats_east2  = nav_lat.isel(x=slice(east2Box[0],east2Box[1]),y=slice(east2Box[2],east2Box[3]))

def get_poly_corners(lons,lats):
    lons = lons.data
    lats = lats.data
    idx = ((0,0),(0,-1),(-1,-1),(-1,0),(0,0))
    x = np.asarray([lons[i,j] for i,j in idx])
    y = np.asarray([lats[i,j] for i,j in idx])

    return x,y

xx,yy = get_poly_corners(lons_east2,lats_east2)

# projection of map
projj = ccrs.EuroPP()
# projection of data
data_crs = ccrs.PlateCarree()

#Figure 1
#plt.figure()
#ax = plt.axes(projection=projj)

fig, ax = plt.subplots(1,1, subplot_kw={"projection":projj},figsize=(7,12))

#ax.plot(lonsT  ,latsT  ,color="k",linewidth=0.1,transform=data_crs)
#ax.plot(lonsT.T,latsT.T,color="k",linewidth=0.1,transform=data_crs)
#ax.plot(xx,yy,color="r",linewidth=0.5,transform=data_crs)
cf = ax.pcolormesh(nav_lon,nav_lat,alpha,cmap=cmapp,transform=data_crs)

ax.coastlines(resolution='10m')

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl.xlabels_top = False
gl.ylabels_left = False

cb = fig.colorbar(cf,ax=ax,orientation="horizontal",pad=0.05)
cb.set_label(r"$\alpha$",fontsize=15)

# create smaller axes inside
#ax_inset = inset_axes(ax,width="30%",height="30%",loc="upper left",axes_class=cartopy.mpl.geoaxes.GeoAxes,axes_kwargs=dict(projection=projj))

ax_inset = fig.add_axes([0.17, 0.5, 0.3, 0.3], projection=projj)
ax_inset.set_extent([12,14,54,56])

#ax_inset.plot(nav_lon2  ,nav_lat2  ,color="k",linewidth=0.1,transform=data_crs)
#ax_inset.plot(nav_lon2.T,nav_lat2.T,color="k",linewidth=0.1,transform=data_crs)
#ax_inset.plot(xx,yy,color="r",linewidth=0.5,transform=data_crs)

ax_inset.pcolormesh(nav_lon2,nav_lat2,alpha2,cmap=cmapp,transform=data_crs)

ax_inset.coastlines(resolution='10m')

gl2 = ax_inset.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
    linewidth=1, color='gray', alpha=0.5, linestyle='--')

gl2.xlabels_top     = False
gl2.ylabels_left    = False
gl2.xlabels_bottom  = False
gl2.ylabels_right   = False

figname = "alpha_neatl36_and_east2"



plt.savefig(figname,bbox_inches="tight",dpi=600)
