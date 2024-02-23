import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import cmocean
import os

"""
	Note!!!
	changing the projection can be quite detrimental for this script!!!

	I suggest not to change the projection (leave PlateCarree). When I tried
	to set EuroPP I got a lot of errors from quiver function. Possibly there 
	are problem with vector transformation.

"""

outFolder = "maps"

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

depth_level	= 0

ds_glo = xr.open_dataset("../cmems_glob12/glo_anfc_subset.nc").isel(depth=depth_level)
ds_bal = xr.open_dataset("../cmems_baltic/bal_anfc_subset2.nc").isel(depth=depth_level)

u_glo = ds_glo["uo"]
v_glo = ds_glo["vo"]
u_bal = ds_bal["uo"]
v_bal = ds_bal["vo"]

u_glo_season 	= u_glo.groupby("time.season").mean(dim="time")
v_glo_season 	= v_glo.groupby("time.season").mean(dim="time")
uv_mod_glo		= np.sqrt(u_glo_season**2 + v_glo_season**2)

u_bal_season 	= u_bal.groupby("time.season").mean(dim="time")
v_bal_season 	= v_bal.groupby("time.season").mean(dim="time")
uv_mod_bal		= np.sqrt(u_bal_season**2 + v_bal_season**2)

lons_glo	= ds_glo["longitude"]
lats_glo	= ds_glo["latitude"]
lons_bal	= ds_bal["lon"]
lats_bal	= ds_bal["lat"]

lons_glo,lats_glo	= np.meshgrid(lons_glo,lats_glo)
lons_bal,lats_bal	= np.meshgrid(lons_bal,lats_bal)

proj	= ccrs.PlateCarree()

nrows = 2
ncols = 4

#fig,ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,figsize=(8,5),subplot_kw={'projection': proj})
fig,ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(17,8),subplot_kw={'projection': proj})
#fig,ax = plt.subplots(2,4,sharex=True,sharey=True,subplot_kw={'projection': ccrs.PlateCarree()},figsize=(13,8))
fig.subplots_adjust(wspace=0.1, hspace=0.02)


for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
	#ax[0,i].streamplot(lons_glo,lats_glo,u_glo_season.sel(season=season),v_glo_season.sel(season=season),transform=ccrs.PlateCarree(),linewidth=2,density=2,color=uv_mod_glo.sel(season=season))
	#ax[1,i].streamplot(lons_bal,lats_bal,u_bal_season.sel(season=season),v_bal_season.sel(season=season),transform=ccrs.PlateCarree(),linewidth=2,density=2,color=uv_mod_bal.sel(season=season))
	#ax[0,i].streamplot(lons_glo.data,lats_glo.data,u_glo_season.sel(season=season),v_glo_season.sel(season=season),transform=ccrs.PlateCarree(),linewidth=2,density=2)
	#ax[1,i].streamplot(lons_bal.data,lats_bal.data,u_bal_season.sel(season=season),v_bal_season.sel(season=season),transform=ccrs.PlateCarree(),linewidth=2,density=2)
	ax[0,i].quiver(lons_glo,lats_glo,u_glo_season.sel(season=season),v_glo_season.sel(season=season),zorder=2,width=0.005,scale=0.5)
	cf1 = ax[0,i].pcolormesh(lons_glo,lats_glo,uv_mod_glo.sel(season=season),vmin=0,vmax=0.2,cmap=cmocean.cm.thermal)
	ax[1,i].quiver(lons_bal,lats_bal,u_bal_season.sel(season=season),v_bal_season.sel(season=season),zorder=2,scale=0.8)
	cf2 = ax[1,i].pcolormesh(lons_bal,lats_bal,uv_mod_bal.sel(season=season),vmin=0,vmax=0.2,cmap=cmocean.cm.thermal)


#ax[0,1].quiver(lons_glo,lats_glo,u_glo_season.sel(season="MAM"),v_glo_season.sel(season="MAM"),transform=ccrs.PlateCarree())

for j in range(nrows):
	for i in range(ncols):
		#ax[j,i].set_facecolor('grey')

		#ax[j,i].set_ylabel("")
		#ax[j,i].set_ylabel("")

		ax[j,i].set_title("")
		ax[j,i].coastlines(resolution='10m')
		gl = ax[j,i].gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
    		linewidth=1, color='gray', alpha=0.5, linestyle='--')
		
		gl.xlabels_top = False
		gl.ylabels_left = False
		gl.xlabels_bottom = False
		gl.ylabels_right = False

		ax[j,i].set_extent([13,14,54.5,55.5])


for i,season in enumerate(("DJF", "MAM", "JJA", "SON")):
	ax[0,i].set_title(season)

cbar = fig.colorbar(cf1, ax=ax.ravel().tolist(),location="right", cmap = cmocean.cm.thermal, label = "Currents [m/s]"  , shrink = 0.8)

#ax[0,0].set_ylabel("GLO (CMEMS) 2022-2023")
#ax[1,0].set_ylabel("BAL (CMEMS) 2022-2023")
#ax[0,0].set_ylabel("ciao",color="r")

# place a text box in upper left in axes coords (I cannot make set_ylabel to work....)
ax[0,0].text(-0.2, 0.5, "GLO (CMEMS) 2022-2023", transform=ax[0,0].transAxes, fontsize=10,
        verticalalignment='center', bbox=props, rotation=90)

ax[1,0].text(-0.2, 0.5, "BAL (CMEMS) 2022-2023", transform=ax[1,0].transAxes, fontsize=10,
        verticalalignment='center', bbox=props, rotation=90)

figname = "glo_bal_vel_level_%d_climatology" % depth_level

plt.savefig(os.path.join(outFolder,figname),bbox_inches="tight")

#plt.show()
