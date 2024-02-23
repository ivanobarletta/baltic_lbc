import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)



depth_level	= 0

varname 	= "thetao"
cmap 		= "Spectral_r"
cmin,cmax 	= 0,20

"""
varname 	= "so"
cmap		= "hot_r"
cmin,cmax 	= 5,15
"""


ds_glo = xr.open_dataset("../cmems_glob12/glo_anfc_subset.nc")
ds_bal = xr.open_dataset("../cmems_baltic/bal_anfc_subset2.nc")

temp_glo = ds_glo[varname]
temp_bal = ds_bal[varname]

temp_glo_season = temp_glo.groupby("time.season").mean(dim="time")
temp_bal_season = temp_bal.groupby("time.season").mean(dim="time")

proj	= ccrs.EuroPP()

nrows = 2
ncols = 4

#fig,ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,figsize=(8,5),subplot_kw={'projection': proj})
fig,ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(11,7),subplot_kw={'projection': proj})
fig.subplots_adjust(wspace=0.0, hspace=0.05)


for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
	cf1 = temp_glo_season.isel(depth=depth_level).sel(season=season).plot(ax=ax[0,i],vmin=cmin,vmax=cmax,extend="both",cmap=cmap,add_colorbar=False,transform=ccrs.PlateCarree())
	cf2	= temp_bal_season.isel(depth=depth_level).sel(season=season).plot(ax=ax[1,i],vmin=cmin,vmax=cmax,extend="both",cmap=cmap,add_colorbar=False,transform=ccrs.PlateCarree())



for j in range(nrows):
	for i in range(ncols):
		ax[j,i].set_facecolor('grey')

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

		ax[j,i].set_extent([12.8,14.2,54.3,55.7])


for i,season in enumerate(("DJF", "MAM", "JJA", "SON")):
	ax[0,i].set_title(season)

cbar = fig.colorbar(cf1, ax=ax.ravel().tolist(),location="right", cmap = cmap, label = "%s [%s]" % (varname,temp_bal.units), shrink = 1)

#ax[0,0].set_ylabel("GLO (CMEMS) 2022-2023")
#ax[1,0].set_ylabel("BAL (CMEMS) 2022-2023")
#ax[0,0].set_ylabel("ciao",color="r")

# place a text box in upper left in axes coords (I cannot make set_ylabel to work....)
ax[0,0].text(-0.2, 0.5, "GLO (CMEMS) 2022-2023", transform=ax[0,0].transAxes, fontsize=10,
        verticalalignment='center', bbox=props, rotation=90)

ax[1,0].text(-0.2, 0.5, "BAL (CMEMS) 2022-2023", transform=ax[1,0].transAxes, fontsize=10,
        verticalalignment='center', bbox=props, rotation=90)

figname = "glo_bal_%s_level_%d_climatology" % (varname,depth_level)

plt.savefig(figname,bbox_inches="tight")

plt.show()
