import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import os

"""
	Note!!!
	changing the projection can be quite detrimental for this script!!!

	I suggest not to change the projection (leave PlateCarree). When I tried
	to set EuroPP I got a lot of errors from quiver function. Possibly there 
	are problem with vector transformation.

"""

outFolder = "maps"


props   = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
proj	= ccrs.EuroPP()

varname 	= "so"
cmap 		= "Spectral_r"
cmin,cmax 	= 0,60

ds_glo      = xr.open_dataset("../cmems_glob12/glo_anfc_subset.nc")
ds_bal      = xr.open_dataset("../cmems_baltic/bal_anfc_subset2.nc")

# load datasets and make season average
var_glo     = ds_glo[varname].groupby("time.season").mean(dim="time")
var_bal     = ds_bal[varname].groupby("time.season").mean(dim="time")

depth_glo   = ds_glo["depth"].data
depth_bal   = ds_bal["depth"].data

lons_glo	= ds_glo["longitude"]
lats_glo	= ds_glo["latitude"]
lons_bal	= ds_bal["lon"]
lats_bal	= ds_bal["lat"]

lons_glo,lats_glo	= np.meshgrid(lons_glo,lats_glo)
lons_bal,lats_bal	= np.meshgrid(lons_bal,lats_bal)


proj	= ccrs.PlateCarree()

nrows = 2
ncols = 4

fig,ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(17,8),subplot_kw={'projection': proj},layout="constrained")
fig.subplots_adjust(wspace=0.1, hspace=0.02)

clevels = [0,5,10,15,20,25,30,35,40,45,50,55,60]


for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):

	# here field glo is the differentiation along z of field
	field_glo = var_glo.sel(season=season).differentiate(coord="depth").data 
	# fill NaN with 0 (np.nanargmin does not handle arrays all made of Nans)	
	field_glo[np.isnan(field_glo)] = 0 
	# now field glo becomes the depth at where d/dz is max in absolute value
	field_glo = depth_glo[np.argmax(abs(field_glo),axis=0)]
	# mask shallowest values
	field_glo = np.ma.masked_where(field_glo <= depth_glo[0], field_glo )
    #cf1 = ax[0,i].pcolormesh(lons_glo,lats_glo,field_glo,vmin=cmin,vmax=cmax,cmap=cmap)
	cf1 = ax[0,i].contourf(lons_glo,lats_glo,field_glo,cmap=cmap,levels=clevels)

	# do the same for bal
	field_bal = var_bal.sel(season=season).differentiate(coord="depth").data
	field_bal[np.isnan(field_bal)] = 0 
	field_bal = depth_bal[np.argmax(abs(field_bal),axis=0)]
	field_bal = np.ma.masked_where(field_bal <= depth_bal[0],field_bal)

    #cf2 = ax[1,i].pcolormesh(lons_bal,lats_bal,field_bal,vmin=cmin,vmax=cmax,cmap=cmap)
	cf2 = ax[1,i].contourf(lons_bal,lats_bal,field_bal,cmap=cmap,levels=clevels)

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

clabel = r"Halocline Depth [m] $ d[max(\frac{\partial S}{\partial z})]$"

cbar = fig.colorbar(cf1, ax=ax.ravel().tolist(),location="right" , shrink = 0.8)
cbar.set_label(clabel, fontsize=15)

#ax[0,0].set_ylabel("GLO (CMEMS) 2022-2023")
#ax[1,0].set_ylabel("BAL (CMEMS) 2022-2023")
#ax[0,0].set_ylabel("ciao",color="r")

# place a text box in upper left in axes coords (I cannot make set_ylabel to work....)
ax[0,0].text(-0.2, 0.5, "GLO (CMEMS) 2022-2023", transform=ax[0,0].transAxes, fontsize=10,
        verticalalignment='center', bbox=props, rotation=90)

ax[1,0].text(-0.2, 0.5, "BAL (CMEMS) 2022-2023", transform=ax[1,0].transAxes, fontsize=10,
        verticalalignment='center', bbox=props, rotation=90)

figname = "cmems_glo_bal_halocline_climatology" 

plt.savefig(os.path.join(outFolder,figname),bbox_inches="tight")

#plt.imshow(np.argmax(np.abs(var_glo.differentiate(coord="depth").sel(season="JJA")).data,axis=0)[::-1,:])


