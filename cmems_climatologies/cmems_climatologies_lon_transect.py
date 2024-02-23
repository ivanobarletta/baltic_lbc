import xarray as xr
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
import numpy as np
import warnings

"""

	meaning of the box list

	-----------------------------	- box[3]			^
	|							|						|
	|							|						|
	|							|						|
	|							|	lat				lat average
	|							|						|
	|							|						|
	-----------------------------	- box[2]			v
	|							|

	box[0]		lon				box[1]					

	the script computes seasonal climatologies of
	2 cmems datatasets and a spatial average along latitude
	
	the resulting datasets has dimensions

	(season,depth,longitude)

	and transect plots are displayed, showing the input variable
	and the (u,w) velocity field overlapped as quiver

			ds1			ds2
		---------	---------
	DJF	|		|	|		|	depth
		|		|	|		|
		---------	---------
		---------	---------
	MAM	|		|	|		|	depth
		|		|	|		|
		---------	---------
		---------	---------
	JJA	|		|	|		|	depth
		|		|	|		|
		---------	---------
		---------	---------
	SON	|		|	|		|	depth
		|		|	|		|
		---------	---------
			lon			lon

"""

outFolder 		= "transects"

wAmpliFactor	= 1000.0	# amplification factor for w component of velocity

box				= [13, 14, 54.6, 55]

ds_glo = xr.open_dataset("../cmems_glob12/glo_anfc_subset.nc")
ds_bal = xr.open_dataset("../cmems_baltic/bal_anfc_subset2.nc")

# I do the climatology here 
print("Computing Climatologies and average along lat")
#ds_glo = ds_glo.groupby("time.season").mean(dim="time")
#ds_bal = ds_bal.groupby("time.season").mean(dim="time")

ds_glo = ds_glo.groupby("time.season").mean(dim="time").sel(latitude=slice(box[2],box[3])).mean(dim="latitude",keep_attrs=True)
ds_bal = ds_bal.groupby("time.season").mean(dim="time").sel(lat=slice(box[2],box[3])).mean(dim="lat",keep_attrs=True)


def makeLonTransectPlot(varname="varname",cmap="jet",lat_interp=0,cmin=0,cmax=30):

	var_glo = ds_glo[varname]
	var_bal = ds_bal[varname]

	var_glo_season = var_glo.groupby("time.season").mean(dim="time")
	var_bal_season = var_bal.groupby("time.season").mean(dim="time")

	nrows = 4
	ncols = 2

	fig,ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,figsize=(10,13))
	#fig,ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(7,11))
	fig.subplots_adjust(wspace=0.1, hspace=0.05)

	for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
		cf1 	= var_glo_season.interp(latitude=lat_interp).sel(season=season).plot(ax=ax[i,0],vmin=cmin,vmax=cmax,extend="both",cmap=cmap,add_colorbar=False)
		cf2 	= var_bal_season.interp(lat=lat_interp).sel(season=season).plot(ax=ax[i,1],vmin=cmin,vmax=cmax,extend="both",cmap=cmap,add_colorbar=False)

	for j in range(nrows):
		for i in range(ncols):
			ax[j,i].set_facecolor('grey')

			ax[j,i].set_ylabel("")
			ax[j,i].set_xlabel("")

			ax[j,i].set_title("")

	ax[0,0].set_title("GLO (CMEMS) 2022-2023")
	ax[0,1].set_title("BAL (CMEMS) 2022-2023")

	ax[0,0].set_xlim(13,14)
	ax[0,0].invert_yaxis()
	ax[0,0].set_ylim(50,0)

	ax[-1,0].tick_params(labelrotation=45)
	ax[-1,1].tick_params(labelrotation=45)

	ax[-1,0].set_xlabel("Lon")
	ax[-1,1].set_xlabel("Lon")

	for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
		ax[i,0].set_ylabel(season)


	cbar = fig.colorbar(cf1, ax=ax.ravel().tolist(),location="right", cmap = cmap, shrink = 1)
	clabel = "%s [lat = %5.2f][%s]" % (varname,lat_interp,var_bal.units)
	cbar.set_label(clabel, fontsize=15)
	#cbar.ax.tick_params(labelsize=25)

	#ax[0,0].set_ylabel("GLO (CMEMS) 2022-2023")
	#ax[1,0].set_ylabel("BAL (CMEMS) 2022-2023")
	#ax[0,0].set_ylabel("ciao",color="r")

	# place a text box in upper left in axes coords (I cannot make set_ylabel to work....)
	"""
	ax[0,0].text(-0.2, 0.5, "GLO (CMEMS) 2022-2023", transform=ax[0,0].transAxes, fontsize=10,
			verticalalignment='center', bbox=props, rotation=90)

	ax[1,0].text(-0.2, 0.5, "BAL (CMEMS) 2022-2023", transform=ax[1,0].transAxes, fontsize=10,
			verticalalignment='center', bbox=props, rotation=90)
	"""

	figname = "glo_bal_%s_transect_lat_%5.2f_climatology.png" % (varname,lat_interp)

	plt.savefig(os.path.join("transects",figname),bbox_inches="tight")

def makeLonTransectPlot2(varname="varname",cmap="jet",lat_interp=0,cmin=0,cmax=30):
	"""
	the same function but with a small map displaying
	the transect line
	"""

	var_glo = ds_glo[varname]
	var_bal = ds_bal[varname]

	var_glo_season = var_glo.groupby("time.season").mean(dim="time")
	var_bal_season = var_bal.groupby("time.season").mean(dim="time")

	nrows = 4
	ncols = 2

	fig,ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,figsize=(11,13),layout="constrained")
	#fig,ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(7,11))
	fig.subplots_adjust(wspace=0.1, hspace=0.05)

	for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
		cf1 	= var_glo_season.interp(latitude=lat_interp).sel(season=season).plot(ax=ax[i,0],vmin=cmin,vmax=cmax,extend="both",cmap=cmap,add_colorbar=False)
		cf2 	= var_bal_season.interp(lat=lat_interp).sel(season=season).plot(ax=ax[i,1],vmin=cmin,vmax=cmax,extend="both",cmap=cmap,add_colorbar=False)

	for j in range(nrows):
		for i in range(ncols):
			ax[j,i].set_facecolor('grey')

			ax[j,i].set_ylabel("")
			ax[j,i].set_xlabel("")

			ax[j,i].set_title("")

	ax[0,0].set_title("GLO (CMEMS) 2022-2023")
	ax[0,1].set_title("BAL (CMEMS) 2022-2023")

	ax[0,0].set_xlim(13,14)
	ax[0,0].invert_yaxis()
	ax[0,0].set_ylim(50,0)

	ax[-1,0].tick_params(labelrotation=45)
	ax[-1,1].tick_params(labelrotation=45)

	ax[-1,0].set_xlabel("Lon")
	ax[-1,1].set_xlabel("Lon")

	for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
		ax[i,0].set_ylabel(season)



	cbar = fig.colorbar(cf1, ax=ax[0:2,1],location="right", cmap = cmap, shrink = 1)
	clabel = "%s [lat = %5.2f][%s]" % (varname,lat_interp,var_bal.units)
	cbar.set_label(clabel, fontsize=15)
	#cbar.ax.tick_params(labelsize=25)

	#ax[0,0].set_ylabel("GLO (CMEMS) 2022-2023")
	#ax[1,0].set_ylabel("BAL (CMEMS) 2022-2023")
	#ax[0,0].set_ylabel("ciao",color="r")

	# place a text box in upper left in axes coords (I cannot make set_ylabel to work....)
	"""
	ax[0,0].text(-0.2, 0.5, "GLO (CMEMS) 2022-2023", transform=ax[0,0].transAxes, fontsize=10,
			verticalalignment='center', bbox=props, rotation=90)

	ax[1,0].text(-0.2, 0.5, "BAL (CMEMS) 2022-2023", transform=ax[1,0].transAxes, fontsize=10,
			verticalalignment='center', bbox=props, rotation=90)
	"""

	ax_box = fig.add_axes([0.8,0.05,0.2,0.2],projection=ccrs.PlateCarree())
	ax_box.hlines(y=lat_interp,xmin=13,xmax=14,transform=ccrs.PlateCarree(),color="r")
	ax_box.coastlines(resolution="10m")
	ax_box.set_extent([12.5,14.5,54,56])

	figname = "2glo_bal_%s_transect_lat_%5.2f_climatology.png" % (varname,lat_interp)

	plt.savefig(os.path.join("transects",figname),bbox_inches="tight")

def makeLonTransectPlot3(varname="varname",cmap="jet",lat_interp=0,cmin=0,cmax=30,box=[0,1,0,1]):
	"""
	the same as makeLonTransectPlot but with a small map displaying
	the transect line and velocity field (u,w)
	"""

	if box == [0,1,0,1]:
		warnings.warn("Please Set box [lon1,lon2,lat1,lat2] \n You will have unespected results!")		

	var_glo = ds_glo[varname]
	var_bal = ds_bal[varname]

	u_glo	= ds_glo["uo"]
	w_glo	= ds_glo["wo"]
	u_bal	= ds_bal["uo"]
	w_bal	= ds_bal["wo"]

	depth_glo	= ds_glo["depth"]
	depth_bal	= ds_bal["depth"]

	lon_glo		= ds_glo["longitude"]
	lon_bal		= ds_bal["lon"]

	x_glo,y_glo	= np.meshgrid(lon_glo,depth_glo)
	x_bal,y_bal	= np.meshgrid(lon_bal,depth_bal)

	"""
	var_glo_season 	= var_glo.groupby("time.season").mean(dim="time")
	var_bal_season 	= var_bal.groupby("time.season").mean(dim="time")

	u_glo_season	= u_glo.groupby("time.season").mean(dim="time")	
	w_glo_season	= w_glo.groupby("time.season").mean(dim="time")
	u_bal_season	= u_bal.groupby("time.season").mean(dim="time")	
	w_bal_season	= w_bal.groupby("time.season").mean(dim="time")
	"""	
	
	nrows = 4
	ncols = 2

	"""
	print(x_glo.shape)
	print(y_glo.shape)
	print(u_glo_season.sel(season="DJF").shape)
	print(w_glo_season.sel(season="DJF").shape)
	"""
	
	fig,ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,figsize=(11,13),layout="constrained")
	#fig,ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(7,11))
	fig.subplots_adjust(wspace=0.1, hspace=0.05)

	for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
		cf1 	= var_glo.sel(season=season).plot(ax=ax[i,0],vmin=cmin,vmax=cmax,extend="both",cmap=cmap,add_colorbar=False)
		# add currents pattern
		uu = u_glo.sel(season=season)
		ww = w_glo.sel(season=season)
		ax[i,0].quiver(x_glo,y_glo,uu,wAmpliFactor*ww,scale=0.5)
		cf2 	= var_bal.sel(season=season).plot(ax=ax[i,1],vmin=cmin,vmax=cmax,extend="both",cmap=cmap,add_colorbar=False)
		uu = u_bal.sel(season=season)
		ww = w_bal.sel(season=season)
		ax[i,1].quiver(x_bal,y_bal,uu,wAmpliFactor*ww,scale=0.5)

	for j in range(nrows):
		for i in range(ncols):
			ax[j,i].set_facecolor('grey')

			ax[j,i].set_ylabel("")
			ax[j,i].set_xlabel("")

			ax[j,i].set_title("")

	ax[0,0].set_title("GLO (CMEMS) 2022-2023")
	ax[0,1].set_title("BAL (CMEMS) 2022-2023")

	ax[0,0].set_xlim(box[0],box[1])
	ax[0,0].invert_yaxis()
	ax[0,0].set_ylim(50,0)

	ax[-1,0].tick_params(labelrotation=45)
	ax[-1,1].tick_params(labelrotation=45)

	ax[-1,0].set_xlabel("Lon")
	ax[-1,1].set_xlabel("Lon")

	for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
		ax[i,0].set_ylabel(season)


	cbar = fig.colorbar(cf1, ax=ax[0:2,1],location="right", cmap = cmap, shrink = 1)
	clabel = "%s [lat = <%5.2f-%5.2f>][%s]" % (varname,box[2],box[3],var_bal.units)
	cbar.set_label(clabel, fontsize=15)
	#cbar.ax.tick_params(labelsize=25)

	#ax[0,0].set_ylabel("GLO (CMEMS) 2022-2023")
	#ax[1,0].set_ylabel("BAL (CMEMS) 2022-2023")
	#ax[0,0].set_ylabel("ciao",color="r")

	# place a text box in upper left in axes coords (I cannot make set_ylabel to work....)
	"""
	ax[0,0].text(-0.2, 0.5, "GLO (CMEMS) 2022-2023", transform=ax[0,0].transAxes, fontsize=10,
			verticalalignment='center', bbox=props, rotation=90)

	ax[1,0].text(-0.2, 0.5, "BAL (CMEMS) 2022-2023", transform=ax[1,0].transAxes, fontsize=10,
			verticalalignment='center', bbox=props, rotation=90)
	"""
	
	# add box showing the transect line
	ax_box = fig.add_axes([0.8,0.05,0.2,0.2],projection=ccrs.PlateCarree())
	#ax_box.hlines(y=lat_interp,xmin=13,xmax=14,transform=ccrs.PlateCarree(),color="r")
	ax_box.plot(np.asarray(box)[[0,1,1,0,0]],np.asarray(box)[[2,2,3,3,2]],color="r",transform=ccrs.PlateCarree())
	ax_box.coastlines(resolution="10m")
	ax_box.set_extent([12.5,14.5,54,56])

	figname = "3glo_bal_%s_transect_lat1-lat2_%5.2f_%5.2f_climatology.png" % (varname,box[2],box[3])

	plt.savefig(os.path.join("transects",figname),bbox_inches="tight")

#makeLonTransectPlot(varname="thetao",cmap="Spectral_r"	,lat_interp=55,cmin=0,cmax=20)
#makeLonTransectPlot(varname="so"	,cmap="hot_r"		,lat_interp=55,cmin=5,cmax=15)
#makeLonTransectPlot2(varname="thetao"	,cmap="Spectral_r"	,lat_interp=55,cmin=0,cmax=20)
#makeLonTransectPlot2(varname="so"		,cmap="hot_r"		,lat_interp=55,cmin=5,cmax=15)
makeLonTransectPlot3(varname="thetao"	,cmap="Spectral_r"	,lat_interp=55,cmin=0,cmax=20,box=box)
makeLonTransectPlot3(varname="so"		,cmap="hot_r"		,lat_interp=55,cmin=5,cmax=15,box=box)


#plt.show()