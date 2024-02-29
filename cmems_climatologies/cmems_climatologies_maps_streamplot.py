import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import cmocean
import os
from datetime import datetime

"""
    Purpose: 
        make streamplot for circulation ad selected
        time / level

        The dataset are interpolated onto evenly-spaced grid
        to avoid related streamplot error 
"""
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

kw={"fill_value": "extrapolate"}

dates   = [datetime(2022,9,29)]
depths  = [0,40]

proj_plot   = ccrs.PlateCarree()
projj   = ccrs.PlateCarree()

outFolder = "maps"

ds_glo  = xr.open_dataset("../cmems_glob12/glo_anfc_subset.nc")
ds_bal  = xr.open_dataset("../cmems_baltic/bal_anfc_subset2.nc")

# create evenly spaced coords from glo coordinates
xi,yi   = np.linspace(ds_glo.longitude.min(),ds_glo.longitude.max(),ds_glo.longitude.size),np.linspace(ds_glo.latitude.min(),ds_glo.latitude.max(),ds_glo.latitude.size) 
ds_glo = ds_glo.interp(time=dates,method="nearest").interp(depth=depths,kwargs=kw).interp(longitude=xi,latitude=yi)

# create evenly spaced coords from bal coordinates
xi,yi   = np.linspace(ds_bal.lon.min(),ds_bal.lon.max(),ds_bal.lon.size),np.linspace(ds_bal.lat.min(),ds_bal.lat.max(),ds_bal.lat.size) 
ds_bal = ds_bal.interp(time=dates,method="nearest").interp(depth=depths,kwargs=kw).interp(lon=xi,lat=yi)

for t in range(len(dates)):
    #fig,ax = plt.subplots(nrows=len(depths),ncols=2,sharex=True,sharey=True,layout="constrained"
    #                      ,subplot_kw={"projection":proj_plot},figsize=(10,10))

    fig,ax = plt.subplots(nrows=len(depths),ncols=2,sharex=True,sharey=True,subplot_kw={"projection":proj_plot},figsize=(10,10))
    
    fig.subplots_adjust(wspace=0.02, hspace=0.02)

    ax[0,0].set_extent([13,14,54.5,55.5])

    for d in range(len(depths)):

        ax[d,0].coastlines(resolution='10m')
        ax[d,1].coastlines(resolution='10m')

        ds_glo.isel(time=t,depth=d).plot.streamplot(ax=ax[d,0],transform=projj,x="longitude",y="latitude",u="uo",v="vo",density=1,linewidth=0.2)
        ds_bal.isel(time=t,depth=d).plot.streamplot(ax=ax[d,1],transform=projj,x="lon",y="lat",u="uo",v="vo",density=1,linewidth=0.2)

        #ds_glo.isel(time=t,depth=d).plot.streamplot(ax=ax[d,0],x="longitude",y="latitude",u="uo",v="vo"
        #                                            ,density=3,linewidth=0.2,color="k",transform=ccrs.PlateCarree())
        #ds_bal.isel(time=t,depth=d).plot.streamplot(ax=ax[d,1],x="lon"      ,y="lat"     ,u="uo",v="vo"
        #                                            ,density=2,linewidth=0.2,color="k")

        ax[d,0].set_title("")
        ax[d,1].set_title("")


        gl = ax[d,0].gridlines(crs=projj, draw_labels=False,
    		linewidth=1, color='gray', alpha=0.5, linestyle='--')
		
        gl.xlabels_top = False
        gl.ylabels_left = False
        gl.xlabels_bottom = False
        gl.ylabels_right = False

        gl = ax[d,1].gridlines(crs=projj, draw_labels=False,
    		linewidth=1, color='gray', alpha=0.5, linestyle='--')
		
        gl.xlabels_top = False
        gl.ylabels_left = False
        gl.xlabels_bottom = False
        gl.ylabels_right = False    

    for d in range(len(depths)):
        # place a text box in upper left in axes coords (I cannot make set_ylabel to work....)
        ax[d,0].text(-0.2, 0.5, "Depth = %d [m]" % depths[d], transform=ax[d,0].transAxes, fontsize=15,
            verticalalignment='center', bbox=props, rotation=90)

    ax[0,0].set_title("GLO (Cmems)")
    ax[0,1].set_title("BAL (Cmems)")

    fig.suptitle("Currents - %s" % dates[t].strftime("%Y%m%d"))


    figname = "glo_bal_streamplot_%s" % (dates[t].strftime("%Y%m%d"))

    plt.savefig(os.path.join(outFolder,figname))

plt.show()
