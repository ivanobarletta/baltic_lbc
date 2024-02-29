import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from datetime import datetime

"""
    Note:
    Purpose of the script: Create 2 Hovemoller from 2 regular datasets (time,depth,lat,lon)

    Set:
        1) varname = name of variable in dataset
        2) cmap = type of colormap
        3) cmin,cmax = limits of colorbar
        4) lon0,lat0 = coordinates of point where hovemoller is extracted
        
        5) set:
            5a) outFolder 
            5b) outFolderProfiles      

        6) profile_dates = list of datetimes where profile are extracted
            and saved to specific files.
            If you don't want profiles just set the list as empty ([])

        7) set method for the halocline or thermocline: 
            "max_grad", depth where the maximum(abs( \partial s / \partial z)) is located
            "threshold", depth where T[z]-T[z-dz] > 0.2
            
"""

def findFirstIndex(profile):
    """
    not really elegant function but I cannot have the first occurrence of T[k+1]-T[k]<-0.2 
    found using numpy. 
    """
    nz = profile.size
    for zlev in range(1,nz):
        if (profile.isel(depth=zlev).data - profile.isel(depth=zlev-1).data) < -0.2:
            break    
    return zlev



outFolder           = "hovemoller"
outFolderProfiles   = "profiles"

varname     = "so"
cmap        = "hot_r"    
cmin,cmax   = 5,20

varname     = "thetao"
cmap        = "Spectral_r"    
cmin,cmax   = 3,20

lon0,lat0   = 13.5,55

method      = "threshold"

max_depth   = 50

# dates where profiles are extracted and plotted in specific files
profile_dates = [datetime(2022,4,1),datetime(2022,7,29),datetime(2022,9,28)]
profile_dates = [datetime(2022,8,25),datetime(2023,8,7)]


ds_glo      = xr.open_dataset("../cmems_glob12/glo_anfc_subset.nc")
ds_bal      = xr.open_dataset("../cmems_baltic/bal_anfc_subset2.nc")

depth_glo   = ds_glo["depth"].data
depth_bal   = ds_bal["depth"].data

var_glo     = ds_glo[varname].interp(longitude=lon0,latitude=lat0)
var_bal     = ds_bal[varname].interp(lon=lon0,lat=lat0)

fig,ax      = plt.subplots(2,1,sharex=True,sharey=True,layout="constrained",figsize=(13,7))

# plot hovemollers
var_glo.plot(ax=ax[0],x="time",y="depth",cmap=cmap,vmin=cmin,vmax=cmax)
var_bal.plot(ax=ax[1],x="time",y="depth",cmap=cmap,vmin=cmin,vmax=cmax)

# calculate halocline with maximum of ds/dz 
if method == "max_grad":
    halocline_glo = depth_glo[np.nanargmax(abs(var_glo.differentiate(coord="depth")).data,axis=1)]
    halocline_bal = depth_bal[np.nanargmax(abs(var_bal.differentiate(coord="depth")).data,axis=1)]
if method == "threshold":
    print("Calculating Thermocline Depth ...")
    thermocline_glo = depth_glo[np.array([findFirstIndex(profile) for profile in var_glo])]
    thermocline_bal = depth_bal[np.array([findFirstIndex(profile) for profile in var_bal])]
    # mask too high values
    thermocline_glo = np.ma.masked_where(thermocline_glo > max_depth, thermocline_glo)
    thermocline_bal = np.ma.masked_where(thermocline_bal > max_depth, thermocline_bal)


if method == "max_grad":
    ax[0].plot(var_glo.time,halocline_glo,color="b")
    ax[1].plot(var_bal.time,halocline_bal,color="b")
if method == "threshold":    
    ax[0].plot(var_glo.time,thermocline_glo,color="0.3")
    ax[1].plot(var_bal.time,thermocline_bal,color="0.3")

ax[0].invert_yaxis()
ax[0].set_ylim(max_depth,0)


# projection of map
projj = ccrs.EuroPP()
ax_inset = fig.add_axes([0.03, 0.78, 0.1, 0.12], projection=projj)
ax_inset.set_extent([13,14,54.3,55.6])
ax_inset.coastlines(resolution='10m')
ax_inset.scatter(lon0,lat0,transform=ccrs.PlateCarree(),color="r",s=2)

figname = "%s_hovemoller_glo_bal_lon_%5.2f_lat_%5.2f.png" % (varname, lon0,lat0)

# plot vertical black lines corresponding to profiles extracted
for date in profile_dates:
    ax[0].vlines(x=date,ymin=max_depth,ymax=0,color="k")
    ax[1].vlines(x=date,ymin=max_depth,ymax=0,color="k")

plt.savefig(os.path.join(outFolder,figname),bbox_inches="tight",dpi=600)

label2      = r"$\frac{\partial (%s)}{\partial z}$" % varname

# create profile plots corresponding to selected dates in separated files
if method == "max_grad":
    for date in profile_dates:
        print(date)
        fig,ax = plt.subplots(figsize=(5,8))
        ax2 = ax.twiny()

        #plot global profile and ds/dz
        profile_glo = var_glo.interp(time=date,method="nearest")
        dpdz_glo    = -profile_glo.differentiate(coord="depth")

        ax.plot(profile_glo,profile_glo.depth,color="k",label="GLO");ax.scatter(profile_glo,profile_glo.depth,color="k")
        ax2.plot(dpdz_glo,dpdz_glo.depth,color="b");ax2.scatter(dpdz_glo,dpdz_glo.depth,color="b")

        #plot baltic profile and ds/dz
        profile_bal = var_bal.interp(time=date,method="nearest")
        dpdz_bal    = -profile_bal.differentiate(coord="depth")

        ax.plot(profile_bal,profile_bal.depth,color="k",label="BAL",linestyle="-.");ax.scatter(profile_bal,profile_bal.depth,color="k")
        ax2.plot(dpdz_bal,dpdz_bal.depth,color="b",linestyle="-.");ax2.scatter(dpdz_bal,dpdz_bal.depth,color="b")

        ax.invert_yaxis()

        ax2.set_xlabel(label2,color="b",fontsize=15)
        ax2.tick_params(axis="x",labelcolor="b")
        ax.set_xlabel("%s [%s]" % (varname,var_bal.units),color="k",fontsize=15)
        ax.tick_params(axis="x",labelcolor="k")
        ax.set_ylabel("Depth [m]", fontsize = 15)
        ax.grid()

        halocline_depth_glo = dpdz_glo.depth[np.nanargmax(abs(dpdz_glo))].data
        halocline_depth_bal = dpdz_bal.depth[np.nanargmax(abs(dpdz_bal))].data

        xminn, xmaxx = min(np.nanmin(profile_glo),np.nanmin(profile_bal)), max(np.nanmax(profile_glo),np.nanmax(profile_bal))
        ax.hlines(y=halocline_depth_glo,color="0.3",xmin=xminn,xmax=xmaxx)
        ax.hlines(y=halocline_depth_bal,color="0.3",xmin=xminn,xmax=xmaxx,linestyle="-.")

        ax.set_ylim(max_depth,0)    

        ax.legend(loc="lower left")

        fig.suptitle("Profiles - %s - (%5.2fE,%5.2fN)" % (date.strftime("%Y%m%d"),lon0,lat0))

        figname = "profile_%s_%s_lon_%5.2f_lat_%5.2f.png" % (varname,date.strftime("%Y%m%d"),lon0,lat0)
        plt.savefig(os.path.join(outFolderProfiles,figname))

if method == "threshold":
    for date in profile_dates:
        print(date)
        fig,ax = plt.subplots(figsize=(5,8))

        #plot global profile and ds/dz
        profile_glo = var_glo.interp(time=date,method="nearest")

        ax.plot(profile_glo,profile_glo.depth,color="k",label="GLO");ax.scatter(profile_glo,profile_glo.depth,color="k")

        #plot baltic profile and ds/dz
        profile_bal = var_bal.interp(time=date,method="nearest")

        ax.plot(profile_bal,profile_bal.depth,color="k",label="BAL",linestyle="-.");ax.scatter(profile_bal,profile_bal.depth,color="k")

        ax.invert_yaxis()

        ax.set_xlabel("%s [%s]" % (varname,var_bal.units),color="k",fontsize=15)
        ax.tick_params(axis="x",labelcolor="k")
        ax.set_ylabel("Depth [m]", fontsize = 15)
        ax.grid()

        thermocline_depth_glo = depth_glo[findFirstIndex(profile_glo)]
        thermocline_depth_bal = depth_bal[findFirstIndex(profile_bal)]

        xminn, xmaxx = min(np.nanmin(profile_glo),np.nanmin(profile_bal)), max(np.nanmax(profile_glo),np.nanmax(profile_bal))
        ax.hlines(y=thermocline_depth_glo,color="0.3",xmin=xminn,xmax=xmaxx)
        ax.hlines(y=thermocline_depth_bal,color="0.3",xmin=xminn,xmax=xmaxx,linestyle="-.")

        ax.set_ylim(max_depth,0)    

        ax.legend(loc="lower left")

        fig.suptitle("Profiles - %s - (%5.2fE,%5.2fN)" % (date.strftime("%Y%m%d"),lon0,lat0))

        figname = "profile_%s_%s_lon_%5.2f_lat_%5.2f.png" % (varname,date.strftime("%Y%m%d"),lon0,lat0)
        plt.savefig(os.path.join(outFolderProfiles,figname))


#plt.show()