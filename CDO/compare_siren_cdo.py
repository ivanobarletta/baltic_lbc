import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy as np

varName = "votemper"
#varName = "vozocrtx"
iLevel = 10
method = "conservative"

pathSiren = "../neatl_lbc_files/NEATL36_obcdta_east_2_20220115P01_R20220123.nc"

pathCDO = "test_out4.nc"
pathCDO = "test_out_cons.nc"
pathCDO = "test_out_cons2.nc"

pathCMEMS = "../cmems_baltic/BAL-NEMO_PHY-DailyMeans-20220115.nc"

dsSiren = xr.open_dataset(pathSiren)
dsCDO   = xr.open_dataset(pathCDO)
dsCMEMS = xr.open_dataset(pathCMEMS)




def guess_cmems_varname(varName):
    varName2 = ""
    if varName == "votemper":
        varName2 = "thetao"
    if varName == "vosaline":
        varName2 = "so"
    if varName == "vozocrtx":
        varName2 = "uo"
    if varName == "vomecrty":
        varName2 = "vo"
    
    return varName2

lons = dsSiren.nav_lon
lats = dsSiren.nav_lat
deptht = dsSiren.deptht


varSiren = dsSiren[varName].isel(T=0,Z=iLevel)
varCDO  = dsCDO[varName].isel(T=0,Z=iLevel)

varName2 = guess_cmems_varname(varName=varName)
kw = {"fill_value":"extrapolate"}
varCMEMS = dsCMEMS[varName2].isel(time=0).interp(depth=deptht.isel(Z=iLevel),kwargs=kw)
lons_cmems = varCMEMS.lon
lats_cmems = varCMEMS.lat

def guess_clabel(varName):
    clabel = " "
    if varName == "votemper": 
        clabel = "C"
    if varName == "vosaline":
        clabel = "PSU"
    if varName in ["vozocrtx","vomecrty"]:
        clabel = "m/s"
    if varName == "sossheig":
        clabel = "m"    
    clabel = "%s [%s]" % (varName,clabel)                
    return clabel

def guess_cmap(varName):
    cmap = "jet"
    if varName in ["vozocrtx","vomecrty"]:
        cmap = "RdBu"
    return cmap    

def calc_min_max(v1,v2):
    min1 = np.nanmin(v1)
    max1 = np.nanmax(v1)
    min2 = np.nanmin(v2)
    max2 = np.nanmax(v2)
    # round to nearest integer
    min12 = np.round(min(min1,min2))
    if abs(min12) < 1:
        min12 = np.around(min(min1,min2))
    max12 = np.round(max(max1,max2))
    if abs(max12) < 1:
        max12 = np.around(max(max1,max2))
    return min12,max12

vminn,vmaxx = calc_min_max(varSiren,varCDO)
cmapp = guess_cmap(varName=varName)
clabel = guess_clabel(varName=varName)
# polygon plot
def get_poly_corners(lons,lats):
    lons = lons.values
    lats = lats.values
    idx = ((0,0),(0,-1),(-1,-1),(-1,0),(0,0))
    x = np.asarray([lons[i,j] for i,j in idx])
    y = np.asarray([lats[i,j] for i,j in idx])

    return x,y

xx,yy = get_poly_corners(lons,lats)

print(xx)
print (yy)

fig,ax = plt.subplots(1,3,sharex=True,sharey=True,figsize=(14,5))
varSiren.plot(ax=ax[0],cmap=cmapp,x="nav_lon",y="nav_lat",vmin=vminn,vmax=vmaxx)
varCDO.plot(ax=ax[1],cmap=cmapp,x="nav_lon",y="nav_lat",vmin=vminn,vmax=vmaxx)
varCMEMS.plot(ax=ax[2],cmap=cmapp,vmin=vminn,vmax=vmaxx)
#cf0 = ax[0].pcolormesh(lons,lats,varSiren.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)
#cf1 = ax[1].pcolormesh(lons,lats,varCDO.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)
#cf2 = ax[2].pcolormesh(lons_cmems,lats_cmems,varCMEMS.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)
ax[2].plot(xx,yy,color="k",zorder=1)
ax[0].set_xlim(13,14)
ax[0].set_ylim(54.4,55.5)
date_str = np.datetime_as_string(dsCDO.time_counter.isel(T=0),unit="h")

figtitle = "%s - depth %s - time: %s - method: %s" % (varName,deptht.isel(Z=iLevel).values,date_str,method)
plt.suptitle(figtitle)
ax[0].set_title("Siren")
ax[1].set_title("CDO")
ax[2].set_title("CMEMS")

#cb = fig.colorbar(cf2, ax=ax.ravel().tolist())
#cb.set_label(clabel)

filename = "%s_level_%d_%s" % (varName,iLevel,method)

plt.savefig(filename,dpi=600,bbox_inches="tight")



