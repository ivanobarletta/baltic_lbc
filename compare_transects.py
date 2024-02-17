import xarray as xr  
import matplotlib.pyplot as plt  
import numpy as np  
from lbc_utils import get_nemo_transect_coords

xidx = 7

varName = "vozocrtx"
varName = "votemper"

depth_list = [0.4940254, 1.541375, 2.645669, 3.819495, 5.078224, 6.440614, 
    7.92956, 9.572997, 11.405, 13.46714, 15.81007, 18.49556, 21.59882, 
    25.21141, 29.44473, 34.43415, 40.34405, 47.37369, 55.76429, 65.80727, 
    77.85385, 92.32607, 109.7293, 130.666, 155.8507, 186.1256, 222.4752, 
    266.0403, 318.1274, 380.213, 453.9377, 541.0889, 643.5668, 763.3331, 
    902.3393, 1062.44, 1245.291, 1452.251, 1684.284, 1941.893, 2225.078, 
    2533.336, 2865.703, 3220.82, 3597.032, 3992.484, 4405.224, 4833.291, 
    5274.784, 5727.917 ]   

sirenLBCPath = "neatl_lbc_files/NEATL36_obcdta_east_2_20220115P01_R20220123.nc"
cdoLBCPath = "cdo_create_LBC/test_out_cons2.nc"
xesmfLBCPath = "xesmf_create_LBC/test_out_xesmf.nc"
tmaskPath = "static_files/mask_gridT_east2.nc"
cmemsPath = "cmems_baltic/BAL-NEMO_PHY-DailyMeans-20220115.nc"

# open datasets
dsSIREN = xr.open_dataset(sirenLBCPath)
dsCDO   = xr.open_dataset(cdoLBCPath)
dsXESMF = xr.open_dataset(xesmfLBCPath)
dstmask = xr.open_dataset(tmaskPath)
dsCMEMS = xr.open_dataset(cmemsPath)

lats = dsSIREN.nav_lat.isel(X=xidx)

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

x_int, y_int = get_nemo_transect_coords(sirenLBCPath,xidx=xidx,transect_name="lat_along_transect")

varSIREN    = dsSIREN[varName].isel(T=0,X=xidx)
varCDO      = dsCDO[varName].isel(T=0,X=xidx)
varXESMF    = dsXESMF[varName].isel(T=0,X=xidx)
tmask       = dstmask["mask"].isel(time_counter=0,x=xidx)
varCMEMS    = dsCMEMS[guess_cmems_varname(varName=varName)].isel(time=0).interp(lon=x_int,lat=y_int,method="nearest")
depthCMEMS  = dsCMEMS["depth"]

# some stats
print("%s,min,max " % varName)
print("siren    ",np.nanmin(varSIREN),np.nanmax(varSIREN))
print("CDO      ",np.nanmin(varCDO)  ,np.nanmax(varCDO))
print("xesmf    ",np.nanmin(varXESMF),np.nanmax(varXESMF))
print("bal-cmems",np.nanmin(varCMEMS),np.nanmax(varCMEMS))

deptht = dsSIREN.deptht
date_str = np.datetime_as_string(dsCDO.time_counter.isel(T=0),unit="h")

def guess_cmap(varName):
    cmap = "jet"
    if varName in ["vozocrtx","vomecrty"]:
        cmap = "RdBu_r"
    return cmap    

def calc_min_max(v1,v2):
    min1 = np.nanmin(v1)
    max1 = np.nanmax(v1)
    min2 = np.nanmin(v2)
    max2 = np.nanmax(v2)
    print(min1,max1)
    print(min2,max2)

    # make shift in case of negative values
    if np.any(np.array([min1,min2])):
        print("ciao")
        shift = np.abs(min(min1,min2))
    else:
        shift = 0    
    print("shift:", shift)

    min1 += shift
    max1 += shift
    min2 += shift
    max2 += shift

    # round to nearest integer
    min12 = np.round(min(min1,min2))
    print()
    print(min12)
    if abs(min12) < 1:
        min12 = np.around(min(min1,min2))
    max12 = np.round(max(max1,max2))
    if abs(max12) < 1:
        max12 = np.around(max(max1,max2))

    print(max12)

    min12 -= shift
    max12 -= shift

    return min12,max12

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

lats2,depth2 = np.meshgrid(lats,np.asarray(depth_list))

lats3,depth3 = np.meshgrid(y_int,depthCMEMS)

print(varSIREN.shape)
vminn,vmaxx = calc_min_max(varSIREN,varCDO)



print("Vmin,vmax:",vminn,vmaxx)

cmapp = guess_cmap(varName=varName)
clabel = guess_clabel(varName=varName)

vminn,vmaxx = calc_min_max(varSIREN,varCDO)

if varName in ["vozocrtx","vomecrty"]:
    vminn = -1
    vmaxx = 1

fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(10,6))

#cf = varSIREN.plot(ax=ax[0],y=depth2,x=lats2,cmap=cmapp,vmin=vminn,vmax=vmaxx)
cf0 = ax[0].pcolormesh(lats2,depth2,varSIREN.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)
cf1 = ax[1].pcolormesh(lats2,depth2,varCDO.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)
cf2 = ax[2].pcolormesh(lats2,depth2,varXESMF.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)
cf3 = ax[3].pcolormesh(lats3,depth3,varCMEMS.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)

cbar = fig.colorbar(cf0, ax=ax[:4], shrink=0.3, location='bottom',label = clabel)

figtitle = "%s - time: %s" % (varName,date_str)
plt.suptitle(figtitle)

ax[0].set_title("Siren")
ax[1].set_title("CDO")
ax[2].set_title("XESMF")
ax[3].set_title("BAL-CMEMS\n (nearest)")

ax[0].set_ylabel("Depth [m]")
ax[1].set_xlabel("Lat along Transect")

plt.gca().invert_yaxis()
plt.ylim((50,0))

filename = "plots/transect_plots/transect_%s_%s" % (varName,date_str)
plt.savefig(filename,dpi=600,bbox_inches="tight")

#######################################
### same plots with Tmask

fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(10,6))

#cf = varSIREN.plot(ax=ax[0],y=depth2,x=lats2,cmap=cmapp,vmin=vminn,vmax=vmaxx)
cf0 = ax[0].pcolormesh(lats2,depth2,varSIREN.data*tmask.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)
cf1 = ax[1].pcolormesh(lats2,depth2,varCDO.data*tmask.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)
cf2 = ax[2].pcolormesh(lats2,depth2,varXESMF.data*tmask.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)
cf3 = ax[3].pcolormesh(lats3,depth3,varCMEMS.data,cmap=cmapp,vmin=vminn,vmax=vmaxx)

cbar = fig.colorbar(cf0, ax=ax[:4], shrink=0.3, location='bottom',label = clabel)

figtitle = "%s - time: %s" % (varName,date_str)
plt.suptitle(figtitle)

ax[0].set_title("Siren")
ax[1].set_title("CDO")
ax[2].set_title("XESMF")
ax[3].set_title("BAL-CMEMS \n (nearest)")


ax[0].set_ylabel("Depth [m]")
ax[1].set_xlabel("Lat along Transect")

plt.gca().invert_yaxis()
plt.ylim((50,0))

filename = "plots/transect_plots/masked_transect_%s_%s" % (varName,date_str)
plt.savefig(filename,dpi=600,bbox_inches="tight")

### same plots with Tmask
