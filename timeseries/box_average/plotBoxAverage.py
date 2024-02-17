import matplotlib.pyplot as plt 
import xarray as xr 
import numba as np 
import os

varname         = "votemper"
varname         = "vosaline"
indexBox        = [3,13,31,41]
zindex          = 5 

rootFolder      = "../lbc_files_concatenated"

pathSiren       = "NEATL36_obcdta_east2_2022-2023.nc"
pathCDO_GLO     = "NEATL36_east2_CDO_GLO_2022-2023.nc"
pathCDO_BAL     = "NEATL36_east2_CDO_BAL_2022-2023.nc"
pathXESMF_BAL   = "NEATL36_east2_XESMF_BAL_2022-2023.nc"

dsSiren         = xr.open_dataset(os.path.join(rootFolder,pathSiren))
dsCDO_BAL       = xr.open_dataset(os.path.join(rootFolder,pathCDO_BAL))
dsCDO_GLO       = xr.open_dataset(os.path.join(rootFolder,pathCDO_GLO))
dsXESMF_BAL     = xr.open_dataset(os.path.join(rootFolder,pathXESMF_BAL))

deptht          = dsSiren["deptht"]

varSiren        = dsSiren[varname].isel(Z=zindex)
varCDO_BAL      = dsCDO_BAL[varname].isel(Z=zindex)
varCDO_GLO      = dsCDO_GLO[varname].isel(Z=zindex)
varXESMF_BAL    = dsXESMF_BAL[varname].isel(Z=zindex)

timeSiren           = dsSiren["time_counter"]
timeSeriesSiren     = varSiren.isel(X=slice(indexBox[0],indexBox[1]),Y=slice(indexBox[2],indexBox[3])).mean(dim=["X","Y"])

timeCDO_BAL         = dsCDO_BAL["time_counter"]
timeSeriesCDO_BAL   = varCDO_BAL.isel(X=slice(indexBox[0],indexBox[1]),Y=slice(indexBox[2],indexBox[3])).mean(dim=["X","Y"])

timeCDO_GLO         = dsCDO_GLO["time_counter"]
timeSeriesCDO_GLO   = varCDO_GLO.isel(X=slice(indexBox[0],indexBox[1]),Y=slice(indexBox[2],indexBox[3])).mean(dim=["X","Y"])

timeXESMF_BAL         = dsXESMF_BAL["time_counter"]
timeSeriesXESMF_BAL   = varXESMF_BAL.isel(X=slice(indexBox[0],indexBox[1]),Y=slice(indexBox[2],indexBox[3])).mean(dim=["X","Y"])

plt.figure(figsize=(15,7))
plt.plot(timeSiren      ,timeSeriesSiren    ,label="GLO (Siren)"        ,linestyle="solid"  ,color="k")
plt.plot(timeCDO_BAL    ,timeSeriesCDO_BAL  ,label="BAL-CMEMS (CDO)"    ,linestyle="solid"  ,color="b")
plt.plot(timeXESMF_BAL  ,timeSeriesXESMF_BAL,label="BAL-CMEMS (xESMF)"  ,linestyle="dashdot",color="g")
plt.plot(timeCDO_GLO    ,timeSeriesCDO_GLO  ,label="GLO-CMEMS (CDO)"    ,linestyle="dotted" ,color="0.3")
plt.grid()
plt.ylabel("box averaged [%s] - Z %5.2f m" % (varname,deptht.isel(Z=zindex).data), fontsize=15)
plt.legend()

figname = "box_average_%s_z_%s" % (varname,zindex)
plt.savefig(figname,bbox_inches="tight",dpi=600)

plt.show()

