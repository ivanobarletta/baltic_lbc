import xarray as xr
from eofs.xarray import Eof
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from dateutil.parser import parse
import pandas as pd
from scipy.stats import pearsonr

projection  = ccrs.EuroPP()
projGeo     = ccrs.PlateCarree()
props       = dict(boxstyle='round', facecolor='wheat', alpha=0.5)


pathOut     = "dsOut_4EOF_20220101_20230701.nc"

nModes      = 3             # the number of modes to analyze
vminn,vmaxx = -0.2,0.2
box         = [6,16,53,60]
date1,date2 = "2022-1-1","2023-1-1"
lResample   = True
resampleWindow  = "5D"

date1,date2 = parse(date1),parse(date2)

dsOut           = xr.open_dataset(pathOut)

nExperiments    = dsOut.sizes["nExp"]

varianceArray   = np.zeros((nExperiments+1,nModes))     # +1 because I have also TGs
cfArray         = []

# slice dataset
dsOut           = dsOut.sel(time=slice(date1,date2))
# resample 
if lResample:
    dsOut   = dsOut.resample(time=resampleWindow).mean()
    half_window = str(0.5*float(resampleWindow.replace("D","")))+"D"
    dsOut["time"] = dsOut.time + pd.Timedelta(half_window)

experimentNames = ["CONTROL","EXP_BT1","EXP_BT2","TESTRUN"]
experimentColors    = ["b","c","r","g"]

experimentNames = ["CONTROL","TESTRUN"]
experimentColors    = ["b","g"]


# create a solver for TG
solverTG        = Eof(dsOut["SLEV"],weights=None,center=True)

# create a list of solvers for models
#solversModel    = [Eof(model,weights=None,center=True) for model in dsOut["ssh"]]
print("computing solvers")
solversModel    = [Eof(dsOut["ssh"].isel(nExp=iexp),center=True) for iexp in range(nExperiments) ]

# gather variance for TG and models
for iMode in range(nModes):
    varianceArray[0,iMode] = 100*solverTG.varianceFraction().isel(mode=[iMode]).values[0]
    for iExp in range(nExperiments):
        varianceArray[1+iExp,iMode] = 100*solversModel[iExp].varianceFraction().isel(mode=[iMode]).values[0]



print("extracting pcs")
tg_pcs          = solverTG.pcs(npcs=nModes,pcscaling=1)                              # pcscaling = 1 -> PC / variance
models_pcs      = [solver.pcs(npcs=nModes,pcscaling=1) for solver in solversModel]

# calculate pearson correlation among pcs
rPearson_array = np.zeros((nExperiments,nModes,2))
for iExp in range(nExperiments):
    for iMode in range(nModes):
        output = pearsonr(tg_pcs.isel(mode=iMode),models_pcs[iExp].isel(mode=iMode))
        rPearson_array[iExp,iMode,0] = output.statistic
        rPearson_array[iExp,iMode,1] = output.pvalue

print("extracting eofs")
tg_eofs         = solverTG.eofs(neofs=nModes,eofscaling=2)                           # eofscaling = 2 -> EOF * variance
models_eofs     = [solver.eofs(neofs=nModes,eofscaling=2) for solver in solversModel]



print("making plot")
nExperiments    = dsOut.sizes["nExp"]
nCols           = 3
nRows           = 1+1+nExperiments

nSubplots       = nCols * nRows

axes            = []
fig,ax          = plt.subplots(nRows,nCols,layout="constrained",figsize=(10,12),subplot_kw={"projection":projection})
#fig,ax          = plt.subplots(nRows,nCols,figsize=(10,12),subplot_kw={"projection":projection})

# I remove the axes because I dont want projection for the first row of subplots
ax[0,0].remove()
ax[0,1].remove()
ax[0,2].remove()

ax[0,0] = fig.add_subplot(nRows,nCols,1)
ax[0,1] = fig.add_subplot(nRows,nCols,2)
ax[0,2] = fig.add_subplot(nRows,nCols,3)

ax[0,1].sharey(ax[0,0])
ax[0,2].sharey(ax[0,0])

ax[0,0].tick_params(axis='x', rotation=30)
ax[0,1].tick_params(axis='x', rotation=30)
ax[0,2].tick_params(axis='x', rotation=30)


for iMode in range(tg_pcs.sizes["mode"]):
    tg_pcs.isel(mode=iMode).plot(ax=ax[0,iMode],color="k",linewidth=0.6,label="TG")

ax[0,0].set_xlabel("")
ax[0,1].set_xlabel("")
ax[0,2].set_xlabel("")

# add Pearson
for iMode in range(nModes):
    for iExp in range(nExperiments):
        ax[0,iMode].text(0.5, 0.9-float(iExp)/10, "rP %s: %.2f" % (experimentNames[iExp],rPearson_array[iExp,iMode,0]), 
                        transform=ax[0,iMode].transAxes, fontsize=8,
                        verticalalignment='center', bbox=props, rotation=0)
            
for iExp,(name,color) in enumerate(zip(experimentNames,experimentColors)):
    print(iExp,name,color)
    for mode in range(models_pcs[iExp].sizes["mode"]):
        print(mode)
        models_pcs[iExp].isel(mode=mode).plot(ax=ax[0,mode],color=color,linewidth=0.6,label=name)

ax[0,1].legend(ncols=nExperiments+1,fontsize=8)


for mode in range(tg_eofs.sizes["mode"]):
    minval = vminn
    maxval = vmaxx    
    #tg_eofs.isel(mode=mode)
    if mode > 10:
        print("adjusting min/max for mode %d" % mode)
        minval = 0.5*vminn
        maxval = 0.5*vmaxx    
    cf = ax[1,mode].scatter(x=tg_eofs.lon,
                            y=tg_eofs.lat,
                            c=tg_eofs.isel(mode=mode).values,
                            transform=projGeo,
                            vmin=minval,
                            vmax=maxval,
                            cmap="bwr",
                            edgecolors="k",
                            linewidths=0.5)
    cfArray.append(cf)

ax[1,0].set_ylabel("EOF TG")

for iExp,name in enumerate(experimentNames):
    for mode in range(models_eofs[iExp].sizes["mode"]):
        minval = vminn
        maxval = vmaxx    
        factor = 1
        if mode > 10:
            print("adjusting min/max for mode %d" % mode)            
            factor  = 0.5
            minval = 0.5*vminn
            maxval = 0.5*vmaxx    
    
        cf = ax[iExp+2,mode].scatter(x=models_eofs[iExp].lon,
                                    y=models_eofs[iExp].lat,
                                    c=models_eofs[iExp].isel(mode=mode),
                                    transform=projGeo,
                                    vmin=minval,
                                    vmax=maxval,
                                    cmap="bwr",
                                    edgecolors="k",
                                    linewidths=0.5)
        cfArray.append(cf)

    ax[iExp+2,0].set_ylabel("EOF %s" % name)

ax[0,1].set_ylabel("")
ax[0,2].set_ylabel("")

for row in range(1,nRows):
    for col in range(nCols):
        ax[row,col].coastlines(resolution="10m")
        ax[row,col].set_extent(box)

#cb1 = fig.colorbar(cfArray[-3],ax=ax[-1,0],orientation="horizontal")
#cb2 = fig.colorbar(cfArray[-2],ax=ax[-1,1],orientation="horizontal")
#cb3 = fig.colorbar(cfArray[-1],ax=ax[-1,2],orientation="horizontal")

cb1 = fig.colorbar(cf,ax=ax[-1,0],orientation="horizontal")
cb2 = fig.colorbar(cf,ax=ax[-1,1],orientation="horizontal")
cb3 = fig.colorbar(cf,ax=ax[-1,2],orientation="horizontal")


cb1.set_label("m")
cb2.set_label("m")
cb3.set_label("m")

# add y labels (the default do not work! I don't know why!)
ax[1,0].text(-0.2, 0.5, "EOF - %s " % "TG", transform=ax[1,0].transAxes, fontsize=15,
            verticalalignment='center', bbox=props, rotation=90)

for iExp in range(nExperiments):
    ax[iExp+2,0].text(-0.2, 0.5, "EOF - %s " % experimentNames[iExp], transform=ax[iExp+2,0].transAxes, fontsize=15,
            verticalalignment='center', bbox=props, rotation=90)


# add variance fraction in Maps (TG)
for iMode in range(nModes):
    ax[1,iMode].text(0.7, 0.8, "%d%%" % varianceArray[0,iMode], transform=ax[1,iMode].transAxes, fontsize=15,
            verticalalignment='center', bbox=props, rotation=0)
# models
for iMode in range(nModes):
    for iExp in range(nExperiments):
        ax[2+iExp,iMode].text(0.7, 0.8, "%d%%" % varianceArray[1+iExp,iMode], transform=ax[2+iExp,iMode].transAxes, fontsize=15,
            verticalalignment='center', bbox=props, rotation=0)


str1    = date1.strftime("%Y%m%d")
str2    = date2.strftime("%Y%m%d")

plt.savefig("eofTGsummary_%s_%s" % (str1,str2),dpi=600,bbox_inches="tight")