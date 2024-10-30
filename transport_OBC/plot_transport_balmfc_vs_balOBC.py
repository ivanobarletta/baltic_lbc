import matplotlib.pyplot as plt 
import matplotlib as mpl
import xarray as xr 
from datetime import datetime
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os

"""
    Note: the dataset I compare are from
        1)  transports calculated from BALMFC product (or a CMEMS
            product) in a "area". The transports has been computed 
            in all longitude transect in the selected area.

            the file is expected to contain the variables

                a) volume_transport2D(time,depth,latitude,longitude)
                b) xfaces(time,depth,latitude,longitude)    (the name xfaces is misleading,should be yfaces for transports along x)

            by selecting time and longitude, you get a view of the 
            trasport [in Sv] and the area [m2] of the single cells 
            of the transect.        

            ------------------------------------------------> lat
            |
            |  --.                                       .--
            |    |                                       |   
            |    ---.                                    |   
            |       |                              .-----. 
            |       |                              |
           \/       .-------.                .-----.   
                            |________________|    
           depth

        2)  transports calculated from OBC files obtained from the 
            same parent model. The file is structured in the same
            way, containing as many transect as in the relaxation
            zone of the corresponding nested model.      

            
"""

date1   = "2022-1-1"
date2   = "2023-12-31"

pathBALMFC  = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/transports_CMEMS/transports_balmfc_area_arkona.nc"
pathBALOBC  = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/transports_OBC/outputs/test_transport_NEATL36_east2_CDO_area.nc"
pathBALOBC  = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/transports_OBC/outputs/test_transport_NEATL36_east2_CDO_area_intlevel3d.nc"

dsBALMFC    = xr.open_dataset(pathBALMFC)
dsBALOBC    = xr.open_dataset(pathBALOBC)

# shift time ahead by 12h
dsBALMFC["time"] = dsBALMFC["time"] + pd.Timedelta("12h")

# clip datasets by dates
dsBALMFC    = dsBALMFC.sel(time=slice(date1,date2))
dsBALOBC    = dsBALOBC.sel(time=slice(date1,date2))

# select all X but last in dsBALOBC (is totally masked because of Umask)
dsBALOBC    = dsBALOBC.isel(X=slice(None,-1))


########################
#### make map plot #####
########################

fig = plt.figure(figsize=(5,7))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([12.75, 14.25, 53.75, 55.75], crs=ccrs.PlateCarree())

ax.coastlines(resolution="10m")

NY_balmfc = dsBALMFC.sizes["latitude"]
NX_balmfc = dsBALMFC.sizes["longitude"]

# plot transect lines of BALMFC 
for xidx in range(NX_balmfc):
    ax.plot(dsBALMFC["longitude"].isel(longitude=xidx).data*np.ones(NY_balmfc),dsBALMFC["latitude"],transform=ccrs.PlateCarree(),linestyle="dashed",color="k",linewidth=0.5)

# plot transect lines of BALOBC
NY_balobc = dsBALOBC.sizes["Y"]
NX_balobc = dsBALOBC.sizes["X"]

for xidx in range(NX_balobc):
    ax.plot(dsBALOBC["nav_lon"].isel(X=xidx),dsBALOBC["nav_lat"].isel(X=xidx),transform=ccrs.PlateCarree(),color="k",linewidth=0.5)

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)

mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False


#######################################
#### timeseries with shaded area ######
#######################################

balmfc_upper = dsBALMFC["volumeTransport_total"].mean(dim="longitude") + 3*dsBALMFC["volumeTransport_total"].std(dim="longitude")
balmfc_lower = dsBALMFC["volumeTransport_total"].mean(dim="longitude") - 3*dsBALMFC["volumeTransport_total"].std(dim="longitude")

balobc_upper = dsBALOBC["volume_transport_total"].mean(dim="X") + 3*dsBALOBC["volume_transport_total"].std(dim="X")
balobc_lower = dsBALOBC["volume_transport_total"].mean(dim="X") - 3*dsBALOBC["volume_transport_total"].std(dim="X")

plt.figure(figsize=(14,6))
plt.fill_between(dsBALMFC["time"],balmfc_upper,balmfc_lower,color="k",alpha=0.15,linestyle="dashed",label="balmfc")
dsBALMFC["volumeTransport_total"].mean(dim="longitude").plot(color="k",linewidth=0.3,linestyle="dashed")
plt.fill_between(dsBALOBC["time"],balobc_upper,balobc_lower,color="b",alpha=0.15,label="balobc")
dsBALOBC["volume_transport_total"].mean(dim="X").plot(color="b",linewidth=0.15)
plt.xlabel("")
plt.ylim(-0.3,0.3)
leg = plt.legend()
for line in leg.get_lines():
    line.set_linewidth(4.0)
plt.title("Volume Transports at Arkona")


##########################
#### all timeseries ######
##########################


plt.figure(figsize=(14,6))
for lon in range(dsBALMFC.sizes["longitude"]):
    dsBALMFC["volumeTransport_total"].isel(longitude=lon).plot(label="balmfc - lon: %d" % lon, linestyle="dashed",linewidth=0.3)
for X in range(dsBALOBC.sizes["X"]):
    dsBALOBC["volume_transport_total"].isel(X=X).plot(label="balobc X: %d" % X, linewidth=0.3)
plt.xlabel("")
plt.ylim(-0.3,0.3)
plt.legend(ncols=7)
plt.title("Volume Transports at Arkona")
