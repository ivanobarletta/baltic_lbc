import pandas as pd
import xarray as xr 
import numpy as np
from scipy.spatial import cKDTree
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import os

pathListTG          = "target_inventory_Sensor1_202403.csv"
pathCoords          = "coordinates.nc"      # path to un-masked coordinates
pathMask            = "mask_gridT.nc"       # path to maskd file

outFile             = "inventoryTG_FilteredDistanceTime.csv"

distanceThreshold   = 1.5 * (1./36)           # 1.5 times the model resolution 
hugeValue           = 1e12

makePlot            = False

# read the csv as DataFrame
df  = pd.read_csv(pathListTG,header=0,sep=";")

nTG = df.shape[0]

mask2D  = xr.open_dataset(pathMask).isel(time_counter=0,deptht=0)["mask"]

# create a tree from model coordinates
print ("creating Tree from NEMO coords")
nav_lon = xr.open_dataset(pathCoords,decode_times=False)["nav_lon"]
nav_lat = xr.open_dataset(pathCoords,decode_times=False)["nav_lat"]
xy_NEMO = np.vstack((nav_lon.data.flatten(),nav_lat.data.flatten())).T
tree    = cKDTree(xy_NEMO)

# masked coordinates
print ("creating Tree from masked NEMO coords")
nav_lon_ma  = np.where(np.isnan(mask2D),hugeValue,nav_lon.data)
nav_lat_ma  = np.where(np.isnan(mask2D),hugeValue,nav_lat.data)
xy_NEMO_ma  = np.vstack((nav_lon_ma.flatten(),nav_lat_ma.flatten())).T
tree_ma     = cKDTree(xy_NEMO_ma)

# create query points from TG coordinates
xy_TG       = np.vstack((df["lon"],df["lat"])).T

# make query 
print ("making query from coordinates")
dist,idx        = tree.query(xy_TG)
print ("making query from masked coordinates")
dist_ma,idx_ma  = tree_ma.query(xy_TG)

for tg in range(nTG):
    print("TG coords: %6.3f, %6.3f" % (xy_TG[tg,0],xy_TG[tg,1]))
    print()
    print("Closest nav_lon,nav_lat:")
    print(nav_lon.data[np.unravel_index(idx[tg],nav_lon.shape)],nav_lat.data[np.unravel_index(idx[tg],nav_lon.shape)])
    print("Closest masked nav_lon,nav_lat:")
    print(nav_lon.data[np.unravel_index(idx_ma[tg],nav_lon.shape)],nav_lat.data[np.unravel_index(idx_ma[tg],nav_lon.shape)])

# store unravel indexes 
unravel_idx     = np.unravel_index(idx,nav_lon.shape)
unravel_idx_ma  = np.unravel_index(idx_ma,nav_lon.shape)

# make plot 
if makePlot:
    fig,ax = plt.subplots(subplot_kw={"projection":ccrs.PlateCarree()},layout="constrained")

    ax.scatter(df["lon"],df["lat"],s=5,transform=ccrs.PlateCarree(),color="b",marker="o");ax.coastlines(resolution="10m")   # plot TG scatter
    ax.scatter(nav_lon.data[unravel_idx],nav_lat.data[unravel_idx],marker="x",s=5,transform=ccrs.PlateCarree())             # closest unmasked coords
    ax.scatter(nav_lon.data[unravel_idx_ma],nav_lat.data[unravel_idx_ma],marker=">",s=5,transform=ccrs.PlateCarree())       # closest unmasked coords
    ax.pcolormesh(nav_lon,nav_lat,mask2D,alpha=0.2,transform=ccrs.PlateCarree())

    plt.show()

# subset TideGauge DataFrame by distance
df2 = df[dist_ma < distanceThreshold]

# subset by time
df3 = df2[ (df2["ini_year"] < 2023) & (df2["end_year"] > 2022) ]

# save to CSV
df3.to_csv(outFile)