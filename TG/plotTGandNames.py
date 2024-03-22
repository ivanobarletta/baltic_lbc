import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import os

inventoryPath   = "inventoryTG_Filtered.csv"

df = pd.read_csv(inventoryPath,header=0,sep=",")

# get  names
names   = df["#ID"] 

# make plot
fig,ax = plt.subplots(subplot_kw={"projection":ccrs.PlateCarree()},layout="constrained")

ax.scatter(df["lon"],df["lat"],s=5,transform=ccrs.PlateCarree(),color="b",marker="o");ax.coastlines(resolution="10m")   
for lon,lat,id in zip(df["lon"],df["lat"],names):
    ax.text(lon,lat,id,transform=ccrs.PlateCarree(),fontsize=5)        

ax.coastlines(resolution="10m")

plt.savefig("TG_locations_and_names",dpi=600,bbox_inches="tight")



