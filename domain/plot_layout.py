import numpy as np 
import cartopy.crs as ccrs 
import xarray as xr 
import matplotlib.pyplot as plt 
import pandas as pd
from shapely.geometry.polygon import Polygon
import sys
from os.path import isfile

print("Warning! this script is under development!")
print("Might not give reasonable results with global meshes")

if len(sys.argv) < 3:
    raise("Error: you must provide path of layout.dat and path of coordinates.nc")

path_layout = sys.argv[1]
path_coords = sys.argv[2]

if isfile(path_layout) == False:
    raise("Error: path of layout.dat not existing")

if isfile(path_coords) == False:
    raise("Error: path of coordinates.nc not existing")

df = pd.read_csv(path_layout,header=2,sep="\s+",index_col=0)

nav_lon = xr.open_dataset(path_coords,decode_times=False)["nav_lon"]
nav_lat = xr.open_dataset(path_coords,decode_times=False)["nav_lat"]

# lower left corner x,y coordinate
llcx,llcy    = [],[]
# lower right corner x,y coordinate
lrcx,lrcy    = [],[]
# upper right corner x,y coordinate
urcx,urcy    = [],[]
# upper left corner x,y coordinate
ulcx,ulcy    = [],[]

for i,(nimpp,njmpp,nlci,nlcj,nldi,nldj) in enumerate(zip(df["nimpp"],df["njmpp"],df["nlci"],df["nlcj"],df["nldi"],df["nldj"])):

    # LowerLeftCorner (x,y)_idx
    llcx_idx    = nimpp - 1
    llcy_idx    = njmpp - 1 

    # LowerRightCorner (x,y)_idx
    lrcx_idx    = nimpp + nlci - nldi - 1
    lrcy_idx    = njmpp - 1 

    # UpperRightCorner (x,y)_idx
    urcx_idx    = nimpp + nlci - nldi - 1
    urcy_idx    = njmpp + nlcj - nldj - 1

    # UpperLeftCorner (x,y)_idx
    ulcx_idx    = nimpp - 1
    ulcy_idx    = njmpp + nlcj - nldj - 1

    llcx.append(nav_lon.isel(x=llcx_idx,y=llcy_idx).values)
    llcy.append(nav_lat.isel(x=llcx_idx,y=llcy_idx).values)

    lrcx.append(nav_lon.isel(x=lrcx_idx,y=lrcy_idx).values)
    lrcy.append(nav_lat.isel(x=lrcx_idx,y=lrcy_idx).values)

    urcx.append(nav_lon.isel(x=urcx_idx,y=urcy_idx).values)
    urcy.append(nav_lat.isel(x=urcx_idx,y=urcy_idx).values)

    ulcx.append(nav_lon.isel(x=ulcx_idx,y=ulcy_idx).values)    
    ulcy.append(nav_lat.isel(x=ulcx_idx,y=ulcy_idx).values) 


# make plot of sub-domains
proj    = ccrs.EuroPP()
proj    = ccrs.PlateCarree()

plt.figure(figsize=(12,12))
ax = plt.axes(projection=proj)
bounds = [(-30., 20., 25., 65.)]
ax.set_extent(*bounds, crs=ccrs.PlateCarree())
ax.coastlines(resolution='10m', linewidth=.5, color='black') # add map

for i,(x1,y1,x2,y2,x3,y3,x4,y4) in enumerate(zip(llcx,llcy,lrcx,lrcy,urcx,urcy,ulcx,ulcy)):
    #print(i,x1,y1,x2,y2,x3,y3,x4,y4)
    print("Area : %d" % i)
    print("Lower Left  Corner: (%.2f,%.2f)" % (x1,y1))
    print("Lower Right Corner: (%.2f,%.2f)" % (x2,y2))
    print("Upper Right Corner: (%.2f,%.2f)" % (x3,y3))
    print("Upper Left  Corner: (%.2f,%.2f)" % (x4,y4))

    #print(x1,y1,x2,y2,x3,y3,x4,y4)
    pgon = Polygon(((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x1,y1)))
    ax.add_geometries([pgon], crs=ccrs.PlateCarree(), facecolor=None, edgecolor='red', alpha=0.4)

    ax.text(x1,y1,str(i),transform=ccrs.PlateCarree(),fontsize=5,color="w")


plt.show()