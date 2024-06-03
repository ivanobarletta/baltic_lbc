import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

dfPath  = sys.argv[1]

# read dataframe with rPearson data
df      = pd.read_csv(dfPath,header=0,sep=",")

# filter out tg with missing coefficient
df      = df[df["tgOK"] == True]

print(df)

useproj = ccrs.EuroPP()
#useproj = ccrs.PlateCarree()
nonproj = ccrs.PlateCarree()

fig, ax = plt.subplots(1,1, subplot_kw={"projection":useproj},figsize=(7,12))

ax.coastlines(resolution="10m")
#ax.add_feature(cfeature.LAND,edgecolor=None)
#ax.add_feature(cfeature.OCEAN,edgecolor=None)
#ax.add_feature(cfeature.COASTLINE)


#box = [9,13.5,54,59]
box = [11,13.5,54,56]
ax.set_extent(box)

experiments = ["rP_CONTROL","rP_EXP_BT1","rP_EXP_BT2","rP_TESTRUN"]

def inBox(lon,lat,box):
    inBoxx = False
    # check if (lon,lat) point is in box
    if (lon > box[0]) and (lon < box[1]) and (lat > box[2]) and (lat < box[3]):
        inBoxx = True
    return inBoxx    



for nameTG,lon1,lat1,y1,y2,y3,y4 in zip(df["#ID"],df["lon"],df["lat"],df["rP_CONTROL"],df["rP_EXP_BT1"],df["rP_EXP_BT2"],df["rP_TESTRUN"]):

    if inBox(lon1,lat1,box) == False:
        print("skipping %s, not in selected box" % nameTG)
        continue

    # convert geographical coordinates to target projection
    xy = useproj.transform_point(lon1, lat1, nonproj)
    print("lon: %.2f ,lat: %.2f converted to %.4e,%.4e " % (lon1,lat1,xy[0],xy[1]))

    xvals=["a","b","c","d"]
    yvals=[y1,y2,y3,y4]
    fcolors=["b","c","r","g"]

    ax.scatter(lon1,lat1,color="k",marker="o",transform=ccrs.PlateCarree())

    ax_h = inset_axes(ax,
                      width=0.25,
                      height=0.6,
                      loc=3,
                      bbox_to_anchor=(xy[0],xy[1]),
                      bbox_transform=ax.transData,
                      borderpad=0,
                      axes_kwargs={'alpha': 0.35, 'visible': True})

    # plot unity bars
    for x,y,c in zip(xvals, yvals, fcolors):
        ax_h.bar(x, 1, label=str(x), fc="0.7")
    # plot bars
    for x,y,c in zip(xvals, yvals, fcolors):
        ax_h.bar(x, y, label=str(x), fc=c)

    ax_h.set_ylim(0.5,1)    
    ax_h.axis('off')
    #ax_h.grid()
    #ax_h.set_xticks([])
    ax_h.set_title(nameTG,fontsize=6,color="r")

plt.savefig("rP_map_barplot_20220101_20220401",dpi=600,bbox_inches="tight")
plt.show()


