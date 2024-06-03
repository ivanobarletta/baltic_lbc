import cartopy.crs as ccrs
import pandas as pd 
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    raise Exception("You must provide path of df")

dfPath  = sys.argv[1]

df = pd.read_csv(dfPath,header=0,sep=",")

df2 = df.copy()
# set the ID as index
df2.index = df["#ID"]

# account only for OK tg
df2 = df2[df2["tgOK"] == True]

# sort by distance
df2 = df2.sort_values(by="tgDist",ascending=False)

# set color and names of experiments

experiments = ["rP_CONTROL","rP_EXP_BT1","rP_EXP_BT2"]
colors = ["b","c","r"]

experiments = ["rP_CONTROL","rP_EXP_BT1","rP_EXP_BT2","rP_TESTRUN"]
colors = ["b","c","r","g"]

# set type of plot ["bar" for vertical bars, "barh" for horizontal bars]
plot_type = "barh"

if plot_type == "bar":

    fig, ax = plt.subplots(2,1,figsize=(11,7))
    df2["tgDist"].to_frame().plot(ax=ax[0],kind=plot_type,color="k",width=0.2)
    df2[experiments].plot(ax=ax[1],kind=plot_type,color=colors)

    ax.set_ylim(0.5,1)
    ax.set_ylabel("rPearson")
    ax[1].grid()
else:
    #fig, ax = plt.subplots(1,2,figsize=(7,11),sharey=True,width_ratios=[3, 1])
    fig, ax = plt.subplots(1,3,figsize=(7,11),sharey=True,width_ratios=[3,1,1])

    df2[experiments].plot(ax=ax[0],kind=plot_type,color=colors,zorder=1)
    df2["tgDist"].to_frame().plot(ax=ax[1],kind=plot_type,color="k",width=0.2,zorder=1,legend=None)

    diff_r_bcl  = (df2["rP_TESTRUN"]-df2["rP_CONTROL"])
    diff_r_bt   = (df2["rP_EXP_BT2"]-df2["rP_EXP_BT1"])
    df_diff     = pd.DataFrame()
    df_diff["rP_diff_bcl"]  = pd.DataFrame(data=diff_r_bcl,columns=["rP_diff_bcl"])
    df_diff["rP_diff_bt"]   = pd.DataFrame(data=diff_r_bt ,columns=["rP_diff_bt"])

    df_diff[["rP_diff_bt","rP_diff_bcl"]].plot(ax=ax[2],kind=plot_type,color=colors[2:])

    ax[0].grid(zorder=0)
    ax[1].grid(zorder=0) 
    ax[0].set_xlim(0.5,1)
    ax[0].set_xlabel("rPearson")

    ax[1].set_xlabel("Dist from OBC [km]")

    ax[0].set_title("Correlation against Daily Averaged Tide Gauge Sea Level")
    ax[1].set_xlim(0,600)
    ax[2].set_xlim(-0.1,0.1)


plt.savefig("%s_plot3_%s" % (plot_type,dfPath.replace(".csv","")),bbox_inches="tight",dpi=600)


# make barplot on map

#proj = ccrs.Orthographic(central_longitude=12,central_latitude=55.5,globe=None)
#fig = plt.figure()
#ax = fig.add_axes([0,0,1,1],projection=proj)
#ax.coastlines(resolution="50m")
#ax.set_extent([10,14,54,57])

