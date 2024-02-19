import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pandas as pd
import os
import string
import warnings

"""

    Note 1)
        Input files are NEMO LBC files like. Dimensions are
        (T,Z,X,Y) not (T,Z,Y,X) as more conventionally.
        The script might work also in the 2nd case but
        there's no guarantee of reasonable results...

    Note 2)
        static files must have same spatial dimensions of NEMO
        lbc files!!

        - coordinates,mask files: same x,y
        - mesh z files : same z,y,x

    Note 3)
        Be careful to the ordering of points in NEMO LBC files:

        The points along the normal direction to the edge are ordered
        growing towards the interior.

        in the case of east boundary (like east2 boundary of NEATL36 configuration)
        x is ordered in reverse

            longitude
       -------------->
       lon lon lon lon
    ..  3   2   1   0
       <--------------
            index

        The same applies to North LBC files accordingly. (Northmost latitude of
        NEMO LBC files has index y=0 (-- > .isel(y=0) )

        So, to make calculations involving static files (mask,scale factors), ordered in the conventional fashion, this
        must be accounted

        example:
            Multiplication by mask for the last east index.

            idx = 0     for NEMO LBC file
            idx = -1    for mask file

            tmasked = varLBC.isel(x=0) * mask.isel(x=-1)

        The same applies for scale factors (e2u,e2t,e3u_0..)

"""
indexBox        = [3,13,31,41]

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

outFolder   = "plots"

# static files
umaskPath   = "../static_files/mask_gridU_east2.nc"     # umask
tmaskPath   = "../static_files/mask_gridT_east2.nc"     # tmask
coordsPath  = "../static_files/coordinates_east2.nc"    # nav_lon,nav_lat,e2u

# datasets
cdoPath     = "../lbc_files_concatenated/NEATL36_east2_CDO_BAL_2022-2023.nc"
xesmfPath   = "../lbc_files_concatenated/NEATL36_east2_XESMF_BAL_2022-2023.nc"
sirenPath   = "../lbc_files_concatenated/NEATL36_obcdta_east2_2022-2023.nc"
cdoGLOPath  = "../lbc_files_concatenated/NEATL36_east2_CDO_GLO_2022-2023.nc"

dsCDO       = xr.open_dataset(cdoPath)
dsXESMF     = xr.open_dataset(xesmfPath)
dsSIREN     = xr.open_dataset(sirenPath)
dsCDOGlo    = xr.open_dataset(cdoGLOPath)

dsUMASK     = xr.open_dataset(umaskPath)
dsTMASK     = xr.open_dataset(tmaskPath)
dsCoords    = xr.open_dataset(coordsPath,decode_times=False)

timeLevel   = 100
xidx        = 0     # 0 - maximum longitude, -1 minimum longitude

varname     = "votemper"

def makeBoxAverageProfile(dscoords=None
                        ,ds_list=[]             # list of datasets to process
                        ,varname="varname"
                        ,date_str="00000000"
                        ,dsmask=None            # ds with 3D mask
                        ,dminn=0                # depth limits
                        ,dmaxx=10
                        ,label_list=[]          # list of labels marking the datasets (same order of ds_list!)
                        ,box=[]):               # box index list []
    
    """
    creates 1x2 subplots. Inputs are NEMO LBC files with dimensions (T,Z,X,Y) (not (T,Z,Y,X) as more conventionally)

        -------- ---------------
        |      | |              |
        |      | |              | 
        |      | |              |
        |      | |              |
     lat| ax[0]| |  ax[1]       | depth
        |      | |              |                     
        |      | |              | 
        |      | |              |
        -------- ----------------
        lon         [units] 
    """

    # some checks
    if len(ds_list) != len(label_list):
        raise("ds_list ad label_list must have same size")

    n_datasets = len(ds_list)

    # create subplots layout
    fig = plt.figure(figsize=(8,9))
    gs  = fig.add_gridspec(1,2,width_ratios=[1,2],hspace=0.1, wspace=0.2)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])

    lons,lats = dscoords["nav_lon"],dscoords["nav_lat"]
    lonsBox, latsBox = lons.isel(x=slice(box[0],box[1]),y=slice(box[2],box[3])), lats.isel(x=slice(box[0],box[1]),y=slice(box[2],box[3]))

    # plot coordinates and transect line (LEFT axes)
    def get_poly_corners(lons,lats):
        lons = lons.data
        lats = lats.data
        idx = ((0,0),(0,-1),(-1,-1),(-1,0),(0,0))
        x = np.asarray([lons[i,j] for i,j in idx])
        y = np.asarray([lats[i,j] for i,j in idx])

        return x,y

    xx,yy = get_poly_corners(lonsBox,latsBox)

    ax0.scatter(lons,lats,s=0.5,color="k")
    ax0.plot(xx,yy,color="r",linewidth=1.2)


    # determine time index to load from date
    date2   = "%s-%s-%sT12" % (date_str[:4],date_str[4:6],date_str[6:])
    # nearest search for time index.
    # indexes might possibly different from different datasets
    timeidx = []
    for ds in ds_list:
        tidx = np.argmin(np.abs(ds["time_counter"].data - np.datetime64(date2)))
        print("fetching [date,idx] [%s,%s] from ds" % (ds["time_counter"].isel(T=tidx).data, tidx ))    
        timeidx.append(tidx)
    
    print(timeidx)

    # retrieve variables 
    variables = []
    try:
        for ds,tidx in zip(ds_list,timeidx):
            variables.append(ds[varname].isel(T=tidx))
    except:
        raise Exception("problems with fetching the time index")

    # retrieve depth from 1st dataset
    deptht = ds_list[0]["deptht"]

    # retrieve mask, reorder shape and invert x axis

                                                       #swap x,y axes  # reverse along x         
    mask3D = np.swapaxes(dsmask["mask"].isel(time_counter=0).data,1,2)[:,::-1,:]

    # use the mask to mask variables
    for v in variables:
        v.data = v.data * mask3D

    for v in variables:
        print(v.shape)

    #plt.figure()
    #variables[0].isel(Z=13).plot(x="nav_lon",y="nav_lat")     
    #plt.show()

    # calculate box averaged profiles
    profiles = []

    #print(variables)

    for v in variables:
        profiles.append(v.isel(X=slice(box[0],box[1]),Y=slice(box[2],box[3])).mean(dim=["X","Y"]))

    styles = ["solid","dotted","dashdot","dashed"]

    #profiles[2].data = profiles[2].data + 0.2

    # plot profiles
    for prof,label,style in zip(profiles,label_list,styles):
        ax1.plot(prof.data,deptht.data,label=label,linestyle=style)

    ax1.invert_yaxis()
    ax1.set_ylim(dmaxx,dminn)
    ax1.grid(linestyle="-",linewidth=1)
    ax1.legend()
    ax1.set_title(date2)
    ax1.set_xlabel("%s [%s]" % (varname,variables[0].units),fontsize=15)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax1.set_ylabel("Depth[m]")

    root = "box_average_profiles_%s_%s" 
    outfile = os.path.join(outFolder,root % (date_str,varname))

    plt.savefig(outfile,bbox_inches="tight",dpi=600)


    
listOfDs    = [dsSIREN,dsCDO,dsXESMF,dsCDOGlo]
listOfLabels = ["GLO(Siren)","BAL(CDO)","BAL(xESMF","GLO(CDO)"]

makeBoxAverageProfile(dscoords=dsCoords,ds_list=listOfDs,varname="votemper",date_str="20220105",dsmask=dsTMASK,dminn=0,dmaxx=50,label_list=listOfLabels,box=indexBox)

makeBoxAverageProfile(dscoords=dsCoords,ds_list=listOfDs,varname="votemper",date_str="20220701",dsmask=dsTMASK,dminn=0,dmaxx=50,label_list=listOfLabels,box=indexBox)
makeBoxAverageProfile(dscoords=dsCoords,ds_list=listOfDs,varname="vosaline",date_str="20220701",dsmask=dsTMASK,dminn=0,dmaxx=50,label_list=listOfLabels,box=indexBox)

makeBoxAverageProfile(dscoords=dsCoords,ds_list=listOfDs,varname="votemper",date_str="20230701",dsmask=dsTMASK,dminn=0,dmaxx=50,label_list=listOfLabels,box=indexBox)
makeBoxAverageProfile(dscoords=dsCoords,ds_list=listOfDs,varname="vosaline",date_str="20230701",dsmask=dsTMASK,dminn=0,dmaxx=50,label_list=listOfLabels,box=indexBox)

