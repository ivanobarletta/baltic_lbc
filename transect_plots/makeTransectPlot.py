import xarray as xr 
import matplotlib.pyplot as plt 
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

def make2x1plot(ds1,ds2,varname="varname",timeLevel=0,xidx=0,dsmask=None,vminn=0,vmaxx=10):
    """
    creates a 2x1 subplots from 2 NEMO LBC dataset with dimensions (T,Z,X,Y)
    
    -------------------------------------
    |                                   |
    |  ax[0]                            | depth
    |                                   |
    -------------------------------------
    -------------------------------------
    |                                   |
    |  ax[1]                            | depth
    |                                   |
    -------------------------------------
            latitude along transect
    """
    
    fig,ax = plt.subplots(2,1,sharex=True,sharey=True)
    date_str = pd.to_datetime(ds1["time_counter"].isel(T=timeLevel).data).strftime("%Y%m%d-%H")	
    var1 = ds1[varname].isel(T=timeLevel,X=xidx) 
    var2 = ds2[varname].isel(T=timeLevel,X=xidx) 
    if dsmask != None:
        maskzy = dsmask["mask"].isel(time_counter=0,x=-1)
        var1.data = var1.data * maskzy.data
        var2.data = var2.data * maskzy.data
    deptht = ds1["deptht"]
    nav_lat = ds1["nav_lat"].isel(X=xidx)
    xx,yy = np.meshgrid(nav_lat,deptht)   
    cmap = "jet" 
    if varname in ["vozocrtx","vomecrty"]:
        cmap = "RdBu"	
    ax[0].invert_yaxis()
    ax[0].set_ylim(60,0)	
    cf1 = ax[0].pcolormesh(xx,yy,var1,label="var1",cmap=cmap,vmin=vminn,vmax=vmaxx)
    cf2 = ax[1].pcolormesh(xx,yy,var2,label="var2",cmap=cmap,vmin=vminn,vmax=vmaxx)

    #cbar = fig.colorbar(cf1, ax=ax[:2], shrink=0.3, location='bottom',label = "%s [%s]" % (varname,var1.units))
    cbar = fig.colorbar(cf1, ax=ax[:2], shrink=0.3, location='right',label = "%s [%s]" % (varname,var1.units))

    ax[0].set_title(date_str)

    ax[0].set_ylabel("depth [m]")
    ax[1].set_ylabel("depth [m]")
    

    # place a text box in upper left in axes coords
    ax[0].text(0.05, 0.95, "CDO", transform=ax[0].transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

    ax[1].text(0.05, 0.95, "xESMF", transform=ax[1].transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

def make_1_2x1plot(dscoords=None,ds1=None,ds2=None,varname="varname",timeLevel=0,xidx=0,dsmask=None,vminn=0,vmaxx=10):
    """
    creates a 1 plus 2x1 subplots from 1 dataset with NEMO coordinates and 2 NEMO LBC datasets with dimensions (T,Z,X,Y)
    
    -------- -------------------------------------
    |      | |                                   |
    |      | |  ax[1]                            | depth
    |      | |                                   |
    |      | -------------------------------------
    | ax[0]| -------------------------------------
    |      | |                                   |
    |      | |  ax[2]                            | depth
    |      | |                                   |
    -------- -------------------------------------
                    latitude along transect
    """
    
    # create subplots layout
    fig = plt.figure(figsize=(13,6))
    gs  = fig.add_gridspec(2,2,width_ratios=[1,2.5],hspace=0.1, wspace=0.1)
    ax0 = fig.add_subplot(gs[:,0])
    ax1 = fig.add_subplot(gs[0,1])
    ax2 = fig.add_subplot(gs[1,1],sharex=ax1,sharey=ax1)

    ax = [ax0,ax1,ax2]

    nav_lon,nav_lat = dsCoords["nav_lon"],dsCoords["nav_lat"]

    # plot coordinates and transect line (LEFT axes)
    ax[0].scatter(nav_lon,nav_lat,s=0.5,color="k")
    ax[0].plot(nav_lon.isel(x=-1),nav_lat.isel(x=-1),color="b",linewidth=1)

    date_str = pd.to_datetime(ds1["time_counter"].isel(T=timeLevel).data).strftime("%Y%m%d-%H")	
    var1 = ds1[varname].isel(T=timeLevel,X=xidx) 
    var2 = ds2[varname].isel(T=timeLevel,X=xidx) 
    if dsmask != None:
        maskzy = dsmask["mask"].isel(time_counter=0,x=-1)
        var1.data = var1.data * maskzy.data
        var2.data = var2.data * maskzy.data
    deptht = ds1["deptht"]
    latTransect = ds1["nav_lat"].isel(X=xidx)
    xx,yy = np.meshgrid(latTransect,deptht)   

    # determine colormap
    cmap = "jet" 
    if varname in ["vozocrtx","vomecrty"]:
        cmap = "RdBu"	
    if varname == "vosaline":
        cmap = "hot_r"

    ax[1].invert_yaxis()
    ax[1].set_ylim(60,0)	

    # do transect plots
    cf1 = ax[1].pcolormesh(xx,yy,var1,label="var1",cmap=cmap,vmin=vminn,vmax=vmaxx)
    cf2 = ax[2].pcolormesh(xx,yy,var2,label="var2",cmap=cmap,vmin=vminn,vmax=vmaxx)

    cbar = fig.colorbar(cf1, ax=ax[1:], shrink=0.3, location='right', boundaries=np.linspace(vminn, vmaxx, 6))

    cbar.set_label(label = "mask * %s [%s]" % (varname,var1.units),size=15,weight='bold')

    ax[1].set_title(date_str)
    ax[1].yaxis.set_label_position("right")
    ax[1].set_ylabel("depth [m]")
    ax[2].yaxis.set_label_position("right")
    ax[2].set_ylabel("depth [m]")
    
    # eliminate xticks for upper twin axis
    #ax[1].set_xticklabels([]) this eliminates for box axes sharing x
    ax[1].get_xaxis().set_visible(False)

    # place a text box in upper left in axes coords
    ax[1].text(0.05, 0.95, "Bal (CDO)", transform=ax[1].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

    ax[2].text(0.05, 0.95, "Bal (xESMF)", transform=ax[2].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

    # save figure
    root    = "transect_%s_%s_coord_cdo_xesfm" 
    outfile = os.path.join(outFolder,root % (date_str,varname))
    plt.savefig(outfile,bbox_inches="tight",dpi=600)

def make_1_2x1plot2(dsCoords=None
                    ,ds1=None
                    ,ds2=None
                    ,varname="varname"
                    ,date_str="00000000"
                    ,xidx=0
                    ,dsmask=None
                    ,vminn=0
                    ,vmaxx=10
                    ,labels=[]):
    """
    creates a 1 plus 2x1 subplots from 1 dataset with NEMO coordinates and 2 NEMO LBC datasets with dimensions (T,Z,X,Y)
    
    -------- -------------------------------------
    |      | |                                   |
    |      | |  ax[1]                            | depth
    |      | |                                   |
    |      | -------------------------------------
    | ax[0]| -------------------------------------
    |      | |                                   |
    |      | |  ax[2]                            | depth
    |      | |                                   |
    -------- -------------------------------------
                    latitude along transect
    """

    if len(labels) != 2:
        raise Exception("you must provide 2 labels -> ( labels=[label1,label2] ) ")

    # create subplots layout
    fig = plt.figure(figsize=(13,6))
    gs  = fig.add_gridspec(2,2,width_ratios=[1,2.5],hspace=0.1, wspace=0.1)
    ax0 = fig.add_subplot(gs[:,0])
    ax1 = fig.add_subplot(gs[0,1])
    ax2 = fig.add_subplot(gs[1,1],sharex=ax1,sharey=ax1)

    ax = [ax0,ax1,ax2]

    nav_lon,nav_lat = dsCoords["nav_lon"],dsCoords["nav_lat"]

    # plot coordinates and transect line (LEFT axes)
    ax[0].scatter(nav_lon,nav_lat,s=0.5,color="k")
    ax[0].plot(nav_lon.isel(x=-1),nav_lat.isel(x=-1),color="b",linewidth=1)

    # determine time index to load from date 
    date2   = "%s-%s-%sT12" % (date_str[:4],date_str[4:6],date_str[6:])
    # indexes might possibly different from different datasets
    tidx_ds1 = np.argmin(np.abs(ds1["time_counter"].data - np.datetime64(date2)))
    tidx_ds2 = np.argmin(np.abs(ds2["time_counter"].data - np.datetime64(date2)))

    print("fetching [date,idx] [%s,%s] from ds1" % (ds1["time_counter"].isel(T=tidx_ds1).data, tidx_ds1 ))    
    print("fetching [date,idx] [%s,%s] from ds2" % (ds2["time_counter"].isel(T=tidx_ds2).data, tidx_ds2 ))

    date_diff = np.abs(ds2["time_counter"].isel(T=tidx_ds2).data - ds1["time_counter"].isel(T=tidx_ds1).data)

    if date_diff != 0:
        warnings.warn("Dates are different. This might possibly lead to inconsistent result!")        

    #date_str = pd.to_datetime(ds1["time_counter"].isel(T=timeLevel).data).strftime("%Y%m%d-%H")
    try:	
        var1 = ds1[varname].isel(T=tidx_ds1,X=xidx) 
        var2 = ds2[varname].isel(T=tidx_ds2,X=xidx) 
    except:
        raise Exception("problems with fetching the time index")

    if dsmask != None:
        maskzy = dsmask["mask"].isel(time_counter=0,x=-1)
        var1.data = var1.data * maskzy.data
        var2.data = var2.data * maskzy.data
    deptht = ds1["deptht"]
    latTransect = ds1["nav_lat"].isel(X=xidx)
    xx,yy = np.meshgrid(latTransect,deptht)   

    # determine colormap
    cmap = "jet" 
    if varname in ["vozocrtx","vomecrty"]:
        cmap = "RdBu"	
    if varname == "vosaline":
        cmap = "hot_r"

    ax[1].invert_yaxis()
    ax[1].set_ylim(60,0)	

    # do transect plots
    cf1 = ax[1].pcolormesh(xx,yy,var1,label="var1",cmap=cmap,vmin=vminn,vmax=vmaxx)
    cf2 = ax[2].pcolormesh(xx,yy,var2,label="var2",cmap=cmap,vmin=vminn,vmax=vmaxx)

    cbar = fig.colorbar(cf1, ax=ax[1:], shrink=0.3, location='right', boundaries=np.linspace(vminn, vmaxx, 6))

    cbar.set_label(label = "mask * %s [%s]" % (varname,var1.units),size=15,weight='bold')

    ax[1].set_title(date_str)
    ax[1].yaxis.set_label_position("right")
    ax[1].set_ylabel("depth [m]")
    ax[2].yaxis.set_label_position("right")
    ax[2].set_ylabel("depth [m]")
    
    # eliminate xticks for upper twin axis
    #ax[1].set_xticklabels([]) <-- this eliminates for both axes sharing x
    ax[1].get_xaxis().set_visible(False)

    # place a text box in upper left in axes coords
    ax[1].text(0.05, 0.95, labels[0], transform=ax[1].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

    ax[2].text(0.05, 0.95, labels[1], transform=ax[2].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

    # eliminate blank spaces
    str1    = labels[0].translate({ord(c): None for c in string.whitespace})
    str2    = labels[1].translate({ord(c): None for c in string.whitespace})

    # save figure
    root    = "2transect_%s_%s_coord_%s_%s" 
    outfile = os.path.join(outFolder,root % (date_str,varname,str1,str2))
    plt.savefig(outfile,bbox_inches="tight",dpi=600)

#make2x1plot(dsCDO,dsXESMF,varname="votemper",timeLevel=timeLevel,xidx=xidx,dsmask=dsTMASK,vminn=4,vmaxx=6)    

#make_1_2x1plot(dsCoords=dsCoords,ds1=dsCDO,ds2=dsXESMF,varname="votemper",timeLevel=timeLevel,xidx=0,dsmask=dsTMASK,vminn=4,vmaxx=6)
#make_1_2x1plot(dsCoords=dsCoords,ds1=dsCDO,ds2=dsXESMF,varname="vozocrtx",timeLevel=timeLevel,xidx=0,dsmask=dsTMASK,vminn=-0.3,vmaxx=0.3)

#make_1_2x1plot(dsCoords=dsCoords,ds1=dsCDO,ds2=dsXESMF,varname="votemper",timeLevel=180,xidx=0,dsmask=dsTMASK,vminn=8,vmaxx=18)
#make_1_2x1plot(dsCoords=dsCoords,ds1=dsCDO,ds2=dsXESMF,varname="vozocrtx",timeLevel=180,xidx=0,dsmask=dsTMASK,vminn=-0.3,vmaxx=0.3)

#make_1_2x1plot2(dsCoords=dsCoords,ds1=dsCDO,ds2=dsXESMF,varname="vozocrtx",date_str="20220102",xidx=0,dsmask=dsUMASK,vminn=-0.3,vmaxx=0.3,labels=["BAL (CDO)","BAL (xESMF)"])

#make_1_2x1plot2(dsCoords=dsCoords,ds1=dsCDO,ds2=dsCDOGlo,varname="vozocrtx",date_str="20220105",xidx=0,dsmask=dsUMASK,vminn=-0.3,vmaxx=0.3,labels=["BAL (CDO)","GLO (CDO)"])

#make_1_2x1plot2(dsCoords=dsCoords,ds1=dsSIREN,ds2=dsCDOGlo,varname="vozocrtx",date_str="20220105",xidx=0,dsmask=dsUMASK,vminn=-0.3,vmaxx=0.3,labels=["GLO (Siren)","GLO (CDO)"])

#make_1_2x1plot2(dsCoords=dsCoords,ds1=dsSIREN,ds2=dsCDOGlo,varname="vozocrtx",date_str="20220701",xidx=0,dsmask=dsUMASK,vminn=-0.3,vmaxx=0.3,labels=["GLO (Siren)","GLO (CDO)"])

make_1_2x1plot2(dsCoords=dsCoords,ds1=dsSIREN,ds2=dsCDOGlo,varname="vozocrtx",date_str="20221001",xidx=0,dsmask=dsUMASK,vminn=-0.3,vmaxx=0.3,labels=["GLO (Siren)","GLO (CDO)"])


plt.show()
