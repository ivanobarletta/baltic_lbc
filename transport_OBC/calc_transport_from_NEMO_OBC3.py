import xarray as xr 
import numpy as np
import glob
import os
import pandas as pd
from datetime import datetime
import re
import sys

# last update 2024-10-24 IVB   

"""
    Note 1) GENERAL usage 

        This script calculates transport across sections
        from "NEMO OBC like" files. 

        Files Necessary:

            - pathHGrid (containing horizontal scale factors)
            - pathZGrid (containing vertical scale factors)
            - pathMask  (containing masks)
        
            - rootPath (path where are files)
            - rootFileName (common root string for files) 
            
        the glob function builds the file list with wildcard    

    Note 2) INDEXING   

        The indexing of points for OBC files

        The indexing of point along normal direction
        to the boundary are reversed (0 outmost, -1 is towards inside )

        This is not, instead, for coordinates, mask files (indexes
        with conventional order)

        that is why xidx = 0 for LBC files corresponds to
                    nx-idx-1 in static files


    Note 3) STATIC FILES

        The STATIC files:

        3a) mask.nc (set pathMask)
            that contains (tmask,umask,vmask,fmask)
                
        3b) mesh_zgr.nc (set pathZGrid)
            that contains vertical scale factors (e3t_0,e3u_0,e3v_0..)

        3c) mesh_hgr.nc (set pathHGrid)
            that contains horizontal info on the mesh, like
            horizontal scale factors

            (e1u,e2u,e1t,e2t,e1v,e2v ...)

        You can obtain these files with SIREN create_meshmask.exe providing the 
            o ) bathy_meter.nc
            oo) coordinates.nc

        and selecting the option in_msh = 3 in SIREN namelist

        of corresponding simulation.       

        BUT! NOTE WELL!!

        In order to obtain these files consistently with the OBC domain,
        consider that OBC are NOT located necessarily on the border of 
        the model domain but they can be also located internally (usually 1 
        grid point inside). 

        That means that a good practice to obtain the above mesh files
        from the TOTAL domain (i.e. total bathy_meter.nc and coordinates.nc)
        and AFTER make the subset in the OBC region.

        If you apply SIREN create_meshmask.exe directly in the OBC domain you will
        get the borders with mask = 0, and it's not necessarily what you want. 

    Note 4)            

        l_rolling is for folders containing datasets coming from
        production systems with a date and production date.

            root_file_name_date_date_prod
            root_file_name_date_date_prod2
            root_file_name_date_date_prod3
            root_file_name_date_date_prod4

        In this case the file list contains duplicates that must
        be removed. Here the date with the latest production date
        is retained

        If the source folder contains duplicates l_rolling must be
        set to True
        
"""




# XESMF produced LBC files from CMEMS-BALtic
rootPath        = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/XESMF/outputs_bal"
rootFileName    = "NEATL36_east2_BAL_XESMF_????????.nc"
outFolder       = "outputs"
outFileName     = os.path.join(outFolder,"transport_NEATL36_east2_XESMF.nc")
l_rolling       = False
varname         = "vozocrtx"


# CDO produced LBC files from CMEMS-GLOBAL
rootPath        = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/CDO/outputs_glo"
rootFileName    = "NEATL36_east2_CDO_GLO_????????.nc"
outFolder       = "outputs"
outFileName     = os.path.join(outFolder,"transport_NEATL36_east2_CDO_GLO.nc")
varname         = "vozocrtx"
l_rolling       = False
computeS        = True         # compute salt transport
varnameS        = "vosaline"
parentModel     = "GLOMFC"
sourcePath      = rootPath


# NEATL36 Siren LBC files (best analysis)
rootPath        = "/mnt/netapp1/Store_Puertos/Store_IBIop/OPERATIVAS/STORE/GLOBAL/GLO12V4/BC/TMP/"
rootFileName    = "NEATL36_obcdta_east_2_*.nc"
outFolder       = "outputs"
outFileName     = os.path.join(outFolder,"transport_NEATL36_east2_SIREN_best_analysis.nc")
l_rolling       = True
varname         = "uo"
computeS        = True
varnameS        = "so"
parentModel     = "GLO12V4"
sourcePath      = rootPath



# CDO produced LBC files from CMEMS-BALtic
rootPath        = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/CDO/outputs_bal"
rootFileName    = "NEATL36_east2_BAL_CDO_202201??.nc"
rootFileName    = "NEATL36_east2_BAL_CDO_????????.nc"
l_rolling       = False
varnameU        = "vozocrtx"
computeS        = True         # compute salt transport
varnameS        = "vosaline"
parentModel     = "BALMFC"
sourcePath      = rootPath
rootPathStatic  = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/static_files/mask_testrun/create_meshmask/NEATL36_EAST2/"
pathZGrid       = rootPathStatic + "mesh_zgr.nc"                # I take e3u_0 from here
pathHGrid       = rootPathStatic + "mesh_hgr.nc"                # I take e2u from here
pathMask        = rootPathStatic + "mask.nc"                    # I take umask from here
xidx            = 0 # outmost index is 0
outFolder       = "outputs"
outFileName     = os.path.join(outFolder,"test_transport_NEATL36_east2_CDO_xidx_%s.nc" % str(xidx).zfill(2))
sshVarName      = "sossheig"
computeHeat     = True          # compute heat transport
varnameT        = "votemper"
obcType         = "E"

# GLOBAL NEATL36 Siren LBC files
rootPath        = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/GLOBAL/PSY4V3R1/BC/TMP/"
rootFileName    = "NEATL36_obcdta_east_2_*.nc"
l_rolling       = False
varnameU        = "vozocrtx"
computeS        = True         # compute salt transport
varnameS        = "vosaline"
parentModel     = "GLOBAL_PSY4V3R1"
sourcePath      = rootPath
rootPathStatic  = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/static_files/mask_controlrun/create_meshmask/NEATL36_EAST2/"
pathZGrid       = rootPathStatic + "mesh_zgr.nc"                # I take e3u_0 from here
pathHGrid       = rootPathStatic + "mesh_hgr.nc"                # I take e2u from here
pathMask        = rootPathStatic + "mask.nc"                    # I take umask from here
xidx            = 1 # outmost index is 0
outFolder       = "outputs"
outFileName     = os.path.join(outFolder,"test_transport_NEATL36_east2_GLOBAL_PSY4V3R1_SIREN_xidx_%s.nc" % str(xidx).zfill(2))
sshVarName      = "sossheig"
computeHeat     = True          # compute heat transport
varnameT        = "votemper"
obcType         = "E"


def calcYFaces(dsIn,sshVarName="sossheig",pathZGrid="",pathHGrid="",pathMask="",xidx=xidx,obcType="E"):
    """
        returns the areas in m**2 of the faces along a transect

        indexing of data (ONLY for data contained in OBC files!)
        if obcType is E
        xidx ---> ... 10 9 8 7 6 5 4 3 2 1 0 
        if obcType is W
        xidx ---> ... 0 1 2 3 4 5 6 7 8 9 10 
        
        for static data numbering is always the "natural" one

        I do the transformation to access to the STATIC data
        if obcType = E
            yfaxes = umask[nx-xidx-1] * e3u[nx-xidx-1] * e2u[nx-xidx-1]
    """

    if obcType not in ["W","E"]:
        raise Exception("Error: obcType (%s) not considered or not implemented" % obcType)

    # I need this to 
    ssh = dsIn[sshVarName]
    nt  = ssh.shape[0]      # 1st dim must be the time

    # pathZGrid must contain e3u_0
    # pathHGrid must contain e2u
    # pathUMask must contain mask
    e3Name  = "e3u_0"
    e2Name  = "e2u"

    dsZGrid = xr.open_dataset(pathZGrid)
    dsHGrid = xr.open_dataset(pathHGrid)
    dsMask  = xr.open_dataset(pathMask)

    # mesh_zgr.nc file (I found it in some paths in cesga)         
    nz = dsZGrid.sizes["Z"]
    nx = dsZGrid.sizes["X"]

    # cells thickness (using value at T point, should be okay)
    if obcType == "E":
        # dz
        dz = dsZGrid[e3Name].isel(X=nx-xidx-1)
        # ds
        dy = dsHGrid[e2Name].isel(X=nx-xidx-1)
    else:
        raise Exception("Error: obcType not implemented yet")
         
    # repeat ds along z (first I do isel and then I create z dim again)
    # xrray does not have "repeat" option like in numpy 
    dy = dy.expand_dims(dim={"Z":nz},axis=0)        # (Z,Y)
    # get u-mask
    umask = dsMask["umask"].isel(X=nx-xidx-1)       # (Z,Y) 
    
    # convert to float
    umask = umask.astype(float)
    # fill 0 with Nan
    umask.data[umask.data==0] = np.nan
    
    # extract ssh (to calculate jacobian)
    try: 
        ssh = ssh.isel(x=xidx)  # (t,y)
    except:
        ssh = ssh.isel(X=xidx)  # (t,y)    

    # calculate depth H 
    H   = np.nansum(dz.data * umask.data, axis=0)                   # (y)
    H   = np.repeat(np.expand_dims(H,axis=0),axis=0,repeats=nt)     # (t,y)

    Jacobian = (H+ssh.data) / H                                     # (t,y) 

    Jacobian = np.repeat(np.expand_dims(Jacobian,axis=1),axis=1,repeats=nz)     #(t,z,y)

                      #  (t,z,y)               (z,y)    (z,y)     (z,y)
    faces_data = Jacobian.compute().data * dz.data * dy.data * umask.data       #(t,z,y)
        
    #print(faces.shape)
    #area = np.nansum(faces)
    #print("area [m**2]: ",area)

    return faces_data

def loadDatasets(rootPath="rootPath",rootFileName="rootFileName"):
    print("Loading from folder:")
    print(rootPath)

    fileList = sorted(glob.glob(os.path.join(rootPath,rootFileName)))
    if len(fileList) == 0:
        raise Exception("loadDatasets: No Files Found!")

    # load mfdataset
    ds = xr.open_mfdataset(fileList,combine="nested",concat_dim="T")

    return ds

def loadDatasetsNoDuplicates(rootPath="rootPath",rootFileName="rootFileName"):
    """
        this function removes duplicates of dates

        there can be more than one file with referring to the 
        same date but with different production date.

        this function removes duplicates from the file list
        keeping the most recent production date

    """
    print("Loading from folder:")
    print(rootPath)

    fileList = sorted(glob.glob(os.path.join(rootPath,rootFileName)))
    if len(fileList) == 0:
        raise Exception("loadDatasets: No Files Found!")

    dates = []
    dates_prod = []
    # search for date and date_prod in file path string

    for file in fileList:
        match1 = re.search(r"\d{8}",file)
        match2 = re.search(r"R\d{8}",file)
        #print(file)
        #print("match1",match1)
        #print("match2",match2)
        date = datetime.strptime(match1.group(),"%Y%m%d")
        date_prod = datetime.strptime(match2.group()[1:],"%Y%m%d")
        dates.append(date)
        dates_prod.append(date_prod)

    nFiles = len(dates)
    dates = np.asarray(dates)
    dates_prod = np.asarray(dates_prod)
    
    if dates.shape != dates_prod.shape:
        raise Exception("problem with extraction of dates")

    # create Pandas dataframe with both dates
    df = pd.DataFrame(data = np.vstack((dates,dates_prod)).T, index=np.arange(nFiles),columns=["date","date_prod"])

    # dataFrame with dropped duplicates
    df2 = df.sort_values("date_prod").drop_duplicates("date",keep="last").sort_values("date")

    # subset fileList 
    fileListSubset = np.asarray(fileList)[df2.index].tolist()

    #f = open("subset.dat","w")
    #for file in fileListSubset:
    #    f.write(file+"\n")
    #f.close()


    # load mfdataset
    ds = xr.open_mfdataset(fileListSubset,combine="nested",concat_dim="T")

    #print(ds.compute())
    #print()
    #for dd in ds.compute().time_counter:
    #    print(dd)

    return ds

def calcXTransport(ds,facesArea,xidx=xidx,varnameU="uo"):

    volumeConversion    = 1e-6  # m^3/s -> Sv
    # vozo (T,Z,X=xidx,Y) -> (T,Z,Y) 
    vozo = ds[varnameU].isel(X=xidx)

    factor  = volumeConversion

    # faces has also shape (T,Z,Y)
    #              u                 dA        [m/s]*[m**2] 
    #           (T,Z,Y)            (T,Z,Y)
    prod =  vozo.compute().data * facesArea    # is a numpy array

    transport_direction1    = factor * np.nansum(prod,axis=(1,2),where=np.array(prod>=0))
    transport_direction2    = factor * np.nansum(prod,axis=(1,2),where=np.array(prod<0))
    transport_total         = factor * np.nansum(prod,axis=(1,2))

    return transport_direction1,transport_direction2,transport_total,factor*prod

def calcXTransportS(ds,facesArea,xidx=xidx,varnameU="uo",varnameS="so",obcType="E"):

    rhoWater            = 1025  # kg/m^3
    saltConversion      = rhoWater * 1e-12  # [g/Kg]*[m^3/s] -> 10^9 Kg/s

    factor  = 1.
    if varnameS == "vosaline":
        factor = saltConversion     

    # vozo (T,Z,X=xidx,Y) -> (T,Z,Y) 
    vozo = ds[varnameU].isel(X=xidx)
    if obcType == "E":
        voso = 0.5 * (ds[varnameS].isel(X=xidx) + ds[varnameS].isel(X=xidx-1))
    else:
        voso = 0.5 * (ds[varnameS].isel(X=xidx) + ds[varnameS].isel(X=xidx+1))

    # faces has also shape (T,Z,Y)
    #             S                      u                   dA           # [g/Kg]*[m/s]*[m**2] 
    #           (T,Z,Y)               (T,Z,Y)             (T,Z,Y)            
    prod =  voso.compute().data * vozo.compute().data * facesArea       # is a numpy array

    transport_direction1    = factor * np.nansum(prod,axis=(1,2),where=np.array(prod>=0))
    transport_direction2    = factor * np.nansum(prod,axis=(1,2),where=np.array(prod<0))
    transport_total         = factor * np.nansum(prod,axis=(1,2))

    return transport_direction1,transport_direction2,transport_total,factor*prod 

def calcXTransportH(ds,facesArea,xidx=xidx,varnameU="uo",varnameT="thetao",obcType="E"):

    rhoWater            = 1025  # kg/m^3
    cp                  = 3996  # J/Kg/K
    tref                = 0 
    heatConversion      = 1e-15  # units of 10^15 Watts

    factor  = 1.
    if varnameS == "vosaline":
        factor = heatConversion     

    # vozo (T,Z,X=xidx,Y) -> (T,Z,Y) 
    vozo = ds[varnameU].isel(X=xidx)
    if obcType == "E":
        voto = 0.5 * (ds[varnameT].isel(X=xidx) + ds[varnameT].isel(X=xidx-1))
    else:    
        voto = 0.5 * (ds[varnameT].isel(X=xidx) + ds[varnameT].isel(X=xidx+1))

    # faces has also shape (T,Z,Y)
    #                           T                          u                  dA           # [J/K/Kg]*[Kg/m**3] [K]*[m/s]*[m**2] 
    #                      (T,Z,Y)                     (T,Z,Y)             (T,Z,Y)            
    prod = cp * rhoWater * (voto.compute().data - tref) * vozo.compute().data * facesArea       # is a numpy array

    transport_direction1    = factor * np.nansum(prod,axis=(1,2),where=np.array(prod>=0))
    transport_direction2    = factor * np.nansum(prod,axis=(1,2),where=np.array(prod<0))
    transport_total         = factor * np.nansum(prod,axis=(1,2))

    return transport_direction1,transport_direction2,transport_total,factor*prod 

def saveDataset(ds,
                transports,             # it's a tuple (direction1,direction2,total,3D)
                facesArea,                  # numpy array
                outFileName="outname.nc",
                computeS=False,
                computeHeat=False,
                transportsS=None,       # it's a tuple (direction1,direction2,total,3D)
                transportsH=None,       
                parentModel="",
                sourcePath="",
                xidx=xidx):

    time    = ds["time_counter"].compute()
    depth   = ds["deptht"].isel(T=0).values
    nav_lat = ds["nav_lat"].isel(X=xidx).values
    nav_lon = ds["nav_lon"].isel(X=xidx).values

    daOutDirection1 = xr.DataArray(data=transports[0],dims=["time"],coords=dict(time=pd.DatetimeIndex(time)))
    daOutDirection1 = daOutDirection1.assign_attrs({"long_name":"volume_transport_direction1","units":"Sv"})

    daOutDirection2 = xr.DataArray(data=transports[1],dims=["time"],coords=dict(time=pd.DatetimeIndex(time)))
    daOutDirection2 = daOutDirection2.assign_attrs({"long_name":"volume_transport_direction2","units":"Sv"})

    daOutTotal = xr.DataArray(data=transports[2],dims=["time"],coords=dict(time=pd.DatetimeIndex(time)))
    daOutTotal = daOutTotal.assign_attrs({"long_name":"volume_transport_total","units":"Sv"})

    daOut3D     = xr.DataArray(data=transports[3],dims=["time","depth","lat"],coords=dict(time=pd.DatetimeIndex(time),depth=depth,lat=nav_lat))
    daOut3D     = daOut3D.assign_attrs({"long_name":"volume_transport","units":"Sv"})

    daFaces     = xr.DataArray(data=facesArea, dims=["time","depth","lat"],coords=dict(time=pd.DatetimeIndex(time),depth=depth,lat=nav_lat))
    daFaces     = daFaces.assign_attrs({"name":"faces_areas","long_name":"faces_areas_along_y","units":"m2"})

    # create and save dataset
    dsOut = xr.Dataset()

    dsOut["volume_transport_direction1"]    = daOutDirection1
    dsOut["volume_transport_direction2"]    = daOutDirection2
    dsOut["volume_transport_total"]         = daOutTotal
    dsOut["volume_transport"]               = daOut3D

    dAnav_lat = xr.DataArray(data = nav_lat, dims=["lat"])
    dAnav_lon = xr.DataArray(data = nav_lon, dims=["lat"])

    dsOut["nav_lat"] = dAnav_lat
    dsOut["nav_lon"] = dAnav_lon

    dsOut   = dsOut.assign_attrs({"parent_model":parentModel})
    dsOut   = dsOut.assign_attrs({"source_path":sourcePath})
    dsOut   = dsOut.assign_attrs({"obcType":obcType})    
    dsOut   = dsOut.assign_attrs({"xidx":str(xidx)})

    # store faces areas
    dsOut["faces"]  = daFaces

    if computeS:
        # check if transportS is not None
        if type(transportsS) == None:
            raise Exception("Error: saveDataset: transportS is None")
        
        daOutSDirection1 = xr.DataArray(data = transportsS[0],dims=["time"],coords=dict(time=pd.DatetimeIndex(time)))
        daOutSDirection1 = daOutSDirection1.assign_attrs({"long_name":"salt_transport_direction1","units":"10^9 Kg/s"})

        daOutSDirection2 = xr.DataArray(data = transportsS[1],dims=["time"],coords=dict(time=pd.DatetimeIndex(time)))
        daOutSDirection2 = daOutSDirection2.assign_attrs({"long_name":"salt_transport_direction2","units":"10^9 Kg/s"})

        daOutSTotal = xr.DataArray(data = transportsS[2],dims=["time"],coords=dict(time=pd.DatetimeIndex(time)))
        daOutSTotal = daOutSTotal.assign_attrs({"long_name":"salt_transport_total","units":"10^9 Kg/s"})

        daOutS3D    = xr.DataArray(data = transportsS[3],dims=["time","depth","lat"],coords=dict(time=pd.DatetimeIndex(time),depth=depth,lat=nav_lat))
        daOutS3D    = daOutS3D.assign_attrs({"long_name":"salt_transport","units":"10^9 Kg/s"})

        dsOut["salt_transport_direction1"]  = daOutSDirection1
        dsOut["salt_transport_direction2"]  = daOutSDirection2
        dsOut["salt_transport_total"]       = daOutSTotal
        dsOut["salt_transport"]             = daOutS3D

    if computeHeat:
        # check if transportS is not None
        if type(transportsH) == None:
            raise Exception("Error: saveDataset: transportH is None")
        
        daOutHDirection1 = xr.DataArray(data = transportsH[0],dims=["time"],coords=dict(time=pd.DatetimeIndex(time)))
        daOutHDirection1 = daOutHDirection1.assign_attrs({"long_name":"heat_transport_direction1","units":"10^15 Watt"})

        daOutHDirection2 = xr.DataArray(data = transportsH[1],dims=["time"],coords=dict(time=pd.DatetimeIndex(time)))
        daOutHDirection2 = daOutHDirection2.assign_attrs({"long_name":"heat_transport_direction2","units":"10^15 Watt"})

        daOutHTotal = xr.DataArray(data = transportsH[2],dims=["time"],coords=dict(time=pd.DatetimeIndex(time)))
        daOutHTotal = daOutHTotal.assign_attrs({"long_name":"heat_transport_total","units":"10^15 Watt"})

        daOutH3D    = xr.DataArray(data = transportsH[3],dims=["time","depth","lat"],coords=dict(time=pd.DatetimeIndex(time),depth=depth,lat=nav_lat))
        daOutH3D    = daOutH3D.assign_attrs({"long_name":"heat_transport","units":"10^15 Watt"})

        dsOut["heat_transport_direction1"]  = daOutHDirection1
        dsOut["heat_transport_direction2"]  = daOutHDirection2
        dsOut["heat_transport_total"]       = daOutHTotal
        dsOut["heat_transport"]             = daOutH3D

    dsOut.to_netcdf(outFileName,encoding={"time" : {"dtype": "i4"} },)


def main(l_rolling = False, computeS=False):

    if l_rolling:
        print("Loading Dataset with No Duplicates...")
        ds = loadDatasetsNoDuplicates(rootPath=rootPath,rootFileName=rootFileName)
    else:
        print("Loading Dataset...")
        ds = loadDatasets(rootPath=rootPath,rootFileName=rootFileName)

    print("Calculating Face Areas of transect..")
    facesArea   = calcYFaces(dsIn=ds,sshVarName=sshVarName,pathZGrid=pathZGrid,pathHGrid=pathHGrid,pathMask=pathMask,xidx=xidx)  # returns a np.array
    #facesT  = calcYFaces(dsIn=ds,sshVarName=sshVarName,pathZGrid=pathZGrid,pathHGrid=pathHGrid,pathMask=pathTMask,varType="s",xidx=xidx)  # returns a np.array

    print("nansum faces", np.nansum(facesArea))
    #print("nansum facesT", np.nansum(facesT))

    print("Calculating Transport..")
    # transports is a tuple! transports = (transport_direction1,transport_direction2,transport_total)
    transports   = calcXTransport(ds,facesArea,xidx=xidx,varnameU=varnameU)

    #print("transports")
    #print(transports)
    #print(transports)

    if computeS:
        # transportsS is a tuple     
        transportsS  = calcXTransportS(ds,facesArea=facesArea,xidx=xidx,varnameU=varnameU,varnameS=varnameS,obcType=obcType)

    if computeHeat:
        # transportsS is a tuple     
        transportsH  = calcXTransportH(ds,facesArea=facesArea,xidx=xidx,varnameU=varnameU,varnameT=varnameT,obcType=obcType)


    saveDataset(ds,                         # input
                transports,                 # it's a tuple (direction1,direction2,total)
                facesArea,
                outFileName=outFileName,    
                computeS=computeS,          # logical
                transportsS=transportsS,    # it's a tuple (direction1,direction2,total)
                computeHeat=computeHeat,
                transportsH=transportsH,
                parentModel=parentModel,    
                sourcePath=sourcePath,
                xidx=xidx)


if __name__ == "__main__":
    main(l_rolling=l_rolling,
         computeS=computeS)
