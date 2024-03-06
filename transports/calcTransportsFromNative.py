import warnings
import numpy as np
import xarray as xr
import glob 
import os
import sys

"""
    -------------------------------------------------
    |                                               |
    |       transectType = we                       |    
    |                                               |      
    | yidx o--------------o                         |
    |                                               |
    |    xidx1           xidx2                      |    
    |                                               |      
    |                                               |
    |                                 o  yidx2      |
    |                                 |             |    
    |                                 |             |      
    |                                 |             |
    |                                 |             |
    |                                 |             |    
    |                                 |             |      
    |                                 o  yidx1      |
    |                                               |
    |                                xidx           |    
    |                                               |      
    |                                               |
    |                        transectType = sn      |
    |                                               |    
    |                                               |      
    |                                               |                              
    -------------------------------------------------

    transectType == "we" 
    idxList = [yidx,xidx1,xidx2]
    transectType == "sn" 
    idxList = [xidx,yidx1,yidx2]
    

"""

uMaskPath       = "../static_files/mask_gridU.nc"
vMaskPath       = "../static_files/mask_gridV.nc"
coordsPath      = "../target_grid/coordinates.nc"
zMeshPath       = "../static_files/mesh_zgr.nc"


# east2
rootFolder      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL"
rootFileName    = "NEATL36_1d25h-m_3DU-uo_*.nc"
fileList        = sorted(glob.glob(os.path.join(rootFolder,"*",rootFileName)))
transectType    = "sn"      # south-north ->  sum(u * e2u * e3u * umask )
idxList         = [-3,1478,1478+73]     # east2 open bonudary (type (sn))
outFile         = "transports_East2.nc"
Ntasks          = 50
varName         = "uo"

# Kattegat 
rootFolder      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL"
rootFileName    = "NEATL36_1d25h-m_3DV-vo_*.nc"
fileList        = sorted(glob.glob(os.path.join(rootFolder,"*",rootFileName)))
transectType    = "we"      # south-north ->  sum(u * e2u * e3u * umask )
idxList         = [1616,945,1033]       # Kattegat (type (we))
outFile         = "transports_Kattegat.nc"
Ntasks          = 50
varName         = "vo"

# Skagerrat
rootFolder      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL"
rootFileName    = "NEATL36_1d25h-m_3DU-uo_*.nc"
fileList        = sorted(glob.glob(os.path.join(rootFolder,"*",rootFileName)))
transectType    = "sn"      # south-north ->  sum(u * e2u * e3u * umask )
idxList         = [905,1592,1734]       # Skagerrak (type (sn))
outFile         = "transports_Skagerrat_compute.nc"
Ntasks          = 20
varName         = "uo"


def calcXFaces(coordsPath="",zCoordsPath="",vMaskPath="",idxList=[]):
    """
        Calculate areas of faces along x 

        areas = e1v * e3v_0 * vmask     [m**2]
                 dx *  dz   * mask

        return numpy array with shape (nz,nx)           

    """
    # do some checks
    if not os.path.isfile(coordsPath):
        raise Exception("file: %s not found" % coordsPath )
    
    if not os.path.isfile(zCoordsPath):
        raise Exception("file: %s not found" % zCoordsPath)
    
    if not os.path.isfile(vMaskPath):
        raise Exception("file: %s not found" % vMaskPath )

    dsCoords    = xr.open_dataset(coordsPath,decode_times=False)
    try:
        dsZCoords   = xr.open_dataset(zCoordsPath)
    except:    
        dsZCoords   = xr.open_dataset(zCoordsPath,decode_times=False) 
    try:    
        dsVMask     = xr.open_dataset(vMaskPath)
    except:        
        dsVMask     = xr.open_dataset(vMaskPath,decode_times=False)

    e1v = dsCoords["e1v"].isel(y=idxList[0],x=slice(idxList[1],idxList[2]))

    e3v = dsZCoords["e3v_0"].isel(y=idxList[0],x=slice(idxList[1],idxList[2]))

    # the variable might be umask or just "mask"
    try:
        vmask = dsVMask["mask"].isel(y=idxList[0],x=slice(idxList[1],idxList[2]))
    except:
        vmask = dsVMask["vmask"].isel(y=idxList[0],x=slice(idxList[1],idxList[2]))

    nz = dsZCoords.sizes["z"]

    e1v     = np.squeeze(e1v)                                           # shape (nx)             
    # repeat e2u along times along z
    e1v     = np.repeat(np.expand_dims(e1v,axis=0),repeats=nz,axis=0)   # shape (nz,nx)
    e3v     = np.squeeze(e3v)                                           # shape (nz,nx)
    vmask   = np.squeeze(vmask)                                         # shape (nz,nx)

    facesArea = e1v.data * e3v.data * vmask.data

    return facesArea

def calcYFaces(coordsPath="",zCoordsPath="",uMaskPath="",idxList=[]):
    """
        Calculate areas of faces along y 

        areas = e2u * e3u_0 * umask     [m**2]
                 dy *  dz   * mask

        return numpy array with shape (nz,ny)           

    """

    # do some checks
    if not os.path.isfile(coordsPath):
        raise Exception("file: %s not found" % coordsPath )
    
    if not os.path.isfile(zCoordsPath):
        raise Exception("file: %s not found" % zCoordsPath)
    
    if not os.path.isfile(uMaskPath):
        raise Exception("file: %s not found" % uMaskPath )

    dsCoords    = xr.open_dataset(coordsPath,decode_times=False)
    try:
        dsZCoords   = xr.open_dataset(zCoordsPath)
    except:    
        dsZCoords   = xr.open_dataset(zCoordsPath,decode_times=False) 
    try:    
        dsUMask     = xr.open_dataset(uMaskPath)
    except:        
        dsUMask     = xr.open_dataset(uMaskPath,decode_times=False)

    e2u = dsCoords["e2u"].isel(x=idxList[0],y=slice(idxList[1],idxList[2]))

    e3u = dsZCoords["e3u_0"].isel(x=idxList[0],y=slice(idxList[1],idxList[2]))

    # the variable might be umask or just "mask"
    try:
        umask = dsUMask["mask"].isel(x=idxList[0],y=slice(idxList[1],idxList[2]))
    except:
        umask = dsUMask["umask"].isel(x=idxList[0],y=slice(idxList[1],idxList[2]))

    nz = dsZCoords.sizes["z"]

    e2u     = np.squeeze(e2u)                                           # shape (ny)             
    # repeat e2u along times along z
    e2u     = np.repeat(np.expand_dims(e2u,axis=0),repeats=nz,axis=0)   # shape (nz,ny)
    e3u     = np.squeeze(e3u)                                           # shape (nz,ny)
    umask   = np.squeeze(umask)                                         # shape (nz,ny)

    facesArea = e2u.data * e3u.data * umask.data

    return facesArea

"""
def calcXTransport():
    pass

def calcYTransport():
    pass

def extractXSlice():
    pass

def extractYSlice():
    pass

    
def extractTime(fileList,timeVarname="time_counter"):
    try:
        ds =  xr.open_mfdataset(fileList)
        time = ds[timeVarname].data
        ds.close()
    except:
        raise Exception("Error: time variable not found in fileList ")    

    return time

"""

def chunk_bounds(Ntasks,task,number):
    """
        returns start,end index of each chunk
        of a number divided in Ntasks

        example for: 
            number = 733
            Ntasks = 8


                (i0, i1)
                ---------   
        task 0  (0, 91)
        task 1  (91, 182)
        task 2  (182, 273)
        task 3  (273, 364)
        task 4  (364, 455)
        task 5  (455, 546)
        task 6  (546, 637)
        task 7  (637, 733)

    """

    if Ntasks > number:
        raise Exception("Ntasks cannot be > number")

    int_div = number // Ntasks
    i0 = int_div * task
    i1 = int_div * (task+1)
    if task == Ntasks-1:
        i1 += number % Ntasks
    return i0,i1

def main(fileList = []                  # list of files to process
         ,varName = "u"                 # variable name for velocity component
         ,transectType = "sn"           # type of transect: sn (south-north) / we (west-east)
         ,idxList = []                  # indexes identifying the transect: sn -> [xidx,yidx1,yidx2] or we -> [yidx,xidx1,xidx2]
         ,coordsPath = ""               # path of dataset with horizontal scale factors 
         ,zCoordsPath = ""              # path of dataset with vertical scale factors
         ,maskPath = ""                 # path of file with mask
         ,Ntasks = 1                    # number of tasks to subdivide computations
         ,outFile = "out.nc"):          # output file name
    
    Nfiles = len(fileList)

    # do some checks
    if Nfiles == 0:
        raise Exception("Error: File list is empty. Check the path of files")        

    print("Processing %d files" % Nfiles )

    if transectType not in ["sn","we"]:
        raise Exception("Error: transectType not recognized")

    if len(idxList) != 3:
        raise Exception("Error: len(idxList) must be 3")        

    # calculate areas of faces
    if transectType == "sn":
        faces = calcYFaces(coordsPath=coordsPath,zCoordsPath=zCoordsPath,uMaskPath=maskPath,idxList=idxList)

    if transectType == "we":
        faces = calcXFaces(coordsPath=coordsPath,zCoordsPath=zCoordsPath,vMaskPath=maskPath,idxList=idxList)

    ds_list = []

    for task in range(Ntasks):

        i0,i1 = chunk_bounds(Ntasks=Ntasks,task=task,number=Nfiles)    

        fileListSubset = fileList[i0:i1]     

        print ("task: %s - processing %d files " % (task,i1-i0))

        ds = xr.open_mfdataset(fileListSubset)

        if transectType == "sn":
            # multiply by faces and sum 
            transport = (ds[varName].isel(x=idxList[0],y=slice(idxList[1],idxList[2])) * faces ).sum(dim=["depthu","y"])
            ds_list.append( transport.to_dataset().rename({varName:"transport"}) )

        if transectType == "we":    
            # multiply by faces and sum 
            transport = (ds[varName].isel(y=idxList[0],x=slice(idxList[1],idxList[2])) * faces ).sum(dim=["depthv","x"])
            ds_list.append( transport.to_dataset().rename({varName:"transport"}) )

        ds.close()

    # concatenate datasets from each task
    print("Concatenating datasets")
    dsOut = xr.concat(ds_list, dim="time_counter")
    dsOut["transport"] = dsOut["transport"].assign_attrs({"units":"m3/s","coordinates":"time_counter"})        

    # add transect coordinates to output dataset
    dsCoords = xr.open_dataset(coordsPath,decode_times=False)
    if transectType == "sn":
        nav_lon = dsCoords["nav_lon"].isel(x=idxList[0],y=slice(idxList[1],idxList[2]))
        nav_lat = dsCoords["nav_lat"].isel(x=idxList[0],y=slice(idxList[1],idxList[2]))
    if transectType == "we":
        nav_lon = dsCoords["nav_lon"].isel(y=idxList[0],x=slice(idxList[1],idxList[2]))
        nav_lat = dsCoords["nav_lat"].isel(y=idxList[0],x=slice(idxList[1],idxList[2]))

    dsOut["nav_lon"] = nav_lon
    dsOut["nav_lat"] = nav_lat

    # save dataset
    print("Saving Output")
    #dsOut.to_netcdf(outFile)
    dsOut.compute().to_netcdf(outFile)

if __name__ == "__main__":
    main(fileList = fileList
        ,varName = varName
        ,transectType = transectType
        ,idxList = idxList
        ,coordsPath= coordsPath
        ,zCoordsPath= zMeshPath
        ,maskPath= uMaskPath
        ,Ntasks= Ntasks
        ,outFile= outFile
        )

