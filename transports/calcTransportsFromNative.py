import warnings
import numpy as np
import xarray as xr
import glob 
import os
import sys



"""
    The script calculates transport from NEMO outputs only along
    fixed curvilinear direction (see below):

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

    Two types of transect are considered:

    transectType == "we" 
        Fixed y index (latitude) -> West-East Transect 
        idxList = [yidx,xidx1,xidx2]
    transectType == "sn" 
        Fixed x index (longitude) -> South-North Transect
        idxList = [xidx,yidx1,yidx2]
    
"""

uMaskPath       = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/static_files/mask_testrun/mask_gridU.nc"
vMaskPath       = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/static_files/mask_testrun/mask_gridV.nc"
coordsPath      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/target_grid/coordinates.nc"
zCoordsPath     = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/static_files/mask_testrun/NEATL36_TESTRUN_1d25h-m_e3u_202201.nc"

########################### Gibraltar
rootFolder      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL"
rootFileName    = "NEATL36_1d25h-m_3DU-uo_*.nc"
fileList        = sorted(glob.glob(os.path.join(rootFolder,"*",rootFileName)))
transectType    = "sn"      # south-north ->  sum(u * e2u * e3u * umask )
idxList         = [515,400,425]       # Gibraltar (type (sn))
outFile         = "transports_control_run_Gibraltar.nc"
Ntasks          = 1
varName         = "uo"

########################### Kattegat 
rootFolder      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL"
rootFileName    = "NEATL36_1d25h-m_3DV-vo_*.nc"
fileList        = sorted(glob.glob(os.path.join(rootFolder,"*",rootFileName)))
transectType    = "we"      # south-north ->  sum(u * e2u * e3u * umask )
idxList         = [1616,945,1033]       # Kattegat (type (we))
outFile         = "transports_control_run_Kattegat.nc"
Ntasks          = 1
varName         = "vo"

########################### Skagerrat
rootFolder      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL"
rootFileName    = "NEATL36_1d25h-m_3DU-uo_*.nc"
fileList        = sorted(glob.glob(os.path.join(rootFolder,"*",rootFileName)))
transectType    = "sn"      # south-north ->  sum(u * e2u * e3u * umask )
idxList         = [905,1592,1734]       # Skagerrak (type (sn))
outFile         = "transports_control_run_Skagerrak.nc"
Ntasks          = 1
varName         = "uo"

########################### east2
rootFolder      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh"
rootFileName    = "NEATL36_TESTRUN_1d25h-m_3DU-uo_*.nc_ZNB"
#fileList        = sorted(glob.glob(os.path.join(rootFolder,"*",rootFileName)))
fileList        = sorted(glob.glob(os.path.join(rootFolder,rootFileName)))

uMaskPath       = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh/static/mask_gridU_ZNB.nc"
vMaskPath       = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh/static/mask_gridV_ZNB.nc"
coordsPath      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh/static/coordinates_NEATL36_ZNB.nc"
zCoordsPath     = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh/static/mesh_mask_testrun2_ZNB.nc"


transectType    = "sn"      # south-north ->  sum(u * e2u * e3u * umask )
idxList         = [-3,177,177+73]     # east2 open bonudary (type (sn))
outFile         = "transports_testrun2_run_East2.nc"
Ntasks          = 1
varName         = "uo"
simulationName  = "TESTRUN2"

calcSaltTransp  = True                                                                  # calculate Salt Transport
rootFileNameS   = "NEATL36_TESTRUN_1d25h-m_3DT-so_*.nc_ZNB"                             # root for salinity filenames
fileListS       =  sorted(glob.glob(os.path.join(rootFolder,rootFileNameS)))            # list of salinity files
varNameS        = "so"

calcHeatTransp  = True
rootFileNameT   = "NEATL36_TESTRUN_1d25h-m_3DT-thetao_*.nc_ZNB"                         # root for temperature filenames
fileListT       =  sorted(glob.glob(os.path.join(rootFolder,rootFileNameT)))            # list of salinity files
varNameT        = "thetao"
debug           = True
transectName    = "EAST2"



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

    try:
        e3v = dsZCoords["e3v_0"].isel(y=idxList[0],x=slice(idxList[1],idxList[2]))
    except:      
        e3v = dsZCoords["e3v_0"].isel(Y=idxList[0],X=slice(idxList[1],idxList[2]))

    # the variable might be umask or just "mask"
    try:
        vmask = dsVMask["mask"].isel(y=idxList[0],x=slice(idxList[1],idxList[2]))
    except:
        vmask = dsVMask["vmask"].isel(y=idxList[0],x=slice(idxList[1],idxList[2]))

    try:
        nz = dsZCoords.sizes["z"]
    except:      
        nz = dsZCoords.sizes["Z"]

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

    try:
        e3u = dsZCoords["e3u_0"].isel(x=idxList[0],y=slice(idxList[1],idxList[2]))
    except:      
        e3u = dsZCoords["e3u_0"].isel(X=idxList[0],Y=slice(idxList[1],idxList[2]))

    # the variable might be umask or just "mask"
    try:
        umask = dsUMask["mask"].isel(x=idxList[0],y=slice(idxList[1],idxList[2]))
    except:
        umask = dsUMask["umask"].isel(x=idxList[0],y=slice(idxList[1],idxList[2]))

    try:
        nz = dsZCoords.sizes["z"]
    except:
        nz = dsZCoords.sizes["Z"]

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
         ,transectName = "transectName"
         ,idxList = []                  # indexes identifying the transect: sn -> [xidx,yidx1,yidx2] or we -> [yidx,xidx1,xidx2]
         ,coordsPath = ""               # path of dataset with horizontal scale factors 
         ,zCoordsPath = ""              # path of dataset with vertical scale factors
         ,uMaskPath = ""                # path of file with mask
         ,vMaskPath = ""                
         ,Ntasks = 1                    # number of tasks to subdivide computations
         ,outFile = "out.nc"            # output file name
         ,simulationName = ""           # name of experiment 
         ,calcSaltTransp = False
         ,fileListS = []
         ,varNameS  = "so"
         ,calcHeatTransp = False
         ,fileListT = []
         ,varNameT = "thetao"
         ,debug = False
         ):          
    
    debugFile   = "transport_from_outputs_debug.dat"

    if debug:
        if os.path.isfile(debugFile):
            os.system("rm %s" % debugFile)

    volConversion   = 1e-6  # m^3/s to Sv
    rhoWater        = 1025  # kg / m^3
    cpWater         = 3996  # J/Kg/K
    saltConversion  = rhoWater * 1e-12
    heatConversion  = 1e-15

    Nfiles  = len(fileList)
    NfilesS = len(fileListS)
    NfilesT = len(fileListT)

    # do some checks
    if Nfiles == 0:
        raise Exception("Error: File list is empty. Check the path of files")        

    if calcSaltTransp == True:
        if NfilesS == 0:
            raise Exception("Error: calcSaltTransp is True and Salinity file list is empty. You set Check the path of files")
        if NfilesS != Nfiles:
            raise Exception("Error: the number of Salinity Files (%s) must be the same as Velocity files (%d)" % (NfilesS,Nfiles))

    if calcHeatTransp == True:
        if NfilesT == 0:
            raise Exception("Error: calcHeatTransp is True and Temperature file list is empty. You set Check the path of files")
        if NfilesT != Nfiles:
            raise Exception("Error: the number of Temperature Files (%s) must be the same as Velocity files (%d)" % (NfilesT,Nfiles))

    print("Processing %d files" % Nfiles )
    if debug:
        with open(debugFile,"a") as f:
            f.write("Processing %d files\n" % Nfiles )

    if transectType not in ["sn","we"]:
        raise Exception("Error: transectType not recognized")

    if len(idxList) != 3:
        raise Exception("Error: len(idxList) must be 3")        

    # calculate areas of faces
    if transectType == "sn":
        faces = calcYFaces(coordsPath=coordsPath,zCoordsPath=zCoordsPath,uMaskPath=uMaskPath,idxList=idxList)

    if transectType == "we":
        faces = calcXFaces(coordsPath=coordsPath,zCoordsPath=zCoordsPath,vMaskPath=vMaskPath,idxList=idxList)

    # List to collect 1D arrays
    daListTotalVolume       = []    # volume
    daListDirection1Volume  = []    
    daListDirection2Volume  = []
    daListTotalSalt         = []    # salt
    daListDirection1Salt    = []
    daListDirection2Salt    = []
    daListTotalHeat         = []    # heat
    daListDirection1Heat    = []
    daListDirection2Heat    = []

    # list to collect 3D arrays
    daListVolume3D          = [] # here I save (t,z,y) or (t,z,x) volume transport transect
    daListSalt3D            = [] # here I save (t,z,y) or (t,z,x) salt transport transect
    daListHeat3D            = [] # here I save (t,z,y) or (t,z,x) heat transport transect

    for task in range(Ntasks):

        i0,i1 = chunk_bounds(Ntasks=Ntasks,task=task,number=Nfiles)    

        fileListSubset = fileList[i0:i1]     

        print ("task: %s - processing %d files \n" % (task,i1-i0))

        if debug:
            with open(debugFile,"a") as f:
                f.write("Opening velocity files\n")

        dsVel = xr.open_mfdataset(fileListSubset)
        if debug:
            with open(debugFile,"a") as f:
                f.write("Velocity Dataset Open\n")

        if calcSaltTransp:
            if debug:
                with open(debugFile,"a") as f:
                    f.write("Opening salinity files\n")
            dsSalt = xr.open_mfdataset(fileListS)
            if debug:
                with open(debugFile,"a") as f:
                    f.write("Salinity Dataset Open\n")

        if calcHeatTransp:
            if debug:
                with open(debugFile,"a") as f:
                    f.write("Opening Temperture files\n")
            dsTemp = xr.open_mfdataset(fileListT)
            if debug:
                with open(debugFile,"a") as f:
                    f.write("Temperature Dataset Open\n")


        if transectType == "sn":
            # multiply by faces and sum 
            volumeTransport3D           = (dsVel[varName].isel(x=idxList[0],y=slice(idxList[1],idxList[2])) * faces ) # (t,z,y)
            volumeTransportTotal        = volumeTransport3D.sum(dim=["depthu","y"])
            volumeTransportDirection1   = volumeTransport3D.where(volumeTransport3D>=0,0).sum(dim=["depthu","y"])
            volumeTransportDirection2   = volumeTransport3D.where(volumeTransport3D<0,0).sum(dim=["depthu","y"])
            # Note: the command where(transport>0,0) to calculate positive transport
            # seems counter-intuitive, and actually is. Normally, to calculate positive 
            # transport one should filter the negative values, like:
            #
            #   transport_direction1 = transport3D.where(transport3D<0,0)
            #
            # this way I should get only positive values. Actually what I get is all the 
            # way around. That's why the code lines above are opposite

            if debug:
                with open(debugFile,"a") as f:
                    f.write("Volume Transport Computed\n")

            if calcSaltTransp:                                # (t,z,y)                             (t,z,y)
                saltTransport3D         = volumeTransport3D * dsSalt[varNameS].isel(x=idxList[0],y=slice(idxList[1],idxList[2])).data #(t,z,y)
                saltTransportTotal      = saltTransport3D.sum(dim=["depthu","y"])
                saltTransportDirection1 = saltTransport3D.where(saltTransport3D>=0,0).sum(dim=["depthu","y"])
                saltTransportDirection2 = saltTransport3D.where(saltTransport3D<0,0).sum(dim=["depthu","y"])
                if debug:
                    with open(debugFile,"a") as f:
                        f.write("Salt Transport Computed\n")

            if calcHeatTransp:                                      # (t,z,y)                             (t,z,y)
                heatTransport3D         = rhoWater * cpWater * volumeTransport3D * dsTemp[varNameT].isel(x=idxList[0],y=slice(idxList[1],idxList[2])).data #(t,z,y)
                heatTransportTotal      = heatTransport3D.sum(dim=["depthu","y"])
                heatTransportDirection1 = heatTransport3D.where(heatTransport3D>=0,0).sum(dim=["depthu","y"])
                heatTransportDirection2 = heatTransport3D.where(heatTransport3D<0,0).sum(dim=["depthu","y"])
                if debug:
                    with open(debugFile,"a") as f:
                        f.write("Heat Transport Computed\n")

        if transectType == "we":    
            # multiply by faces and sum 
            volumeTransport3D = (dsVel[varName].isel(y=idxList[0],x=slice(idxList[1],idxList[2])) * faces )          #(t,z,x)
            volumeTransportTotal      = volumeTransport3D.sum(dim=["depthv","x"])
            volumeTransportDirection1 = volumeTransport3D.where(volumeTransport3D>=0,0).sum(dim=["depthv","x"])       
            volumeTransportDirection2 = volumeTransport3D.where(volumeTransport3D<0,0).sum(dim=["depthv","x"])

            if debug:
                with open(debugFile,"a") as f:
                    f.write("Volume Transport Computed\n")

            if calcSaltTransp:                                 # (t,z,y)                             (t,z,y)
                saltTransport3D         = volumeTransport3D * dsSalt[varNameS].isel(y=idxList[0],x=slice(idxList[1],idxList[2])).data #(t,z,x)
                saltTransportTotal      = saltTransport3D.sum(dim=["depthv","x"])
                saltTransportDirection1 = saltTransport3D.where(saltTransport3D>=0,0).sum(dim=["depthv","x"])
                saltTransportDirection2 = saltTransport3D.where(saltTransport3D<0,0).sum(dim=["depthv","x"])

                if debug:
                    with open(debugFile,"a") as f:
                        f.write("Salt Transport Computed\n")

            if calcHeatTransp:                                     # (t,z,y)                             (t,z,y)
                heatTransport3D         = rhoWater * cpWater * volumeTransport3D * dsTemp[varNameT].isel(y=idxList[0],x=slice(idxList[1],idxList[2])).data #(t,z,x)
                heatTransportTotal      = heatTransport3D.sum(dim=["depthv","x"])
                heatTransportDirection1 = heatTransport3D.where(heatTransport3D>=0,0).sum(dim=["depthv","x"])
                heatTransportDirection2 = heatTransport3D.where(heatTransport3D<0,0).sum(dim=["depthv","x"])

                if debug:
                    with open(debugFile,"a") as f:
                        f.write("Heat Transport Computed\n")

        volumeTransportTotal.name = "volume_transport_total"
        daListTotalVolume.append(volumeTransportTotal) 
        volumeTransportDirection1.name = "volume_transport_direction1"
        daListDirection1Volume.append(volumeTransportDirection1)
        volumeTransportDirection2.name = "volume_transport_direction2"
        daListDirection2Volume.append(volumeTransportDirection2)
        volumeTransport3D.name = "volume_transport"
        daListVolume3D.append(volumeTransport3D)

        if calcSaltTransp:
            saltTransportTotal.name = "salt_transport_total"
            daListTotalSalt.append(saltTransportTotal)
            saltTransportDirection1.name = "salt_transport_direction1"
            daListDirection1Salt.append(saltTransportDirection1)
            saltTransportDirection2.name = "salt_transport_direction2"
            daListDirection2Salt.append(saltTransportDirection2)
            saltTransport3D.name = "salt_transport"
            daListSalt3D.append(saltTransport3D)

        if calcHeatTransp:
            heatTransportTotal.name = "heat_transport_total"
            daListTotalHeat.append(heatTransportTotal)
            heatTransportDirection1.name = "heat_transport_direction1"
            daListDirection1Heat.append(heatTransportDirection1)
            heatTransportDirection2.name = "heat_transport_direction2"
            daListDirection2Heat.append(heatTransportDirection2)
            heatTransport3D.name = "heat_transport"
            daListHeat3D.append(heatTransport3D)

        dsVel.close()
        if calcSaltTransp:
            dsSalt.close()

        if calcHeatTransp:
            dsTemp.close()

    # concatenate datasets from each task
    if debug:
        with open(debugFile,"a") as f:
            f.write("Concatenating DataArrays\n")
    print("Concatenating DataArrays")
    dsOut = xr.Dataset()
    #dsOut = xr.concat(ds_list, dim="time_counter")
    daOutTotalVolume        = xr.concat(daListTotalVolume,dim="time_counter")
    daOutTotalVolume        = daOutTotalVolume * volConversion     
    daOutTotalVolume        = daOutTotalVolume.assign_attrs({"units":"Sv","coordinates":"time_counter"})

    daOutDirection1Volume   = xr.concat(daListDirection1Volume,dim="time_counter")
    daOutDirection1Volume   = daOutDirection1Volume * volConversion
    daOutDirection1Volume   = daOutDirection1Volume.assign_attrs({"units":"Sv","coordinates":"time_counter"})

    daOutDirection2Volume   = xr.concat(daListDirection2Volume,dim="time_counter")
    daOutDirection2Volume   = daOutDirection2Volume * volConversion
    daOutDirection2Volume   = daOutDirection2Volume.assign_attrs({"units":"Sv","coordinates":"time_counter"})

    daOutVolume             = xr.concat(daListVolume3D,dim="time_counter")
    daOutVolume             = daOutVolume * volConversion
    if transectType == "sn":
        daOutVolume = daOutVolume.assign_attrs({"units":"Sv","coordinates":"time_counter depthu nav_lat"})     
    if transectType == "we":
        daOutVolume = daOutVolume.assign_attrs({"units":"Sv","coordinates":"time_counter depthv nav_lon"})     

    # populate dataset
    dsOut["volume_transport_total"]         = daOutTotalVolume
    dsOut["volume_transport_direction1"]    = daOutDirection1Volume
    dsOut["volume_transport_direction2"]    = daOutDirection2Volume
    dsOut["volume_transport"]               = daOutVolume

    if calcSaltTransp:
        daOutTotalSalt          = xr.concat(daListTotalSalt,dim="time_counter")
        daOutTotalSalt          = daOutTotalSalt * saltConversion 
        daOutTotalSalt          = daOutTotalSalt.assign_attrs({"units":"10^9 Kg/s","coordinates":"time_counter"})

        daOutDirection1Salt     = xr.concat(daListDirection1Salt,dim="time_counter")
        daOutDirection1Salt     = daOutDirection1Salt * saltConversion 
        daOutDirection1Salt     = daOutDirection1Salt.assign_attrs({"units":"10^9 Kg/s","coordinates":"time_counter"})

        daOutDirection2Salt     = xr.concat(daListDirection2Salt,dim="time_counter")
        daOutDirection2Salt     = daOutDirection2Salt * saltConversion 
        daOutDirection2Salt     = daOutDirection2Salt.assign_attrs({"units":"10^9 Kg/s","coordinates":"time_counter"})

        daOutSalt             = xr.concat(daListSalt3D,dim="time_counter")
        daOutSalt             = daOutSalt * saltConversion
        if transectType == "sn":
            daOutSalt = daOutSalt.assign_attrs({"units":"10^9 Kg/s","coordinates":"time_counter depthu nav_lat"})     
        if transectType == "we":
            daOutSalt = daOutSalt.assign_attrs({"units":"10^9 Kg/s","coordinates":"time_counter depthv nav_lon"})     

        # populate dataset
        dsOut["salt_transport_total"]       = daOutTotalSalt
        dsOut["salt_transport_direction1"]  = daOutDirection1Salt
        dsOut["salt_transport_direction2"]  = daOutDirection2Salt
        dsOut["salt_transport"]             = daOutSalt
 
    if calcHeatTransp:
        daOutTotalHeat          = xr.concat(daListTotalHeat,dim="time_counter")
        daOutTotalHeat          = daOutTotalHeat * heatConversion 
        daOutTotalHeat          = daOutTotalHeat.assign_attrs({"units":"10^15 Watt","coordinates":"time_counter"})

        daOutDirection1Heat     = xr.concat(daListDirection1Heat,dim="time_counter")
        daOutDirection1Heat     = daOutDirection1Heat * heatConversion 
        daOutDirection1Heat     = daOutDirection1Heat.assign_attrs({"units":"10^15 Watt","coordinates":"time_counter"})

        daOutDirection2Heat     = xr.concat(daListDirection2Heat,dim="time_counter")
        daOutDirection2Heat     = daOutDirection2Heat * heatConversion 
        daOutDirection2Heat     = daOutDirection2Heat.assign_attrs({"units":"10^15 Watt","coordinates":"time_counter"})

        daOutHeat             = xr.concat(daListHeat3D,dim="time_counter")
        daOutHeat             = daOutHeat * heatConversion
        if transectType == "sn":
            daOutHeat = daOutHeat.assign_attrs({"units":"10^15 Watt","coordinates":"time_counter depthu nav_lat"})     
        if transectType == "we":
            daOutHeat = daOutHeat.assign_attrs({"units":"10^15 Watt","coordinates":"time_counter depthv nav_lon"})     

        # populate dataset
        dsOut["heat_transport_total"]       = daOutTotalHeat
        dsOut["heat_transport_direction1"]  = daOutDirection1Heat
        dsOut["heat_transport_direction2"]  = daOutDirection2Heat
        dsOut["heat_transport"]             = daOutHeat
   
    if debug:
        with open(debugFile,"a") as f:
            f.write("Output Dataset Population Complete\n")

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

    # add transect type attribute
    dsOut   = dsOut.assign_attrs({"transect_type":transectType})
    dsOut   = dsOut.assign_attrs({"source":"NEMO Model Outputs"})
    dsOut   = dsOut.assign_attrs({"simulation_name":simulationName})
    dsOut   = dsOut.assign_attrs({"transect_name":transectName})

    # save dataset
    print("Saving Output")
    if debug:
        with open(debugFile,"a") as f:
            f.write("Saving Output (the .compute() phase is long..)\n")
    #dsOut.to_netcdf(outFile)
    dsOut.compute().to_netcdf(outFile)

if __name__ == "__main__":
    main(fileList = fileList
        ,varName = varName
        ,transectType = transectType
        ,transectName = transectName
        ,idxList = idxList
        ,coordsPath= coordsPath
        ,zCoordsPath= zCoordsPath
        ,uMaskPath = uMaskPath
        ,vMaskPath = uMaskPath
        ,Ntasks = Ntasks
        ,outFile = outFile
        ,simulationName = simulationName
        ,calcSaltTransp = calcSaltTransp
        ,fileListS = fileListS
        ,varNameS = varNameS
        ,calcHeatTransp= calcHeatTransp
        ,fileListT= fileListT
        ,varNameT= varNameT
        ,debug = debug
        )

