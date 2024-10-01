import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
from math import radians, cos, sin, asin, sqrt
#from index_extraction import createTransectLine,distance,haversine
from transport_tools import calcNEMOTransport
from os.path import join



transectName = "arkona2_test"        # name of transect (for outFile)
P1 = (13.6,54.54)
P2 = (13.6,55.45)





transectName = "SamsoBelt"
P1 = (10.70,56.28)
P2 = (11.42,55.77)

transectName = "East2"
P1 = (13.46,54.43)
P2 = (13.93,55.49)

transectName = "East2_rev"
P1 = (13.93,55.49)
P2 = (13.46,54.43)

transectName = "SamsoBelt_rev"
P1 = (11.42,55.77)
P2 = (10.70,56.28)


transectName = "GreatBelt_rev"             # name of transect (for outFile)
P1 = (11.25,55.50)
P2 = (10.70,55.35)

transectName = "kattegat_rev"
P1 = (12.65,56.85)
P2 = (10.15,56.85) 

transectName = "Skagerrak_rev"
P1 = (10.10,59.15)  
P2 = (10.10,57.50)  

transectName = "SOUND"
P1 = (12.50,55.85)
P2 = (12.95,55.85)

transectName = "SOUND_rev"
P1 = (12.95,55.85)
P2 = (12.50,55.85)

transectName = "arkona"             # name of transect (for outFile)
P1 = (13.33,54.45)
P2 = (13.33,55.45)

transectName = "arkona_rev"             # name of transect (for outFile)
P1 = (13.33,55.45)
P2 = (13.33,54.45)

transectName = "DARSS"
P1 = (12.00,54.80  )
P2 = (12.30,54.20 ) 

transectName = "DARSS_rev"
P1 = (12.30,54.20) 
P2 = (12.00,54.80)

transectName = "kattegat"
P1 = (10.15,56.85) 
P2 = (12.65,56.85)

transectName = "GreatBelt"             # name of transect (for outFile)
P1 = (10.70,55.35)
P2 = (11.25,55.50)

transectName = "Skagerrak"
P1 = (10.10,57.50)  
P2 = (10.10,59.15)  

root            = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh/"
rootStatic      = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh/static/"

resolution      = 0.4* (1./36)      # it's better to have resolution higher than the native mesh 
pathCoords2D    = rootStatic + "coordinates_NEATL36_ZNB.nc"                                 # path of coordinates.nc
pathZMesh       = rootStatic + "mesh_mask_testrun2_ZNB.nc"                                  # path of file with vertical scale factors   #e3t
pathMaskU       = rootStatic + "mask_gridU_ZNB.nc"                
pathMaskV       = rootStatic + "mask_gridV_ZNB.nc"                        
outFileRoot     = "test_out_transport_%s.nc" # outfile root

pathU           = root + "NEATL36_TESTRUN_1d25h-m_3DU-uo_2021????-????????.nc_ZNB"        # path with U velocity files (you can use wildcards)
pathV           = root + "NEATL36_TESTRUN_1d25h-m_3DV-vo_2021????-????????.nc_ZNB"        # path with V velocity files (you can use wildcards)  
pathS           = root + "NEATL36_TESTRUN_1d25h-m_3DT-so_2021????-????????.nc_ZNB"        # path of Salinity files (you can use wildcards)

pathU           = root + "NEATL36_TESTRUN_1d25h-m_3DU-uo_20220[1-6]??-????????.nc_ZNB"        # path with U velocity files (you can use wildcards)
pathV           = root + "NEATL36_TESTRUN_1d25h-m_3DV-vo_20220[1-6]??-????????.nc_ZNB"        # path with V velocity files (you can use wildcards)  
pathS           = root + "NEATL36_TESTRUN_1d25h-m_3DT-so_20220[1-6]??-????????.nc_ZNB"        # path of Salinity files (you can use wildcards)

pathU           = root + "NEATL36_TESTRUN_1d25h-m_3DU-uo_????????-????????.nc_ZNB"        # path with U velocity files (you can use wildcards)
pathV           = root + "NEATL36_TESTRUN_1d25h-m_3DV-vo_????????-????????.nc_ZNB"        # path with V velocity files (you can use wildcards)  
pathS           = root + "NEATL36_TESTRUN_1d25h-m_3DT-so_????????-????????.nc_ZNB"        # path of Salinity files (you can use wildcards)

computeS        = False               # compute salinity transport (if True you must set pathS)       
varNameS        = "vosaline"           # name of salinity variable in files 
verboseLevel    = 1
fileType        = "nemo_out"           
makePlot        = False
obcType         = "E"   # to be removed

calcNEMOTransport(                                  # Wrapper to the main function
                        P1 = P1,                            # start point (lon,lat) tuple
                        P2 = P2,                            # end point (lon,lat) tuple
                        resolution = resolution,            # resolution (in Deg) of user-defined transect line (tip: multiply the mesh resolution by a factor 0.4 -> resolution = 0.4*mesh_resolution)
                        pathCoords2D = pathCoords2D,        # path of coordinates.nc
                        pathZMesh = pathZMesh,              # path of file with vertical scale factors   #e3t,e3u,e3v..
                        pathMaskU = pathMaskU,              # path of Umask file   
                        pathMaskV = pathMaskV,              # path of Vmask file             
                        pathU = pathU,                      # path with U velocity files (you can use wildcards)
                        pathV = pathV,                      # path with V velocity files (you can use wildcards)  
                        pathS = pathS,                      # path of Salinity files (you can use wildcards)
                        outFileRoot = outFileRoot,          # outfile root
                        transectName = transectName,        # name of transect (for outFile)
                        computeS = computeS,                # compute salinity transport (if True you must set pathS)   
                        varNameS = varNameS,                # name of salinity variable in fileS    
                        verboseLevel = verboseLevel,        # print more info
                        fileType = fileType,                # file type is one of "nemo_out" or "nemo_obc"
                        makePlot= makePlot,                 # make a plot of mesh geometry, transect line and segments that approximate the latter
                        obcType = obcType                   # in case of "nemo_obc" fileType, specify W/E (in the latter case, coordinates are reversed along x dim)
                    )

