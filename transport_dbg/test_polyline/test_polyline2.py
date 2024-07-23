import xarray as xr 
import numpy as np 
import matplotlib.pyplot as plt 
import cartopy.crs as ccrs
from math import radians, cos, sin, asin, sqrt
#from index_extraction import createTransectLine,distance,haversine
from transport_tools import calcNEMOTransport
from os.path import join

transectName = "arkona"        # name of transect (for outFile)
P1 = (13.33,54.45)
P2 = (13.33,55.45)

transectName = "arkona2_test"        # name of transect (for outFile)
P1 = (13.6,54.54)
P2 = (13.6,55.45)

root = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/CDO/outputs_bal/"
rootStatic = "/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/static_files/"

resolution  = 0.4* (1./36)      # it's better to have resolution higher than the native mesh 
pathCoords2D = rootStatic + "coordinates_east2.nc"              # path of coordinates.nc
pathZMesh = rootStatic+"mask_testrun/create_meshmask/mesh_mask_east2.nc"                 # path of file with vertical scale factors   #e3t
pathU = root + "NEATL36_east2_BAL_CDO_????????.nc"                     # path with U velocity files (you can use wildcards)
pathV = root + "NEATL36_east2_BAL_CDO_????????.nc"                     # path with V velocity files (you can use wildcards)  
pathS = root + "NEATL36_east2_BAL_CDO_????????.nc"                     # path of Salinity files (you can use wildcards)
pathMaskU = rootStatic + "mask_testrun/mask_gridU_east2.nc"                
pathMaskV = rootStatic + "mask_testrun/mask_gridV_east2.nc"                        
outFileRoot = "out_transport_%s.nc" # outfile root

computeS = True               # compute salinity transport (if True you must set pathS)       
varNameS = "vosaline"           # name of salinity variable in files 
verbose = True
fileType = "nemo_obc"           
makePlot = True
obcType = "E"


calcNEMOTransport(                                  # Wrapper to the main function
                        P1 = P1,                            # start point (lon,lat) tuple
                        P2 = P2,                            # end point (lon,lat) tuple
                        resolution = resolution,            # resolution (in Deg) of user-defined transect line (tip: multiply the mesh resolution by a factor 0.4 -> resolution = 0.4*mesh_resolution)
                        pathCoords2D = pathCoords2D,        # path of coordinates.nc
                        pathZMesh = pathZMesh,              # path of file with vertical scale factors   #e3t,e3u,e3v..
                        pathU = pathU,                      # path with U velocity files (you can use wildcards)
                        pathV = pathV,                      # path with V velocity files (you can use wildcards)  
                        pathS = pathS,                      # path of Salinity files (you can use wildcards)
                        pathMaskU = pathMaskU,              # path of Umask file   
                        pathMaskV = pathMaskV,              # path of Vmask file             
                        outFileRoot = outFileRoot,                  # outfile root
                        transectName = transectName,        # name of transect (for outFile)
                        computeS = computeS,                # compute salinity transport (if True you must set pathS)   
                        varNameS = varNameS,                # name of salinity variable in fileS    
                        verbose = verbose,                  # print more info
                        fileType = fileType,                # file type is one of "nemo_out" or "nemo_obc"
                        makePlot= makePlot,                 # make a plot of mesh geometry, transect line and segments that approximate the latter
                        obcType=obcType                     # in case of "nemo_obc" fileType, specify W/E (in the latter case, coordinates are reversed along x dim)
                    )

