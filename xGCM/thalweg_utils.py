import numpy as np 
import xarray as xr
from os.path import isfile

def open_thalweg_coords(pathIn=""):
    # read file with thalweg coords (lon lat)
    # file format is with space between lon and lat

    # lon1 lat1
    # lon2 lat2
    # lon3 lat3
    #...
    # lonN,latN

    if not isfile(pathIn):
        raise Exception("Error! file with thalweg coords not existing")

    thalweg_coords = np.genfromtxt(pathIn)

    lonsThalweg = xr.DataArray(thalweg_coords[:,0],dims="thalwegPoint")
    latsThalweg = xr.DataArray(thalweg_coords[:,1],dims="thalwegPoint")

    return (lonsThalweg,latsThalweg)

def find_thalweg_idx(lon2d,lat2d,lons_thalweg,lats_thalweg,mask=None,verbose=False):

    if len(lons_thalweg) != len(lats_thalweg):
        raise Exception("Error! lons_thalweg must have same size of lats_thalweg")

    ny,nx = lon2d.shape 

    points_list = []

    for lon,lat in zip(lons_thalweg,lats_thalweg):
        dist = np.sqrt((lon-lon2d)**2 + (lat2d-lat)**2)
        # set land points to huge distance
        if mask != None:
            dist[mask] = 1e12
        try:     
            (yidx,xidx) = np.unravel_index(np.argmin(dist.values),(ny,nx))
        except:     # case of numpy array    
            (yidx,xidx) = np.unravel_index(np.argmin(dist)       ,(ny,nx))
        if verbose:     
            print(yidx,xidx)

        points_list.append((yidx,xidx))

    if len(points_list) == 0:
        raise Exception("Error! found 0 thalweg points. Check coordinates")
        
    # check for equal indices
    for k in range(len(points_list)-1):
        if points_list[k+1] == points_list[k]:
            print("found equal indices %d,%d,%s,%s" % (k,k+1,points_list[k],points_list[k+1]))

    return points_list

def extract_thalweg(da=xr.DataArray(),pointsIndexes=[],verbose=False):
    # this function creates a DataArray of a thalweg.

    # Arguments:

    # 1) da with dims (t,z,y,x), (z,y,x) or (y,x)
    #
    # 2) list of tuples representing indexes on the mesh
    #   [(yidx0,xidx0),
    #    (yidx1,xidx1),
    #    (yidx2,xidx2),
    #       ....,
    #    (yidxN,xidxN)]     
    #
    # creates a new DataArray with 
    # shape: 
    #   (z,thalwegPoint)
    # or 
    #   (thalwegPoint)
    #

    # nP (points along thalweg)
    nPoints = len(pointsIndexes)
    if verbose:
        print("Extracting %d thalweg points " % nPoints)

    # convert list to 1d arrays
    pointsListY = np.array([point[0] for point in pointsIndexes])
    pointsListX = np.array([point[1] for point in pointsIndexes])

    # check dimension of input DataArray (It has to be (z,y,x), no
    # time dimension allowed)

    nDims = len(da.sizes)
    if verbose:
        print("# of dimensions of input array: %d \n sizes: %s" % (nDims,da.sizes) )

    if nDims not in [2,3,4]:
        raise Exception("Error! the number of dimension of da must be 4 (t,z,y,x), \
                        3 (z,y,x) or 2 (y,x). found: %s" % da.sizes)

    # create new shape
    oldShape = da.shape
    ny,nx = oldShape[-2:]
    newShape = oldShape[:-2] + (oldShape[-2]*oldShape[-1],)
    dataReshaped = da.values.reshape(newShape)

    """
    # reshape to flatten last 2 dims
    if nDims == 3:
        nz,ny,nx = da.shape
        dataReshaped = da.values.reshape( nz, ny * nx )
    if nDims == 2:        
        ny,nx = da.shape
        dataReshaped = da.values.reshape( ny * nx )
    """ 

    yxIndices = np.ravel_multi_index((pointsListY,pointsListX),(ny,nx))

    if nDims == 4:
        selectedValues = dataReshaped[:,:, yxIndices]        

        # find name of 1d depth coordinate (gdept_1d, gdepu_1d...)
        dep1dName = None
        for key in da.coords:
            if "gdep" in key and "_1d" in key:
                dep1dName = key
        if dep1dName == None:
            raise Exception("Error! I could not find the name of 1d depth coord")

        # create new DataArray (take z coordinate from input dA)
        newDa = xr.DataArray(data=selectedValues,
                    dims=(da.dims[0],da.dims[1],"thalwegPoint"),  #(z,thalwegPoint)
                    coords={da.dims[0]:da[da.dims[0]],
                            da.dims[1]:da[da.dims[1]],
                            "thalwegPoint": np.arange(nPoints),
                            dep1dName:da[dep1dName]}
                            )

    if nDims == 3:
        selectedValues = dataReshaped[:, yxIndices]
        # find name of 1d depth coordinate (gdept_1d, gdepu_1d...)
        dep1dName = None
        for key in da.coords:
            if "gdep" in key and "_1d" in key:
                dep1dName = key
        if dep1dName == None:
            raise Exception("Error! I could not find the name of 1d depth coord")
        
        # create new DataArray (take z coordinate from input dA)
        newDa = xr.DataArray(data=selectedValues,
                    dims=(da.dims[0],"thalwegPoint"),  #(z,thalwegPoint)
                    coords={da.dims[0]:da[da.dims[0]],
                            "thalwegPoint": np.arange(nPoints),
                            dep1dName:da[dep1dName]}
                            )

    if nDims == 2:        
        selectedValues = dataReshaped[yxIndices]

        # create new DataArray (take z coordinate from input dA)
        newDa = xr.DataArray(data=selectedValues,
                    dims=("thalwegPoint"),  
                    coords={
                            "thalwegPoint": np.arange(nPoints),
                            }
                            )

    return newDa

def calc_tangents_to_thalweg(lons_thalweg,lats_thalweg,sign=1.):

    #  P1
    #    o  \ 
    #     \  \ v1
    #      \ _\|      v2
    #       \   ----------->
    #        o ------------ o            
    #      P2              P3 \  \
    #                          \  \ v3
    #                           \ _\|
    #                            o 
    #                            P4

    # I use the coordinate of thalweg to calculate 
    # components of tangential vectors. 

    # I use downwind method, so from NPoints I get NPoints-1 vectors
    # The tangent in point N is extrapolated from point N-1.

    if len(lons_thalweg) != len(lats_thalweg):
        raise Exception("Error! lons_thalweg and lats_thalweg have different length")

    NPoints = len(lons_thalweg) 

    # initialize arrays to store tangent vectors
    vx = np.zeros(NPoints)
    vy = np.zeros(NPoints)

    # 
    if type(lons_thalweg) == xr.DataArray:
        lons1D = lons_thalweg.data 
        lats1D = lats_thalweg.data         
    else:
        lons1D = lons_thalweg
        lats1D = lats_thalweg         

    for iP in range(NPoints-1):
        dx = lons_thalweg[iP+1] - lons_thalweg[iP]
        dy = lats_thalweg[iP+1] - lats_thalweg[iP] 
        # normalize
        mod = np.sqrt(dx**2 + dy**2)
        dx,dy = list(map(lambda x: (1./mod)* x, [dx,dy]))
        # apply sign (+/- 1)
        vx[iP] = sign*dx
        vy[iP] = sign*dy

    # apply boundary condition (it's a bit crude ...)
    vx[-1] = vx[-2]
    vy[-1] = vy[-2]

    vx = xr.DataArray(data=vx,dims=["thalwegPoint"])
    vy = xr.DataArray(data=vy,dims=["thalwegPoint"])

    return vx,vy

