import numpy as np
import xarray as xr 

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r

def distance(lon1,lat1,lon2,lat2):
    '''
    This funtion provides the distance calculation on the sphere surface (source: https://www.kompf.de/gps/distcalc.html)
    args:
        lat1,lon1 = pair of values on reference line
        lat2,lat2 = pair of values for which the distance to the reference line should be tested
    returns: 
        distance of (lat1,lon1) and (lat2,lon2) in kilometers
    '''
    lat11 = np.radians(lat1)
    lon11 = np.radians(lon1)
    lat22 = np.radians(lat2)
    lon22 = np.radians(lon2)
    a = 6378.388
    part1 = np.sin(lat11)*np.sin(lat22) 
    part2 = np.cos(lat11)*np.cos(lat22)*np.cos(lon22-lon11)
    term = part1+part2
    # following if/else statement is necessary due to numerical issues
    if type(term) != float:   
        term = np.where(term>1,1,term)
    else:
        if term >1:
            term=1
    distance = a * np.arccos(term)
    return distance

def create_transect_line(P1=(),P2=(),resolution=0.1):

    if P1[0] == P2[0]:  # south north
        yTransect = np.arange(P1[1],P2[1],resolution)
        xTransect = P1[0]*np.ones(len(yTransect))
    elif P1[1] == P2[1]:    #west-east
        xTransect = np.arange(P1[0],P2[0],resolution)
        yTransect = P1[1]*np.ones(len(xTransect))
    else:                   # generic
        xTransect = np.arange(P1[0],P2[0],resolution)
        yTransect = np.arange(P1[1],P2[1],resolution)

    #line = xr.DataArray(data=1,dims=['y','x'],coords=[yTransect,xTransect])

    lons = xr.DataArray(xTransect,dims="transect")
    lats = xr.DataArray(yTransect,dims="transect")

    return lons,lats     

def calc_next_idx(prevIdx=(0,0),direction="n"):
    nextIdx = prevIdx
    # change temporarily to list (tuple cannot be changed)
    nextIdx = list(nextIdx)
    if direction == "n":
        nextIdx[0] += 1
    if direction == "s":
        nextIdx[0] -= 1
    if direction == "e":
        nextIdx[1] += 1
    if direction == "w":
        nextIdx[1] -= 1
    # convert back again to tuple
    nextIdx = tuple(nextIdx)

    return nextIdx

def next_direction(prevDirection="n",
                   prevIdx=(0,0),
                   latsF=xr.DataArray(),            # lats of F points in mesh
                   lonsF=xr.DataArray(),            # lons of F points in mesh
                   latsTransect=xr.DataArray(),     # lats of transect points
                   lonsTransect=xr.DataArray(),     # lons of transect points
                   verbose=False
                   ):

    nx = lonsF.sizes["x"]
    ny = lonsF.sizes["y"]

    directions      = ["s","e","n","w"]
    anti_directions = ["n","w","s","e"]
    # the next direction cannot be the opposite to the previous
    directions.remove(anti_directions[directions.index(prevDirection)])

    if verbose:
        print("directions")
        print(directions)

    # initialize with huge distances
    distances = 1000*np.ones(len(directions)) # 1000 km
    
    lonA = lonsF.isel(y=prevIdx[0],x=prevIdx[1])
    latA = latsF.isel(y=prevIdx[0],x=prevIdx[1])

    # store temporary idx
    idxListTemp = [(0,0),(0,0),(0,0)]

    # loop through the 3 possible directions left
    for i,dir in enumerate(directions):
        nextIdxTemp = calc_next_idx(prevIdx=prevIdx,direction=dir)

        if verbose:
            print("     i,direction",i,dir)        
            print("     nextIdxTemp")
            print("     ",nextIdxTemp)

        idxListTemp[i] = nextIdxTemp

        # check on bounds
        if nextIdxTemp[0] <0 or nextIdxTemp[0] >= ny: 
            continue
        if nextIdxTemp[1] <0 or nextIdxTemp[1] >= nx: 
            continue    

        lonB = lonsF.isel(y=nextIdxTemp[0],x=nextIdxTemp[1])
        latB = latsF.isel(y=nextIdxTemp[0],x=nextIdxTemp[1])
        # calc distances of nextPoint from transect line
        distances[i] = np.min(haversine(lonB,latB,lonsTransect,latsTransect)).item()
    
    # the next grid points is the the one closer to the transect line
    nextIdx         = idxListTemp[np.argmin(distances)]
    newDirection    = directions[np.argmin(distances)]

    if verbose:
        print("   --")
        print("   distances:" ,distances)
        print("   nextIdx:", nextIdx)
        print("   newDirection:", newDirection)

    return nextIdx,newDirection

def extract_f_indices(
                    latsTransect=xr.DataArray(),
                    lonsTransect=xr.DataArray(),
                    latsF=xr.DataArray(),
                    lonsF=xr.DataArray(),
                    verbose=False,
                    nCountMax=1000
                    ):
    
    if verbose:
        print("Extracting indices of F-points")

    # find initial F point in mesh (closest to 1st transect point)
    dist_start  = haversine(lonsTransect.isel(transect=0),latsTransect.isel(transect=0),lonsF,latsF)
    startIdx    = np.unravel_index(np.argmin(dist_start.values),dist_start.shape)

    if verbose:
        print("startIdx")
        print(startIdx)

    # find final F point in mesh (closest to last transect point)
    dist_end    = haversine(lonsTransect.isel(transect=-1),latsTransect.isel(transect=-1),lonsF,latsF)
    endIdx      = np.unravel_index(np.argmin(dist_end.values),dist_end.shape)

    if verbose:
        print("endIdx")
        print(endIdx)

    # the loop should last until the idx_end is found
    prevIdx = startIdx
    prevDirection = "w"

    indices = []
    indices.append(startIdx)

    nCount = 0
    while True:
        nCount += 1
        nextIdx,newDirection = next_direction(prevDirection=prevDirection,
                                prevIdx=prevIdx,
                                latsF=latsF,                    # lats of F points in mesh
                                lonsF=lonsF,                    # lons of F points in mesh
                                latsTransect=latsTransect,      # lats of transect points
                                lonsTransect=lonsTransect       # lons of transect points
                                )

        prevDirection = newDirection
        prevIdx = nextIdx
        #print("next Idx")
        #print(nextIdx)
        if nCount >= nCountMax:
            raise Exception("Problem with search of F indices. Loop too long! you have to make some checks")
        # update indices list 
        indices.append(nextIdx)         
        if nextIdx == endIdx:
            print("I have found the last point")
            break

    return indices

def check_f_indices(indices=[],verbose=False):
    # check integrity of extracted F indices
    # differences of consecutive must be 1
    # in absolute value in either y or x direction

    # indices is a list of tuples
    # [(yidx1,xidx1),(yidx2,xidx2),(yidx1,xidx1),...,(yidxN,xidxN)]
    if verbose:
        print("checking itegrity of F indices")

    if len(indices) == 0:
        raise Exception("Error: check_f_indices: input list has 0 length")

    indicesOK = True

    for i in range(1,len(indices)):
        xdiff = indices[i][1] - indices[i-1][1]
        ydiff = indices[i][0] - indices[i-1][0]
        diff_tot = abs(xdiff) + abs(ydiff)
        if diff_tot != 1:
            indicesOK = False
            break

    if indicesOK == False: 
        print("Error: problem with F indices extracted")             
        print("       there is at least 1 case with a jump")

    return indicesOK         

def extract_uv_indices(indicesF=[],verbose=False):
    # extract indices of U,V segments from indices of F points

    if len(indicesF) == 0:
        raise Exception("Error: extract_uv_indices: input list has 0 length")

    nSegments = len(indicesF)-1
    if verbose:
        print("extract_uv_indices: nSegments: %d" % nSegments)

    # initialize nSegments along U,V
    nSegmentsU  = 0
    nSegmentsV  = 0

    indicesU = []
    indicesV = []

    for i in range(len(indicesF)-1):
        xdiff = indicesF[i+1][1] - indicesF[i][1]
        ydiff = indicesF[i+1][0] - indicesF[i][0]
        if abs(xdiff) == 1: # V-segment
            nSegmentsV += 1
            if xdiff > 0:
                new_v_idx = (indicesF[i][0],indicesF[i][1]+xdiff)
            if xdiff < 0:
                new_v_idx = (indicesF[i][0],indicesF[i][1])
            indicesV.append(new_v_idx)    
        if abs(ydiff) == 1: # U-segment      
            nSegmentsU += 1     
            if ydiff > 0:
                new_u_idx = (indicesF[i][0]+ydiff,indicesF[i][1])
            if ydiff < 0:
                new_u_idx = (indicesF[i][0],indicesF[i][1])
            indicesU.append(new_u_idx)

    # check on number of segmets
    if nSegments != nSegmentsU + nSegmentsV:
        raise Exception("Error: extract_uv_indices: nU+nV != nSegments")

    if verbose:
        print("extract_uv_indices:")
        print("nSegmentsU: " ,nSegmentsU )
        print("nSegmentsV: " ,nSegmentsV )

    return indicesV,indicesU

    # t_indices are time-varying, they depend on sign of velocity field
def extract_t_indices(
                    
                    ):
    pass

"""
def fix_f_indices(indexes):
    new_x_idx = []
    new_y_idx = []
    # add the first pair
    new_x_idx.append(indexes[0,1])
    new_y_idx.append(indexes[0,0])
    for i,(yidx,xidx) in enumerate(zip(indexes[:,0],indexes[:,1])):
        print("i,xidx,yidx %d,%d,%d" %(i,yidx,xidx))
        if i < indexes.shape[0]-1:
            yidx_next = indexes[i+1,0]
            xidx_next = indexes[i+1,1]
            ydiff = yidx_next - yidx
            xdiff = xidx_next - xidx
            tot_idx_diff = abs(xdiff) + abs(ydiff)
            if tot_idx_diff == 0:
                raise Exception("Error: problem! there is still a duplicate!")
            if tot_idx_diff == 1:
                print("Okay")
                new_x_idx.append(xidx_next)
                new_y_idx.append(yidx_next)
            if tot_idx_diff == 2:
                print("there is a single jump. Let's check if is along x or y or diagonal")
                print("xdiff: %d"  % xdiff)
                print("ydiff: %d"  % ydiff)
                if abs(xdiff) == abs(ydiff):
                    #  y+1   o         o      or     o        o 
                    #           .                          .
                    #               .                   .
                    #   y    o         o             o        o
                    #        x        x+1            x       x+1
                    print("it's diagonal single jump. I add a pair of indexes along y")
                    new_x_idx.append(xidx)
                    new_y_idx.append(yidx_next)
                if abs(xdiff) == 2:
                    #  y+1   o         o         o      or     o         o        o
                    #                                     
                    #                                     
                    #   y    o         o         o             o         o        o
                    #        ....................>             <................... 
                    #        x        x+1       x+2            x       x+1       x+2 
                    print("jump along x")
                    new_x_idx.append(xidx + int(0.5*xdiff))
                    new_y_idx.append(yidx_next)
                if abs(ydiff) == 2:
                    print("jump along y")
                    new_x_idx.append(xidx_next)
                    new_y_idx.append(yidx + int(0.5*ydiff))

            if tot_idx_diff > 2:
                raise Exception("Error: I cannot this case. Increase the resolution of transect")


    newIndices = np.vstack((np.asarray(new_y_idx),np.asarray(new_x_idx))).T

    return newIndices
"""
