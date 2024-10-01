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

def geo2stereo(lam,phi):
    # convert from geographical to steregraphic projection
    xStereo = 2*np.cos(np.deg2rad(lam))*np.tan(0.25*np.pi - 0.5*np.deg2rad(phi))
    yStereo = 2*np.sin(np.deg2rad(lam))*np.tan(0.25*np.pi - 0.5*np.deg2rad(phi))

    return xStereo,yStereo

def create_transect_line(P1=(),P2=(),resolution=0.1):

    if P1[0] == P2[0]:  # south north
        yTransect = np.arange(P1[1],P2[1],resolution)
        if P2[1] < P1[1]:
            yTransect = np.arange(P1[1],P2[1],-resolution)
        xTransect = P1[0]*np.ones(len(yTransect))
    elif P1[1] == P2[1]:    #west-east
        xTransect = np.arange(P1[0],P2[0],resolution)
        if P2[0] < P1[0]:
            xTransect = np.arange(P1[0],P2[0],-resolution)
        yTransect = P1[1]*np.ones(len(xTransect))
    else:                   # generic
        # lons,lats must have same size, cannot use np.arange
        deltaLon = abs(P2[0]-P1[0])
        deltaLat = abs(P2[1]-P1[1])
        nPoints = int(max(np.round(deltaLon/resolution),np.round(deltaLat/resolution)))
        xTransect = np.linspace(P1[0],P2[0],nPoints)
        yTransect = np.linspace(P1[1],P2[1],nPoints)

    #line = xr.DataArray(data=1,dims=['y','x'],coords=[yTransect,xTransect])

    lons = xr.DataArray(xTransect,dims="transect")
    lats = xr.DataArray(yTransect,dims="transect")

    if lons.sizes != lats.sizes:
        raise Exception("Error: create_transect_line. dimensions lons.sizes (%d) != lats.sizes (%d)" % (lons.sizes["transect"],lats.sizes["transect"]))

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
                   verboseLevel = 0
                   ):

    nx = lonsF.sizes["x"]
    ny = lonsF.sizes["y"]

    directions      = ["s","e","n","w"]
    anti_directions = ["n","w","s","e"]
    # the next direction cannot be the opposite to the previous
    directions.remove(anti_directions[directions.index(prevDirection)])

    if verboseLevel>5:
        print("directions")
        print(directions)

    # initialize with huge distances
    distancesFromLine   = 1000*np.ones(len(directions)) # 1000 km
    distancesFromEnd    = 1000*np.ones(len(directions)) # 1000 km

    lonA = lonsF.isel(y=prevIdx[0],x=prevIdx[1])
    latA = latsF.isel(y=prevIdx[0],x=prevIdx[1])

    # store temporary idx
    idxListTemp = [(0,0),(0,0),(0,0)]

    # loop through the 3 possible directions left
    for i,dir in enumerate(directions):
        nextIdxTemp = calc_next_idx(prevIdx=prevIdx,direction=dir)

        if verboseLevel>5:
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
        distancesFromLine[i]    = np.min(haversine(lonB,latB,lonsTransect,latsTransect)).item()
        # distances of nextPoint from end point of transect line
        distancesFromEnd[i]     = np.min(haversine(lonB,latB,lonsTransect.isel(transect=-1),latsTransect.isel(transect=-1))).item()
    
    # the next grid points is calculated considering how close it is to the 
    # transect line and the final point of the transect. It seems to work well
    distances = distancesFromLine + distancesFromEnd
    nextIdx         = idxListTemp[np.argmin(distances)]
    newDirection    = directions[np.argmin(distances)]

    if verboseLevel>5:
        print("   --")
        print("   distances:" ,distances)
        print("   nextIdx:", nextIdx)
        print("   newDirection:", newDirection)

    return nextIdx,newDirection

def calcStartDirection(
        startIdx = (0,0),
        latsF=xr.DataArray(),            # lats of F points in mesh
        lonsF=xr.DataArray(),            # lons of F points in mesh
        latsTransect=xr.DataArray(),     # lats of transect points
        lonsTransect=xr.DataArray(),     # lons of transect points

                        ):

    directions = ["s","e","n","w"]



    pass    

def extract_f_indices(
                    latsTransect=xr.DataArray(),
                    lonsTransect=xr.DataArray(),
                    latsF=xr.DataArray(),
                    lonsF=xr.DataArray(),
                    nCountMax=200,
                    verboseLevel = 0
                    ):
    
    if verboseLevel>0:
        print("Extracting indices of F-points")

    # find initial F point in mesh (closest to 1st transect point)
    dist_start  = haversine(lonsTransect.isel(transect=0),latsTransect.isel(transect=0),lonsF,latsF)
    startIdx    = np.unravel_index(np.argmin(dist_start.values),dist_start.shape)

    if verboseLevel>1:
        print("startIdx")
        print(startIdx)
        print("Lat,Lon: %f,%f" % (latsF.isel(y=startIdx[0],x=startIdx[1]),lonsF.isel(y=startIdx[0],x=startIdx[1])))

    # find final F point in mesh (closest to last transect point)
    dist_end    = haversine(lonsTransect.isel(transect=-1),latsTransect.isel(transect=-1),lonsF,latsF)
    endIdx      = np.unravel_index(np.argmin(dist_end.values),dist_end.shape)

    if verboseLevel>1:
        print("endIdx")
        print(endIdx)
        print("Lat,Lon: %f,%f" % (latsF.isel(y=endIdx[0],x=endIdx[1]),lonsF.isel(y=endIdx[0],x=endIdx[1])))

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
                                lonsTransect=lonsTransect,      # lons of transect points
                                verboseLevel = verboseLevel
                                )

        prevDirection = newDirection
        prevIdx = nextIdx
        if verboseLevel > 5:
            print("extract_f_indixces:")
            print("nextDirection %s" % newDirection)
            print("next Idx")
            print(nextIdx)
            print("Lat,Lon: %f,%f" % (latsF.isel(y=nextIdx[0],x=nextIdx[1]),lonsF.isel(y=nextIdx[0],x=nextIdx[1])))
        if nCount >= nCountMax:
            raise Exception("Problem with search of F indices. Loop too long! you have to make some checks")
        # update indices list 
        indices.append(nextIdx)         
        if nextIdx == endIdx:
            print("I have found the last point")
            break

    return indices

def check_f_indices(indices=[],verboseLevel=False):
    # check integrity of extracted F indices
    # differences of consecutive must be 1
    # in absolute value in either y or x direction

    # indices is a list of tuples
    # [(yidx1,xidx1),(yidx2,xidx2),(yidx1,xidx1),...,(yidxN,xidxN)]
    if verboseLevel>0:
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

def extract_uv_indices(indicesF=[],verboseLevel=False):
    # extract indices of U,V segments from indices of F points

    if len(indicesF) == 0:
        raise Exception("Error: extract_uv_indices: input list has 0 length")

    nSegments = len(indicesF)-1
    if verboseLevel>0:
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

    if verboseLevel>0:
        print("extract_uv_indices:")
        print("nSegmentsU: " ,nSegmentsU )
        print("nSegmentsV: " ,nSegmentsV )

    return indicesV,indicesU

    # t_indices are time-varying, they depend on sign of velocity field

def calc_signs(
                indicesU = [],
                indicesV = [],
                lonsF = xr.DataArray(),
                latsF = xr.DataArray(),
                lonsTransect = xr.DataArray(),
                latsTransect = xr.DataArray(),
                transectType = "",
                dsMeshZ = xr.Dataset(),
                verboseLevel=False
               ):
    # extract indices of U,V segments from indices of F points
    if verboseLevel>0:
        print("calc_signs: calculate signs of segments")

    if len(indicesU) == 0:
        raise Exception("Error: calc_signs: input list has 0 length")

    signsU = []
    signsV = []


#               (positive normal for v)    
#                    ^     
#                    |   
#          f         v         f
#      j - .-------------------.       
#          |                   |
#          |                   |
#          |                   |
#      j   |u        o         |u ->  (positive normal for u) 
#          |        (j,i)      |
#          |                   |
#          |         v         |
#    j-1 - .-------------------.
#          f                   f
#          |         i         |       
#         i-1                  i
#


    # P1-P2 segment of transect line

    lamP1 = lonsTransect.isel(transect=0).item()
    phiP1 = latsTransect.isel(transect=0).item()

    lamP2 = lonsTransect.isel(transect=-1).item()
    phiP2 = latsTransect.isel(transect=-1).item()

    xP1Stereo,yP1Stereo = geo2stereo(lamP1,phiP1)
    xP2Stereo,yP2Stereo = geo2stereo(lamP2,phiP2)

    dXP1P2 = xP2Stereo-xP1Stereo
    dYP1P2 = yP2Stereo-yP1Stereo

    normTransect = np.array([dYP1P2,-dXP1P2])   # rotated 90 clockwise

    # PA-PB segment of either single U or V point

    for iSegment in range(len(indicesU)):
        xidxBack    = indicesU[iSegment][1]-1      # north-south segment, xidx is the same
        xidxForward = indicesU[iSegment][1]-1
        yidxBack    = indicesU[iSegment][0]-1
        yidxForward = indicesU[iSegment][0]
        lamPA = lonsF.isel(y=yidxBack   ,x=xidxBack   ) 
        lamPB = lonsF.isel(y=yidxForward,x=xidxForward) 
        phiPA = latsF.isel(y=yidxBack   ,x=xidxBack   ) 
        phiPB = latsF.isel(y=yidxForward,x=xidxForward) 
        # convert to stereographic
        xPAStereo,yPAStereo = geo2stereo(lamPA,phiPA)
        xPBStereo,yPBStereo = geo2stereo(lamPB,phiPB)
        dXAB = xPBStereo-xPAStereo
        dYAB = yPBStereo-yPAStereo
        # rotate 90deg clockwise
        dXABr     = dYAB
        dYABr     = -dXAB
        # do scalar product with transect segment
        normSegment = np.array([dXABr,dYABr])
        prod        = np.dot(normTransect,normSegment)
        signsU.append(np.sign(prod))


    for iSegment in range(len(indicesV)):
        xidxBack    = indicesV[iSegment][1]-1      # north-south segment, xidx is the same
        xidxForward = indicesV[iSegment][1]
        yidxBack    = indicesV[iSegment][0]
        yidxForward = indicesV[iSegment][0]
        lamPA = lonsF.isel(y=yidxBack   ,x=xidxBack   ) 
        lamPB = lonsF.isel(y=yidxForward,x=xidxForward) 
        phiPA = latsF.isel(y=yidxBack   ,x=xidxBack   ) 
        phiPB = latsF.isel(y=yidxForward,x=xidxForward) 
        # convert to stereographic
        xPAStereo,yPAStereo = geo2stereo(lamPA,phiPA)
        xPBStereo,yPBStereo = geo2stereo(lamPB,phiPB)
        dXAB = xPBStereo-xPAStereo
        dYAB = yPBStereo-yPAStereo
        # rotate 90deg anti-clockwise
        dXABr     = -dYAB
        dYABr     = dXAB 
        # do scalar product with transect segment
        normSegment = np.array([dXABr,dYABr])
        prod        = np.dot(normTransect,normSegment)
        signsV.append(np.sign(prod))

    signsU = np.array(signsU)
    signsV = np.array(signsV)

    if verboseLevel>3:
        print("signsU:")
        print(signsU)
        print("signsV:")
        print(signsV)

    nz = dsMeshZ.sizes["Z"]

    # expand dims (nUpoints) -> (z,nUpoints)
    #             (nVpoints) -> (z,nVpoints)  
    if transectType == "general":
        # expand dims for both signs
        signsU  = np.repeat(np.expand_dims(signsU,axis=0),repeats=nz,axis=0)
        signsV  = np.repeat(np.expand_dims(signsV,axis=0),repeats=nz,axis=0)
    if transectType == "pureX":
        signsV  = np.repeat(np.expand_dims(signsV,axis=0),repeats=nz,axis=0)
    if transectType == "pureY":            
        signsU  = np.repeat(np.expand_dims(signsU,axis=0),repeats=nz,axis=0)


    return signsU,signsV


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
