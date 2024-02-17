import xarray as xr
import os
import numpy as np
#import matplotlib.pyplot as plt
import sys
try:
    from cdo import Cdo
except:
    print("CDO is necessary to make horizontal Fill")
    sys.exit(0)    

def seaoverland(data, iterations=1, copy=False):
    
    if copy:
        data = np.ma.copy(data)

    if not isinstance(data, np.ma.masked_array) or not data.mask.any():
        return data

    for _ in range(iterations):
        shifted = []
        ni, nj = data.shape
        for i in range(-1, 2):
            for j in range(-1, 2):
                if i != 0 or j != 0:
                    # Shift the grid by i horizontally and j vertically and
                    # append it to the array. Shifted grids are 2 units smaller
                    # in both dimensions to accomodate the shift.
                    shifted.append(data[1 + i:ni - 1 + i, 1 + j:nj - 1 + j])

        # Calculate the mean value of the shifted grids to obtain the
        # approximated values. Only non-masked entries are taken into account.
        approx = np.ma.mean(shifted, axis=0)

        # Create a view without the outer points (so it is the same size as the
        # shifted grids), then copy the approximated values for the cells that
        # are masked.
        view = data[1:-1, 1:-1]
        np.copyto(view, approx, where=(view.mask & ~approx.mask))

        # Combine the two masks, unmasking values that were masked in view but
        # have been successfully approximated.
        view.mask &= approx.mask

    return data

def makeSOL(ds,nIters=5):
    # perform SeaOveLand for an xarray dataset

    fo = ds.copy()
    arr=fo.to_array()
    maArr=np.ma.masked_array(arr,mask=np.isnan(arr))

    # this for depth time lat lon dims
    n=np.array([seaoverland(tm,nIters)  for var in maArr for dpt in var for tm in dpt]).reshape(arr.shape)

    # this for  time lat lon dims
    #n=np.array([seaoverland(tm,5)  for var in maArr for tm in var ]).reshape(arr.shape)

    # replace dataset with filled variables
    count = 0 
    for var in fo.data_vars:
        #print(fo[var].data.shape)
        ndims = len(fo[var].data.shape)
        if ndims == 3:
            fo[var].data = np.expand_dims(n[count,0,0,:],axis=0)
        else:
            fo[var].data = n[count]
        count += 1

    return fo

def makeCDOSOL(ds):
    # do horizontal fill with CDO (I really cannot have it made by xarray...)
    cdo = Cdo()
    dsOut = cdo.setmisstonn(input=ds,returnXDataset=True)
    return dsOut

#fo.to_netcdf("SOL_%s" % fileIn)

