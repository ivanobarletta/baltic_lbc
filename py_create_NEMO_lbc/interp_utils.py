import xarray as xr
from os.path import isfile
import sys
import matplotlib.pyplot as plt

def interp_variables(dictionary,lons2d,lats2d,zlevels,nml):
    # interpolate all datasets present in dictionary

    # create DataArray for zlevels

    #zlevs = xr.DataArray(data = zlevels, dims = "Z", coords=["deptht"])

    #zLevs = xr.DataArray(data = zlevels, dims = ["Z"], coords=dict(Z=(["Z",],zlevels)))
    zLevs = xr.DataArray(data = zlevels, dims = ["Z"])

    outDataArrays = []
    for iVar in dictionary:
        dsPath = dictionary[iVar]["path"]
        varName = dictionary[iVar]["name"]
        print("    opening dataset: %s" %dsPath)
        print("    for Variable   : %s" %varName)

        if isfile(dsPath):
            print ("    File present: %s" % dsPath)
        else:    
            raise Exception("    File not found: %s" %dsPath )
        
        # open datasets
        try:
            ds = xr.open_dataset(dsPath)
        except:
            print("Problem with Dataset %s" % dsPath)
            sys.exit(0)

        intpMethod = "nearest"
        intpMethod = "linear"

        kw={"fill_value": "extrapolate"}
        # do interpolation (2 steps for 3D vars because extrapolation is possible only in 1D)
        try:     
            #da_out = ds[varName].interp(lon=lons2d,lat=lats2d,depth=zLevs,method=intpMethod,kwargs=kwargs)
            #da_out = ds[varName].interp(lon=lons2d,lat=lats2d,depth=zLevs,method=intpMethod)
            da_out = ds[varName].interp(depth=zLevs,kwargs=kw).interp(lon=lons2d,lat=lats2d)
        except:
            # 2D variables    
            da_out = ds[varName].interp(lon=lons2d,lat=lats2d,method=intpMethod)
        ds.close()
        outDataArrays.append(da_out) 

    # merge datasets
    dsOut = xr.merge(outDataArrays)

    lnExtrap2D = nml["namvar"]["ln_extrap_2d"]

    #print ("---------")
    #print (dsOut)    
    #print ("---------")
    #print ("---------")
    #if lnExtrap2D:
    #    print("Extrapolating the dataset")
    #    dsOut = dsOut.interp(X=dsOut.X,Y=dsOut.Y,kwargs={"fill_value":"extrapolate"})
        
    # time intepolation 
    times = dsOut.time.data
    timeLevs = xr.DataArray(data = times, dims = ["T"])
    dsOut = dsOut.interp(time=timeLevs,method="nearest")

    # Define the new variable
    deptht = xr.DataArray(name="deptht", data=zlevels, dims=("Z"))
    time_counter = xr.DataArray(name="time_counter", data=times, dims=("T") )

    # Add the new variable to the dataset
    dsOut = dsOut.assign(deptht=deptht)
    dsOut = dsOut.assign(time_counter=time_counter)

    return dsOut

