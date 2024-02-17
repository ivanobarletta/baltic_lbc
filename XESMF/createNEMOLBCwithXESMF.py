import xarray as xr 
import xesmf as xe 
import numpy as np
import matplotlib.pyplot as plt
import sys
from utils import slice_dataset,interp_vertical
from utils import rename_dataset,reshape_dataset,project_velocity
from utils import change_time_reference, add_TZXY, change_type
from sol_ocean import makeSOL,makeCDOSOL
try:
    from cdo import Cdo
except:
    print("CDO is necessary to make vertical Fill")
    sys.exit(0)
    

fileIN = "../CDO/BAL-NEMO_PHY-DailyMeans-20220115.nc"
fileSSH = "../CDO/BAL-NEMO_PHY-detided_ssh-20220115.nc"
target_hgrid = "target_hgrid.nc"
target_vgrid = "target_vgrid"
angles_path = "../CDO/angles36_east_2.nc"
method = "conservative"
outPath = "test_out_xesmf.nc"
variables = ["thetao","so","uo","vo"]
SOL = True                              # use SeaOverLand routine to fill dataset
l_CDOSOL = True             # use CDO for horizontal fill
l_dbg = False

# Possible errors:
# 
# Error after applying 3D mask to input dataset:
#   ds["mask"] = xr.where(np.isnan(ds["thetao"].isel(time=0)),1,0)
# ValueError: mask must have the same shape as the 
# latitude/longitudecoordinates, got: mask.shape = (763, 774, 56), lon.shape = (763, 774)
# 
# the mask can only be 2D

# Hints 
# to extrapolate is necessary to skip NaN when using regridder
#
#   ds_out = regridder_extrap(ds,skipna=True)
# 

# Creation of destination grid 
#grid_out = xr.Dataset(
#    {
#        "lat": (["lat"], np.arange(16, 75, 1.0), {"units": "degrees_north"}),
#        "lon": (["lon"], np.arange(200, 330, 1.5), {"units": "degrees_east"}),
#    }
#)

# How to add a 2D mask to a dataset
# ds["mask"] = xr.where(np.isnan(ds["varname"].isel(time=0,depth=0)),0,1)


def main(inputPath=fileIN,
            targetHGridPath=target_hgrid,
            targetVGridPath=target_vgrid,
            method=method,
            fillMethod="nearest_s2d",
            outFile=outPath):


    # open datasets
    print("Opening Datasets:")
    dsIn = xr.open_dataset(inputPath)
    # select only needed variables
    dsIn = dsIn[variables]
    # open ssh file    
    dsSSH = xr.open_dataset(fileSSH)

    # merge ssh from other file
    print("Merging with SSH dataset")
    dsIn = xr.merge([dsIn,dsSSH])
    # open target grid
    targetHGrid = xr.open_dataset(targetHGridPath)

    # create box to slice input dataset to a smaller region
    print("slicing Input Dataset")
    box = [targetHGrid.lon.min(),targetHGrid.lon.max(),targetHGrid.lat.min(),targetHGrid.lat.max()]
    dsIn = slice_dataset(dsIn,box)

    if l_dbg:
        dsIn.to_netcdf("dsSlice.nc")

    if not SOL:
        # use esmf to fill Nan horizontally (does not work...) 
        # add 2D mask to dataset
        print("Adding 2D mask to dataset")
        dsIn["mask"] = xr.where(np.isnan(dsIn["thetao"].isel(time=0,depth=0)),0,1)

        # create grid to fill the dataset (same as Input Dataset)
        
        gridFill = xr.Dataset({
            "lat": (["lat"], dsIn.lat.values, {"units": "degrees_north","standard_name":"latitude"}),
            "lon": (["lon"], dsIn.lon.values, {"units": "degrees_east","standard_name":"longitude"}),
                                })
        
        """
        lon2d, lat2d = np.meshgrid(dsIn.lon.values,dsIn.lat.values)
        gridFill = xr.Dataset({
            "lat": (["y","x"], lat2d, {"units": "degrees_north","standard_name":"latitude"}),
            "lon": (["y","x"], lon2d, {"units": "degrees_east","standard_name":"longitude"}),
                                })
        """
        # create a regridder to fill the input dataset
        print("Creating regridder to fill")
        regridder_fill = xe.Regridder(dsIn,gridFill,method="nearest_s2d")
        print("Filling Dataset")
        #dsFill = regridder_fill(dsIn,skipna=True,na_thres=0)
        dsFill = regridder_fill(dsIn,skipna=False,na_thres=0)

    else:
        if l_CDOSOL:
            print("Makig SeaOverLand on dataset with CDO")
            dsFill = makeCDOSOL(dsIn)
        else:    
            print("Makig SeaOverLand on dataset")     
            dsFill = makeSOL(dsIn)

    if l_dbg:
        dsFill.to_netcdf("dsFill.nc")    
    #sys.exit(0)

    # create regridder object
    print("Creating regridder to interpolate")
    regridder_cons = xe.Regridder(dsFill, targetHGrid, method=method)
    print("Interpolating Dataset")
    dsOut = regridder_cons(dsFill)

    # vertical inerpolation
    print("Doing Vertical interpolation")
    if l_dbg:
        dsOut.to_netcdf("before_vinterpolation.nc")
    dsOut = interp_vertical(dsOut,target_vgrid)

    # do vertical fill with CDO (I really cannot have it made by xarray...)
    cdo = Cdo()
    dsOut = cdo.vertfillmiss(input=dsOut,returnXDataset=True)

    if l_dbg:
        dsOut.to_netcdf("dsVFill.nc")

    # project velocities
    print("Projecting Velocities")
    dsOut = project_velocity(dsOut,xr.open_dataset(angles_path))

    if l_dbg:
        dsOut.to_netcdf("after_proj.nc")
    # replace coordinates (For some reasone after CDO verfillmiss nav_lon,nav_lat
    # disappear ad re-appear after projection of velocity. Mistery...)
    # I use the nav_lon,nav_lat from target grid because are double
    # those from angles_path should be the same but are real        
    #dsOut["nav_lon"] = np.array(targetHGrid.lon.data,dtype="float64")
    #dsOut["nav_lon"] = np.array(targetHGrid.lat.data,dtype="float64")
    dsOut = dsOut.drop_vars(["nav_lon","nav_lat"])
    dsOut["nav_lon"] = targetHGrid.lon
    dsOut["nav_lat"] = targetHGrid.lat
    dsOut = dsOut.drop_vars(["lon","lat"])


    # rename dims/vars and reshape (t,z,y,x) -> (t,z,x,y)
    print("Transposing dimensions")
    dsOut = reshape_dataset(dsOut)
    print("renaming Variables(Dimensions)")
    dsOut = rename_dataset(dsOut)

    dsOut = change_time_reference(dsOut)

    dsOut = add_TZXY(dsOut)

    dsOut = change_type(dsOut)

    # print to file
    print("Saving to File: %s" %outFile)
    dsOut.to_netcdf(outFile)
    
#    plt.figure()
#    dsOut["thetao"].isel(time=0,depth=0).plot(x="lon",y="lat",cmap="jet",vmin=0,vmax=5)
#    plt.show()


if __name__ == "__main__":
    main(inputPath=fileIN,
            targetHGridPath=target_hgrid,
            targetVGridPath=target_vgrid,
            method=method,
            outFile=outPath)