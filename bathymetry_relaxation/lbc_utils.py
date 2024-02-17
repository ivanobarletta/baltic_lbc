import xarray as xr

def get_nemo_transect_coords(ds,xidx=0,transect_name="transect_name"):
    # get the nav_lon,nav_lat from NEMO
    # LBC file (takes the outmost line (X=0))
    # returns 2 DataArrays to use for 
    # interpolation.    

    lons = ds["nav_lon"].data
    lats = ds["nav_lat"].data

    x_int = xr.DataArray(lons,dims=transect_name)
    y_int = xr.DataArray(lats,dims=transect_name)     

    return (x_int,y_int)