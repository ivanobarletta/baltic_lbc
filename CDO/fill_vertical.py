import xarray as xr
import numpy as np
import sys
from os.path import isfile

if len(sys.argv) < 3:
    print("Provide name of input.nc and output dataset!")
    sys.exit(0)

inPath = sys.argv[1]
outPath = sys.argv[2]

if not isfile(inPath):
    print("File: %s: does not exist" % inPath )
    sys.exit(0)


# Load your NetCDF dataset
# Replace 'your_dataset.nc' with the actual file path
ds = xr.open_dataset(inPath)

# Define a function to interpolate missing values along the vertical axis
def interpolate_vertical(ds_var):
    # Iterate over each time, lon, and lat point
    for time in ds_var.time:
        for lon in ds_var.lon:
            for lat in ds_var.lat:
                # Extract the data array for the current (time, lon, lat) point
                data_array = ds_var.sel(time=time, lon=lon, lat=lat)

                # Find the indices of missing values
                missing_indices = np.isnan(data_array)

                # Iterate over each depth level
                for depth in ds_var.depth:
                    # If the current depth level has a missing value
                    if missing_indices.sel(depth=depth):
                        # Find the closest non-missing value along the vertical axis
                        closest_depth = ds_var.sel(time=time, lon=lon, lat=lat, method='nearest', depth=depth, tolerance=0.1)
                        # Update the missing value with the closest non-missing value
                        data_array.loc[dict(depth=depth)] = closest_depth.values

                # Update the original dataset with the interpolated values
                ds_var.loc[dict(time=time, lon=lon, lat=lat)] = data_array

    return ds_var

# Apply the interpolation function to each variable in the dataset
for var_name in ds.data_vars:
    ds[var_name] = interpolate_vertical(ds[var_name])

# Save the updated dataset to a new NetCDF file
ds.to_netcdf(outPath)

