import xarray as xr
import numpy as np

#  Example of vertical distribution of layers.
#  Layer 6 has partial cell thickness so depth_6 is shallower 
#  than in the full cell case.
#  
#      ---       
#       o    e3t_1      - depth_1
#      ---
#       o    e3t_2      - depth_2
#      ---
# 
#       o    e3t_3      - depth_3
#
#      ---
# 
#       o    e3t_4      - depth_4
#
#      ---
#
#
#       o    e3t_5      - depth_5
#
#
#      ---
#
#       o    e3t_6      - depth_6
#
#      ... 
#      /// 
#      ///
#      ---

# the script creates 3D vertical coordinates using e3t and 3D mask. The land points cell thickness is set to 0 and
# the resulting depths of land cells are the same as cells above.

# If cell k=7 is below sea bottom, depth_(k>=7) = depth_6

outFile     = "NEATL36_EAST2_testrun_depth3D.nc"

pathE3T     = "mesh_zgr.nc"

e3t         = xr.open_dataset(pathE3T)["e3t_0"]

# I set land point thickness to 0

# initialize with zeros
depths = xr.zeros_like(e3t)

# cumsum function 
for k in range(e3t.sizes["Z"]):
    if k == 0:
        depths.isel(Z=k)[:] = 0.5 * e3t.isel(Z=k)[:]
    else:
        depths.isel(Z=k)[:] = depths.isel(Z=k-1)[:] + 0.5 * e3t.isel(Z=k-1)[:] + 0.5 * e3t.isel(Z=k)[:]

# Recommended name
depths.name = "deptht" 

# save to netcdf
depths.to_netcdf(outFile)
