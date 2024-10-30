import xarray as xr
import numpy as np
import sys

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

# the script creates 3D vertical coordinates using e3t 


# If cell k=7 is below sea bottom, depth_(k>=7) = depth_6

outFile     = "cmems_balmfc_depths.nc"

pathE3T     = "cmems_mod_bal_phy_anfc_static_e1t-e2t-e3t_9.04E-30.21E_53.01N-65.89N_0.50-712.02m.nc"
e3t         = xr.open_dataset(pathE3T)["e3t"]

# initialize with zeros
depths = xr.zeros_like(e3t)

# cumsum function 
for k in range(e3t.sizes["depth"]):
    if k == 0:
        depths.isel(depth=k)[:] = 0.5 * e3t.isel(depth=k)[:]
    else:
        depths.isel(depth=k)[:] = depths.isel(depth=k-1)[:] + 0.5 * e3t.isel(depth=k-1)[:] + 0.5 * e3t.isel(depth=k)[:]

# Recommended name
depths.name = "deptht" 

# save to netcdf
depths.to_netcdf(outFile)
