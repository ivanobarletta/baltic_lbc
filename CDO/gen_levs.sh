#!/bin/bash

# Create a new netCDF file with dimensions (lon, lat)
ncgen -o target_levels.nc << EOF
netcdf example {
  dimensions:
    lev = 50;
  variables:
    float lev(lev);
}
EOF

TARGET_VGRID_PATH="target_vgrid"
#TARGET_VGRID=0.4940254,1.541375,2.645669,3.819495,5.078224,6.440614,7.92956,9.572997,11.405,13.46714,15.81007,18.49556,21.59882,25.21141,29.44473,34.43415,40.34405,47.37369,55.76429,65.80727,77.85385,92.32607,109.7293,130.666,155.8507,186.1256,222.4752,266.0403,318.1274,380.213,453.9377,541.0889,643.5668,763.3331,902.3393,1062.44,1245.291,1452.251,1684.284,1941.893,2225.078,2533.336,2865.703,3220.82,3597.032,3992.484,4405.224,4833.291,5274.784,5727.917
#TARGET_VGRID="0.49,1.54,2.64,3.82,5.07"
LEVELS=$(cat ${TARGET_VGRID_PATH})

# Set the values for the variable in the dataset
#ncap2 -s 'lev(:)={1,2,3,4,5}' target_levels.nc
#ncap2 -s 'lev={1,2,3,4,5}' target_levels.nc
ncap2 -s "lev={${LEVELS}}" target_levels.nc

# Optional: Check the content of the created netCDF file
ncdump target_levels.nc

