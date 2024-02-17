#!/bin/bash

module load cdo

VARIABLES=e1t,e2t,e3t_0

COORDINATES_PATH="../static_files/coordinates_east2.nc"
MESH_ZGR_PATH="../static_files/mesh_zgr_east2.nc"
TMASK_PATH="../static_files/mask_gridT_east2.nc"

cdo merge $COORDINATES_PATH $MESH_ZGR_PATH $TMASK_PATH merge.nc 

cdo -expr,"volume=e1t*e2t*e3t_0*mask" merge.nc volume.nc 

ncrename -O -d z_2,z volumeCDO.nc volumeCDO.nc

ncatted -O -a units,time,m,c,"seconds since 1900-01-01 00:00:00" volumeCDO.nc volumeCDO.nc 


