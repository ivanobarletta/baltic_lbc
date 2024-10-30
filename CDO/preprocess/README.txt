# the scripts createDepth3DfromCMEMS.py create 
# a file that contains an array named depths3D(depth,y,x)

# to use this file in the vertical interpolation you must regrid
# to the same target mesh of the input file. So you have to remap
# like is done below (you can change to another remap method..)

# 	cdo remapcon,/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/CDO/target_hgrid2 cmems_balmfc_depths.nc balmfc_depth3D_on_NEATL36_EAST2.nc

# the cdo,tgtcoordinate takes a netCDF file that contains 1 single 3D array(Z,Y,X)
# so you must extract this single array with ncks -C 

#	ncks -C -v depth3D balmfc_depth3D_on_NEATL36_EAST2.nc 2balmfc_depth3D_on_NEATL36_EAST2.nc

