#!/bin/bash

#module load cdo
# no need to load NCO module (already imported with CDO)
#
# WARNING!!! 
# make sure that variables in VAR_LIST are present in the dataset
#


# method for horizontal interpolation
METHOD="remapcon"
METHOD="remapbic"

# path of cost,sint to project geographical velocity onto NEMO (i,j) system
ANGLES_PATH="angles36_east_2.nc"
ANGLES_PATH="../static_files/pyMakeAngles/new_angles_NEATL36_east2.nc"
ANGLES_PATH="../static_files/makeAngles/neatl36_angles_east2.nc"

VAR_SSH=zos_detided				# warning!!! glob has zos, not zos_detided
VAR_SSH=zos
# variables to interpolate
VARS_LIST=thetao,so,uo,vo,$VAR_SSH		# 
VARS_LIST=thetao,so,uo,vo,$VAR_SSH		# !!!! WARNING !!! glob has zos, not zos_detided

OUTPATH=outputs_glo

DATEYMD=$1
INPUT_DATASET=$2
INPUT_DATASET_SSH=$3
OUTPUT_DATASET=${OUTPATH}/NEATL36_east2_CDO_GLO_${DATEYMD}.nc
TEMPORARY_FILE="${DATEYMD}_tmp.nc"
TEMPORARY_FILE2="${DATEYMD}_tmp2.nc"
TEMPORARY_FILE3="${DATEYMD}_tmp3.nc"
BOX="13,14,54,55.6"

# target_hgrid2 is ordered in the same way as in NEATL36 OBC files
TARGET_HGRID_PATH="target_hgrid2"
TARGET_VGRID_PATH="target_vgrid"
# target_vgrid file must be like
# z1,z2,z3... with no spaces

#TARGET_VGRID=0.4940254,1.541375,2.645669,3.819495,5.078224,6.440614,7.92956,9.572997,11.405,13.46714,15.81007,18.49556,21.59882,25.21141,29.44473,34.43415,40.34405,47.37369,55.76429,65.80727,77.85385,92.32607,109.7293,130.666,155.8507,186.1256,222.4752,266.0403,318.1274,380.213,453.9377,541.0889,643.5668,763.3331,902.3393,1062.44,1245.291,1452.251,1684.284,1941.893,2225.078,2533.336,2865.703,3220.82,3597.032,3992.484,4405.224,4833.291,5274.784,5727.917
#TARGET_VGRID="0.49,1.54,2.64,3.82,5.07"
LEVELS=$(cat ${TARGET_VGRID_PATH})

### create nc file with target levels

#ncgen -o target_levels.nc << EOF
#netcdf example {
#    dimensions:
#      ilev = 50;
#    variables:
#      float lev(ilev);
#}
#EOF
#
#ncap2 -s "lev={${LEVELS}}" target_levels.nc

export EXTRAPOLATE=1
export REMAP_EXTRAPOLATE="on"

echo "Extapolation:  $REMAP_EXTRAPOLATE" 

echo "INPUT DATASET:        " $INPUT_DATASET
echo "INPUT DATASET: (SSH): " $INPUT_DATASET_SSH

# ------------------ Select a limited area
cdo sellonlatbox,${BOX} ${INPUT_DATASET} ${TEMPORARY_FILE}
cdo sellonlatbox,${BOX} ${INPUT_DATASET_SSH} ${TEMPORARY_FILE2}

# ------------------ merge the datasets
echo "Merging baltic 3D and Detided SSH..."
cdo -O merge ${TEMPORARY_FILE} ${TEMPORARY_FILE2} ${TEMPORARY_FILE3}
mv ${TEMPORARY_FILE3} ${TEMPORARY_FILE}

# ------------------ select only needed variables
echo "Selecting only needed variables..."
cdo select,name=${VARS_LIST} ${TEMPORARY_FILE} ${TEMPORARY_FILE2}
mv ${TEMPORARY_FILE2} ${TEMPORARY_FILE}

# ------------------ fill Nan with Distance Weighted Average
#echo "Filling Nan with DWA (Distance Weighted Average)..."
#cdo setmisstodis ${TEMPORARY_FILE} ${TEMPORARY_FILE2}
echo "Filling Nan with NN (Nearest Neighbour)..."
cdo setmisstonn ${TEMPORARY_FILE} ${TEMPORARY_FILE2}
mv ${TEMPORARY_FILE2} ${TEMPORARY_FILE}

# ------------------ fill NaN along vertical 
# Works ONLY with CDO Version > 1.9!! 
echo "Filling Nan along vertical with Nearest Method..."
cdo vertfillmiss,method=nearest ${TEMPORARY_FILE} ${TEMPORARY_FILE2} 
mv ${TEMPORARY_FILE2} ${TEMPORARY_FILE}

# ------------------ Horizontal Interpolation 
echo "Horizontal Interpolation..."
cdo ${METHOD},${TARGET_HGRID_PATH} ${TEMPORARY_FILE} ${TEMPORARY_FILE2} 
mv ${TEMPORARY_FILE2} ${TEMPORARY_FILE}

# ------------------ Vertical Interpolation
echo "Vertical Interpolation..."
cdo intlevelx,level=${LEVELS} ${TEMPORARY_FILE} ${TEMPORARY_FILE2}
mv ${TEMPORARY_FILE2} ${TEMPORARY_FILE}

# ------------------ Merging angles into dataset
echo "Merging angles files into dataset.."
cdo -O merge ${TEMPORARY_FILE} ${ANGLES_PATH} ${TEMPORARY_FILE2} 
mv ${TEMPORARY_FILE2} ${TEMPORARY_FILE}

# ------------------ Velocity projection
# aexpr does the same as expr but retains other fields
echo "Projecting the velocity components along i,j..."
cdo -aexpr,"uo=uo*gcost-vo*gsint;vo=uo*gsint+vo*gcost" ${TEMPORARY_FILE} ${TEMPORARY_FILE2}
mv ${TEMPORARY_FILE2} ${TEMPORARY_FILE}

# ------------------ set reference time to hours since 1950-01-01..
echo "Setting reference time to hours since 1950-01-01.."
cdo -setreftime,1950-01-01,00:00:00,1hour  ${TEMPORARY_FILE} ${OUTPUT_DATASET}

# ------------------ rename variables
echo "Renaming variables and dimensions..."
ncrename -O -v thetao,votemper ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v so,vosaline ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v uo,vozocrtx ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v vo,vomecrty ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v ${VAR_SSH},sossheig ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v depth,deptht ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v time,time_counter ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# ------------------ convert time_counter to float
ncap2 -O -s "time_counter=float(time_counter)" ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# ------------------ rename dimensions
ncrename -O -d time,T ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -d depth,Z ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -d x,X ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -d y,Y ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# ------------------ rearrange dimension (T,Z,Y,X) -> (T,Z,X,Y)
#                     (T,  Y,X) -> (T,  X,Y)
echo "Re-arranging dimensions (T,Z,Y,X) -> (T,Z,X,Y)..."
ncpdq -O -a T,Z,X,Y ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# ------------------ add X,Y,Z,T variables --------------------
# X = 1,2,3,.....,NX
# Y = 1,2,3,.....,NY
echo "Adding X,Y,Z,T coordinates"
ncap2 -A -s 'X[$X]=int(1)' -s 'X[$X]=array(1,1,$X)' ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncap2 -A -s 'Y[$Y]=int(1)' -s 'Y[$Y]=array(1,1,$Y)' ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncap2 -A -s 'Z[$Z]=int(1)' -s 'Z[$Z]=array(1,1,$Z)' ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncap2 -A -s 'T[$T]=int(1)' -s 'T[$T]=array(1,1,$T)' ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# ------------------ add attributes to coordinates -----------------
echo "Adding attributes to coordinates..."
# X  
ncatted -a standard_name,X,a,c,"projection_x_coordinate" ${OUTPUT_DATASET}
ncatted -a units,X,a,c,"1" ${OUTPUT_DATASET}
ncatted -a axis,X,a,c,"X" ${OUTPUT_DATASET}
ncatted -a grid_point,X,a,c,"T" ${OUTPUT_DATASET}
# Y  
ncatted -a standard_name,Y,a,c,"projection_y_coordinate" ${OUTPUT_DATASET}
ncatted -a units,Y,a,c,"1" ${OUTPUT_DATASET}
ncatted -a axis,Y,a,c,"Y" ${OUTPUT_DATASET}
ncatted -a grid_point,Y,a,c,"T" ${OUTPUT_DATASET}
# Z  
ncatted -a standard_name,Z,a,c,"projection_z_coordinate" ${OUTPUT_DATASET}
ncatted -a units,Z,a,c,"1" ${OUTPUT_DATASET}
ncatted -a axis,Z,a,c,"Z" ${OUTPUT_DATASET}
ncatted -a grid_point,Z,a,c,"T" ${OUTPUT_DATASET}
# T  
ncatted -a standard_name,T,a,c,"projection_t_coordinate" ${OUTPUT_DATASET}
ncatted -a units,T,a,c,"1" ${OUTPUT_DATASET}
ncatted -a axis,T,a,c,"T" ${OUTPUT_DATASET}
ncatted -a grid_point,T,a,c,"T" ${OUTPUT_DATASET}

# ------------------ remove unneeded variables 
# NOTE WELL! nav_lon,nav_lat are those taken from ${ANGLES_PATH}
# not necessarily those I need. I remove them and I rename
# lon,lat variables (the one where I did the interpolation) with nav_lon,nav_lat
echo "Removing Unneeded Variables (lon_bnds, gcost ...)" 
ncks -C -O -x -v gcost,gsint,lon_bnds,lat_bnds,nav_lon,nav_lat ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# ------------------- Here I rename lon,lat with nav_lon,nav_lat
ncrename -O -v lon,nav_lon ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v lat,nav_lat ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# ------------------- Declare Coordinates for variables
echo "Declaring coordinates for variables..."
#             attr    ,    var,mode,type ,attribute_value
ncatted -O -a coordinates,votemper,m,c,"nav_lon nav_lat" ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncatted -O -a coordinates,vosaline,m,c,"nav_lon nav_lat" ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncatted -O -a coordinates,vozocrtx,m,c,"nav_lon nav_lat" ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncatted -O -a coordinates,vomecrty,m,c,"nav_lon nav_lat" ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncatted -O -a coordinates,sossheig,m,c,"nav_lon nav_lat" ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# ------------------ remove intermediate files
echo "Removing intermediate files..."
#rm -f input_sel.nc input_sel_ssh.nc select_fill.nc
#rm -f merge.nc select.nc intp2D.nc intp3D.nc intp3D_fill.nc intp3D_fill_angles.nc 
#rm -f intp3D_fill_angles_proj.nc
if [ -f "${TEMPORARY_FILE}" ]; then
	rm ${TEMPORARY_FILE}
fi	
if [ -f "${TEMPORARY_FILE2}" ]; then
	rm ${TEMPORARY_FILE2}
fi	
if [ -f "${TEMPORARY_FILE3}" ]; then
	rm ${TEMPORARY_FILE3}
fi	


