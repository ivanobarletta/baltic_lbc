#!/bin/bash

module load cdo
# no need to load NCO module (already imported with CDO)

METHOD="remapbic"

VARS_LIST=thetao,so,uo,vo,zos_detided

INPUT_DATASET="BAL-NEMO_PHY-DailyMeans-20220115.nc"
INPUT_DATASET_SSH="BAL-NEMO_PHY-detided_ssh-20220115.nc"
OUTPUT_DATASET="test_out.nc"

TARGET_HGRID_PATH="target_hgrid"
TARGET_VGRID_PATH="target_vgrid"
#TARGET_VGRID=0.4940254,1.541375,2.645669,3.819495,5.078224,6.440614,7.92956,9.572997,11.405,13.46714,15.81007,18.49556,21.59882,25.21141,29.44473,34.43415,40.34405,47.37369,55.76429,65.80727,77.85385,92.32607,109.7293,130.666,155.8507,186.1256,222.4752,266.0403,318.1274,380.213,453.9377,541.0889,643.5668,763.3331,902.3393,1062.44,1245.291,1452.251,1684.284,1941.893,2225.078,2533.336,2865.703,3220.82,3597.032,3992.484,4405.224,4833.291,5274.784,5727.917
#TARGET_VGRID="0.49,1.54,2.64,3.82,5.07"
LEVELS=$(cat ${TARGET_VGRID_PATH})

ANGLES_PATH="angles36_east_2.nc"

export EXTRAPOLATE=1
export REMAP_EXTRAPOLATE="on"

echo "Extapolation:  $REMAP_EXTRAPOLATE" 

#-f nc4

echo "INPUT DATASET:        " $INPUT_DATASET
echo "INPUT DATASET: (SSH): " $INPUT_DATASET_SSH

# merge the datasets
echo "Merging baltic 3D and Detided SSH..."
cdo -O merge $INPUT_DATASET $INPUT_DATASET_SSH merge.nc

# select only needed variables
echo "Selecting only needed variables..."
cdo select,name=${VARS_LIST} merge.nc select.nc

echo "Horizontal Interpolation..."
#cdo remapbic,${TARGET_HGRID_PATH} ${INPUT_DATASET} intp2D.nc 
#cdo remapdis,${TARGET_HGRID_PATH} ${INPUT_DATASET} intp2D.nc 
cdo ${METHOD},${TARGET_HGRID_PATH} select.nc intp2D.nc 

echo "Vertical Interpolation..."
#cdo intlevelx,level=${TARGET_VGRID} intp2D.nc intp3D.nc
cdo intlevelx,level=${LEVELS} intp2D.nc intp3D.nc

echo "Filling Nan with DWA (Distance Weighted Average)..."
cdo setmisstodis intp3D.nc intp3D_fill.nc

echo "Merging angles files into dataset"
cdo -O merge intp3D_fill.nc ${ANGLES_PATH} intp3D_fill_angles.nc 

# Velocity projection
# aexpr does the same as expr but retains other fields
echo "Projecting the velocity components along i,j..."
cdo -aexpr,"uo=(-1)*uo*gcost-vo*gsint;vo=uo*gsint-vo*gcost" intp3D_fill_angles.nc intp3D_fill_angles_proj.nc

# set reference time to hours since 1950-01-01..
cdo -setreftime,1950-01-01,00:00:00,1hour  intp3D_fill_angles_proj.nc ${OUTPUT_DATASET}

echo "Renaming variables and dimensions..."
# rename variables
ncrename -O -v thetao,votemper ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v so,vosaline ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v uo,vozocrtx ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v vo,vomecrty ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v zos_detided,sossheig ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v depth,deptht ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -v time,time_counter ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# rename dimensions
ncrename -O -d time,T ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -d depth,Z ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -d x,X ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncrename -O -d y,Y ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# rearrange dimension (T,Z,Y,X) -> (T,Z,X,Y)
#                     (T,  Y,X) -> (T,  X,Y)
echo "Re-arranging dimensions (T,Z,Y,X) -> (T,Z,X,Y)..."
ncpdq -O -a T,Z,X,Y ${OUTPUT_DATASET} ${OUTPUT_DATASET}

# add X,Y,Z,T variables
# X = 1,2,3,.....,NX
# Y = 1,2,3,.....,NY
ncap2 -A -s 'X[$X]=int(1)' -s 'X[$X]=array(1,1,$X)' ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncap2 -A -s 'Y[$Y]=int(1)' -s 'Y[$Y]=array(1,1,$Y)' ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncap2 -A -s 'Z[$Z]=int(1)' -s 'Z[$Z]=array(1,1,$Z)' ${OUTPUT_DATASET} ${OUTPUT_DATASET}
ncap2 -A -s 'T[$T]=int(1)' -s 'T[$T]=array(1,1,$T)' ${OUTPUT_DATASET} ${OUTPUT_DATASET}


# remove intermediate files
echo "Removing intermediate files..."
rm -f merge.nc select.nc intp2D.nc intp3D.nc intp3D_fill.nc intp3D_fill_angles.nc 
rm -f intp3D_fill_angles_proj.nc



