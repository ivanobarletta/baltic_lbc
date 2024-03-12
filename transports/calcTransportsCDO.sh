#!/bin/bash

#SBATCH -J transpCDO
#SBATCH -o logs/logg%J.out
#SBATCH -e logs/logg%J.err
#SBATCH --time=00:40:00
#SBATCH --mem=10G

#
# this script is still under development!!!
#
#

module load cdo/2.3.0

rootFolder=/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL
rootFileName=NEATL36_1d25h-m_3DU-uo_*.nc

fileList="${rootFolder}/*/${rootFileName}"

#echo ${fileList}

varName=uo

#idxList="517,517,400,425"
idxList="$1,$2,$3,$4"

# check arguments
if [[ $# -lt 4 ]] ; then
    echo "Error: not enough arguments provided"
    echo "provide selindexbox indices"
    exit 1
fi

coordsPath="../target_grid/coordinates.nc"
zMeshPath="../static_files/mesh_zgr.nc"

# same list but with underscores
idxList2=$(echo "${idxList}" | sed "s/,/_/g")

outFile="transports_CDO_${idxList2}.nc"

tempDir="temp_${idxList2}"

if [ ! -d "${tempDir}" ]; then
    mkdir "${tempDir}"
fi


# slice scale factor files
cdo -O selindexbox,${idxList} ${coordsPath} ${tempDir}/box_coords.nc
cdo -O selindexbox,${idxList} ${zMeshPath} ${tempDir}/box_mesh_zgr.nc
# merge in a single file
cdo -O merge ${tempDir}/box_mesh_zgr.nc ${tempDir}/box_coords.nc ${tempDir}/merge.nc
# calculate faces
cdo -O -expr,"faces=e2u*e3u_0" ${tempDir}/merge.nc ${tempDir}/faces.nc

numberLength=5

count=0
for file in ${fileList}; do 
	echo $file; 
	count=$(($count + 1))
	formattedNumber=$(printf "%0${numberLength}d" "${count}")

	echo "cdo selindexbox,$idxList $file ${tempDir}/temp_${formattedNumber}.nc"
	cdo -O selindexbox,$idxList $file ${tempDir}/temp_${formattedNumber}.nc
	# rename depthu dimension -> z 
	ncrename -O -d depthu,z ${tempDir}/temp_${formattedNumber}.nc ${tempDir}/temp_${formattedNumber}.nc
	# merge coords and vel data
	echo "cdo merge ${tempDir}/temp_${formattedNumber}.nc ${tempDir}/faces.nc ${tempDir}/temp2_${formattedNumber}.nc"
	cdo -O merge ${tempDir}/temp_${formattedNumber}.nc ${tempDir}/faces.nc ${tempDir}/temp2_${formattedNumber}.nc
	# calculate transport across single cells
	echo "cdo -expr,"uo_dA=uo*faces" ${tempDir}/temp2_${formattedNumber}.nc ${tempDir}/temp_uo_dA_${formattedNumber}.nc"
	cdo -O -expr,"uo_dA=uo*faces" ${tempDir}/temp2_${formattedNumber}.nc ${tempDir}/temp_uo_dA_${formattedNumber}.nc
	# use ncwa to calculate the sum of d(Transp) along depth,y dimension
	# the -N option prevents ncwa to divide by the denominator that normalizes the integral
	echo "ncwa -O -N -a z,y ${tempDir}/temp_uo_dA_${formattedNumber}.nc ${tempDir}/temp_transp_${formattedNumber}.nc"
	ncwa -O -N -a z,y ${tempDir}/temp_uo_dA_${formattedNumber}.nc ${tempDir}/temp_transp_${formattedNumber}.nc
	# remove x dimension
	ncwa -O -a x ${tempDir}/temp_transp_${formattedNumber}.nc ${tempDir}/temp_transp_${formattedNumber}.nc
done

# concatenate files
ncrcat -O ${tempDir}/temp_transp_*.nc ${outFile}
# rename variable
ncrename -O -v uo_dA,transport ${outFile} ${outFile}
# add units
ncatted -O -a units,transport,c,c,"m**3/s" ${outFile} ${outFile}

# remove temporary directory (if last command was ok)
if [ $? = 0 ]; then
        rm -rf ${tempDir}
fi

