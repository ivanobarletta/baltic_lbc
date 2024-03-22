#!/bin/bash

#SBATCH -J timeseries
#SBATCH -o logs/logg%J.out
#SBATCH -e logs/logg%J.err
#SBATCH --time=01:20:00
#SBATCH --mem=5G

#
#  this script take as argument (LAT,LON) coordinate to extract an hovemoller
#  from a list of files.
#  
#  check that the list of files is made of existing files and that varName is the 
#  same as in the dataset
#  
#  WARNING!!! provide coordinates in this order!!!
#  		1) LAT
#  		2) LON		
#  
#  


module load cdo/2.3.0

varName=so

rootFolder=/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL
rootFileName=NEATL36_1d25h-m_3DT-${varName}_*.nc

fileList=${rootFolder}/*/${rootFileName}

# source coordinates (unmasked) 
# I cannot use the inputfiles mask because their coordinates are partially masked
coordFile0="../target_grid/coordinates_NEATL36.nc"

if [ ! -f ${coordFile} ]; then
    echo "Error: coordinate file: ${coordFile} not found"
    exit 1
fi


# grab path of first file in the list
firstFile=$(ls ${fileList} | head -1)

# check arguments
if [[ $# -lt 2 ]] ; then
    echo "Error: not enough arguments provided"
    echo "provide (lat,lon) where to extract timeseries"
    exit 1
fi

lat0=$1
lon0=$2

# set up paths for temporary files / directories
prefix="${lat0}_${lon0}_${varName}"

tempDir="temp_${prefix}_ncks"
coordFile="${tempDir}/${prefix}_coords.nc"
distFile="${tempDir}/${prefix}_dist.nc"
maskFile="${tempDir}/mask.nc"

# create temporary dir to store temporary files

if [ ! -d "${tempDir}" ]; then
    mkdir "${tempDir}"
fi

outFile=control_run_${varName}_lat_${lat0}_lon_${lon0}.nc

echo "Extracting indexes of minimum distance from lat=${lat0} lon=${lon0}"
#-----------------------------------------------------------------------------------------------
# extract indexes from coords
ncks -O -v nav_lon,nav_lat ${coordFile0} ${coordFile}

# calculate mask from first file
# create mask from file 
cdo -expr,"mask=abs(${varName})/abs(${varName})" ${firstFile} ${maskFile}

# take surface mask
ncks -O -d deptht,0,0 ${maskFile} ${maskFile}

# add lat0,lon0 to coordinates
ncap2 -O -s "lon0=${lon0}" ${coordFile} ${coordFile}
ncap2 -O -s "lat0=${lat0}" ${coordFile} ${coordFile}

# compute distances
ncap2 -O -s 'distance=sqrt((nav_lon-lon0)*(nav_lon-lon0) + (nav_lat-lat0)*(nav_lat-lat0))' ${coordFile} ${distFile}

# multiply distances by the mask
cdo -O -mul ${maskFile} ${distFile} ${distFile}

# this puts the indexes of minimum distance in the attributes of the distance variable
ncap2 -O -s "distance@min=min_index(distance);print(distance@min)" ${distFile} ${distFile}

# extract indexes from attribute (that U is because the indexes are saved as idxULL
# I don't know what that ULL stands for...
yidx=$(ncdump -h ${distFile} | grep "distance:min" | cut -d "=" -f 2  | cut -d "," -f 1 | cut -d "U" -f 1)
xidx=$(ncdump -h ${distFile} | grep "distance:min" | cut -d "=" -f 2  | cut -d "," -f 2 | cut -d "U" -f 1)

# trim blank spaces
yidx=$(echo "${yidx}" | tr -d " ")
xidx=$(echo "${xidx}" | tr -d " ")

echo "I have found indexes   ${yidx} ${xidx}"
echo "for latitude longitude ${lat0} ${lon0}"

# remove intermediate files
#rm -f ${coordFile}
#rm -f ${distFile}
#rm -f ${maskFile}
#-----------------------------------------------------------------------------------------------


echo "Doing loop through file list"

count=0

numberLength=5
# format output numbers with 5 digits 98 -> 00098

#-----------------------------------------------------------------------------------------------
# File Loop
for file in ${fileList}; do
	count=$(($count + 1))
	formattedNumber=$(printf "%0${numberLength}d" "${count}")

	# Check existence of temporary file.
	# The queue might not be long enough so some files are not processed.
	# I retain the temporary folder so I can run again the script and extract
	# only the missing timeseries from non processed files

	if [ ! -f "${tempDir}/temp_${formattedNumber}.nc" ]; then
		echo "ncks -O -d x,${xidx},${xidx} -d y,${yidx},${yidx} -v ${varName} ${file} ${tempDir}/temp_${formattedNumber}.nc"
		ncks -O -d x,${xidx},${xidx} -d y,${yidx},${yidx} -v ${varName} ${file} ${tempDir}/temp_${formattedNumber}.nc 
	else
		echo "file: ${tempDir}/temp_${formattedNumber}.nc exists. Skipping"
	fi	
done
#-----------------------------------------------------------------------------------------------

echo "Concatenating files"
# concatenate all the files
ncrcat -O ${tempDir}/temp_*.nc ${outFile}

# remove temporary directory (if last command was ok)
if [ $? = 0 ]; then
	rm -rf ${tempDir}
fi

