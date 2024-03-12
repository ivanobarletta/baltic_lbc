#!/bin/bash

#this is a test script to estimante the indexes
#of closest mesh point to a given (lat,lon) taking
#into account also of land masked values


sampleFile="sample_ssh.nc"				# source file
coordsFile="../target_grid/coordinates_NEATL36.nc" 	# unmasked coordinates

varName=ssh

# check arguments
if [[ $# -lt 2 ]] ; then
    echo "Error: not enough arguments provided"
    echo "provide (lat,lon) to calculate closest"
    echo "indexes"
    exit 1
fi

lat0=$1
lon0=$2

echo "extracting indices for lat=${lat0} lon=${lon0}"

# extract coords
ncks -O -v nav_lon,nav_lat ${coordsFile} coordinates.nc

# create mask from file 
cdo -expr,"mask=abs(${varName})/abs(${varName})" ${sampleFile} mask.nc

# add lat0,lon0 to coordinates
ncap2 -O -s "lon0=${lon0}" coordinates.nc coordinates.nc
ncap2 -O -s "lat0=${lat0}" coordinates.nc coordinates.nc

# compute distances
ncap2 -O -s 'distance=sqrt((nav_lon-lon0)*(nav_lon-lon0) + (nav_lat-lat0)*(nav_lat-lat0))' coordinates.nc distances.nc

# multiply distances by the mask
cdo -O -mul mask.nc distances.nc distances.nc
## append mask to distances
#ncks -A -v mask mask.nc distances.nc 
#cdo -expr,"distance=mask*distance" distances.nc distances2.nc

ncap2 -O -s "distance@min=min_index(distance);print(distance@min)" distances.nc distances.nc

yidx=$(ncdump -h distances.nc | grep "distance:min" | cut -d "=" -f 2  | cut -d "," -f 1 | cut -d "U" -f 1)
xidx=$(ncdump -h distances.nc | grep "distance:min" | cut -d "=" -f 2  | cut -d "," -f 2 | cut -d "U" -f 1)

# trim blank spaces
yidx=$(echo "${yidx}" | tr -d " ")
xidx=$(echo "${xidx}" | tr -d " ")

echo "I have found indexes   ${yidx} ${xidx}"
echo "for latitude longitude ${lat0} ${lon0}"

# remove intermediate files
#rm -f coordinates.nc
#rm -f distances.nc

ncks -O -d x,${xidx},${xidx} -d y,${yidx},${yidx} -v ${varName} ${sampleFile} test_out.nc
#echo "ncks -O -d x,${xidx},${xidx} -d y,${yidx},${yidx} -v ${varName} ${sampleFile} test.nc"


