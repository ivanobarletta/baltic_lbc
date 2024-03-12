#!/bin/bash

#SBATCH -J timeseries
#SBATCH -o logs/logg%J.out
#SBATCH -e logs/logg%J.err
#SBATCH --time=01:00:00
#SBATCH --mem=5G

#this script does not take into account masked
#value!!!
#use extract_ssh_cdo.sh instead!!!
#



module load cdo/2.3.0

rootFolder=/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL
rootFileName=NEATL36_1h-m_2DT-oce_*.nc 

fileList=${rootFolder}/*/${rootFileName}

lat=58.25
lon=9.98

lat=48
lon=-1.

outFile=control_run_ssh_lat_${lat}_lon_${lon}.dat

count=0

for file in ${fileList}; do
	count=$(($count + 1))


	echo $file
	if [ ${count} -gt 1 ]; then
		# select lines from the 2nd one (skip header)
		cdo --silent -outputtab,date,time,lat,lon,value -remapnn,lon=${lon}_lat=${lat} -select,name=ssh $file | tail -n +2 >> ${outFile}
	else 
		cdo --silent -outputtab,date,time,lat,lon,value -remapnn,lon=${lon}_lat=${lat} -select,name=ssh $file > ${outFile}
	fi
done
