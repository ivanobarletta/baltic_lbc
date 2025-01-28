#!/bin/bash

# The script is a driver to 
# 	extract_from_archive.sh 
# the latter extracts files from tar.gz archives of NEMO files
# 

# the type of file extracted depends on the choice of 
# 1) fileNameRoot
#	examples:		
#	NEATL36_1d25h-m_ 	for 25h averaged files
#	NEATL36_1d-m_ 		for daily averaged files
#	NEATL36_6ts15mi 	for 15 minutes sampled files
#	NEATL36_1h-m_ 		for hourly files

# 2) fileType
#	examples:
#	3DT,2DT for variables on T cell
#	3DU,2DU for variables on U cell
#	3DV,2DV for variables on V cell
#
# 3) varName
#	examples:
#	so for salinity
#	thetao for temperature
#	uo,vo for 3D velocity components
#	oce for other variables	
#
# The complete filename is built by "extract_from_archive.sh" script.
# The general filename is like
#	{fileNameRoot}{fileType}-{varName}_YYYYmmDD-YYYYmmDD.nc

# NOTE!!!!!
# extract_from_archive.sh makes a 
# 	cdo selindexbox,xidx1,xidx2,yidx1,yidx2 inFile.nc reduced.nc
# where xidx1,xidx2,yidx1,yidx2 are defined within the script 
# if you don't want to subset just set the indexes to 
#	idxList="0,100000,0,100000"
#

rootFolder="/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/CONTROL"
dateInit="20211229"	# date of 1st simulated week
nWeeks=105
fileType="3DV"
varName="vo"
fileNameRoot="NEATL36_1d25h-m_"

fileType="2DT"
varName="oce"
fileNameRoot="NEATL36_1d-m_"

fileType="2DT"
varName="oce"
fileNameRoot="NEATL36_1d25h-m_"

# create list of weeks
for (( w = 0; w < ${nWeeks}; ++w )); do
	nDays=$(( 7*w ))		# days to add
	#echo "nDays: $nDays"
    	weekDate=$(date +%Y%m%d -d "$dateInit +${nDays} days")
	archivePath=${rootFolder}/R${weekDate}.tar.gz
	echo "Extracting files from: ${archivePath}"
	# count # of files in the folder. If < 7, submit the job
	nFiles=$(ls R${weekDate}/${fileNameRoot}${fileType}-${varName}*_ZNB | wc -l)	
	echo "Files already extraxted from archive: $nFiles"
	if [ ${nFiles} -lt 7 ]; then
		echo "   submitting the job to extract from archive"	
		sbatch extract_from_archive.sh ${weekDate} ${archivePath} ${fileNameRoot} ${fileType} ${varName}
	else
		echo "   Folder has all the files"
	fi
done

