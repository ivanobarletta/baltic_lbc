#!/bin/bash

#SBATCH -J tarextract
#SBATCH -o logs/logg%J.out
#SBATCH -e logs/logg%J.err
#SBATCH --time=03:00:00
#SBATCH --mem=45G

module load cdo/2.3.0

# extract files from a single archive
weekDate=$1
archivePath=$2
rootFile=$3		# NEATL36_1d25h-m_,NEATL36_1d-m_,...
fileType=$4		# 3DT/3DU/3DV...
varName=$5		# so/thetao/uo/vo variable name used for filename

idxList="850,1093,1300,1800"	# set 0,10000,0,100000 if you don't want to subset

echo "-----------------------"
echo "Archive: ${archivePath}"

for (( d = 0; d < 7; ++d )); do
    	# custom format using +
    	#date +%Y%m%d -d "$dateInit +$i days"
    	newDate=$(date +%Y%m%d -d "${weekDate} +${d} days")
	echo "   Extracting file with date: ${newDate}"
	fileName=${rootFile}${fileType}-${varName}_${newDate}-${newDate}.nc	
	archiveFolder="R${weekDate}"
	fileName2=${archiveFolder}/${fileName}
	echo "   File to extract: ${fileName2}"
	reducedFile=${fileName2}_ZNB
	if [ -f ${reducedFile} ]; then		# check existence of file
		echo "    File ${filename2} Exists"
	else
		echo "   Extracting File"		
		# Do the Extraction
		tar -xvf ${archivePath} ${fileName2}
		echo "   File Extracted"
		# Do the index selection
		cdo selindexbox,${idxList} ${fileName2} ${reducedFile}
		echo "   File Reduced"
		# remove the big file
		rm ${fileName2}
	fi
done

