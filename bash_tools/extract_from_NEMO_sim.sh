#!/bin/bash

#SBATCH -J subset
#SBATCH -o logs/logg%J.out
#SBATCH -e logs/logg%J.err
#SBATCH --time=01:55:00
#SBATCH --mem=10G

# simple script to subset files from a NEMO simulation





module load cdo/2.3.0

year=2023

rootFolder="/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/RAWDATA_NEATL36/PHY/PHY-FC-FRE/TESTRUN2"

varName="so"
rootFile="NEATL36_TESTRUN_1d25h-m_3DT-${varName}_${year}*-*.nc"
varName="uo"
rootFile="NEATL36_TESTRUN_1d25h-m_3DU-${varName}_${year}*-*.nc"
varName="thetao"
rootFile="NEATL36_TESTRUN_1d25h-m_3DT-${varName}_${year}*-*.nc"
varName="vo"
rootFile="NEATL36_TESTRUN_1d25h-m_3DV-${varName}_${year}*-*.nc"

varName="oce"
rootFile="NEATL36_TESTRUN_1h-m_2DU-${varName}_${year}*-*.nc"


idxList="600,1092,1200,1893"
idxList="700,1092,1300,1893"
idxList="850,1092,1350,1780"
idxList="850,1093,1300,1800"

fileList=${rootFolder}/R*/${rootFile}

for file in ${fileList}; do
	ls -l $file
	# it's essential to grab the basename of file (to write output in pwd..)
	file2=$(basename $file)
	# make selection and save in current directory 	
	echo ${file2}
	if [ ! -f "${file2}_ZNB" ]; then	# do this only if temporary file does not exist
		cdo -O selindexbox,${idxList} ${file} ${file2}_ZNB
	else
		echo "temporary file is already existing"
	fi
done

# concatenate
#ncrcat -O ${varName}/cut_*nc concat_${varName}_testrun2_${year}.nc

# remove temporary files
#rm ${varName}/cut_*nc 

