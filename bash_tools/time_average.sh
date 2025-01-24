#!/bin/bash

#SBATCH -J timemean
#SBATCH -o logs/logg%J.out
#SBATCH -e logs/logg%J.err
#SBATCH --time=01:00:00
#SBATCH --mem=45G

module load cdo/2.3.0

# create list of weeks

rootFolder="/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh"
fileType="3DV"
varName="vo"
fileType="3DU"
varName="uo"
fileNameRoot="NEATL36_TESTRUN_1d25h-m_"

fileList=${rootFolder}/R*/${fileNameRoot}${fileType}-${varName}_*-*.nc_ZNB

for file in ${fileList}; do
	ls -l $file
done

ncra ${fileList} mean_${fileNameRoot}${fileType}-${varName}_ZNB
