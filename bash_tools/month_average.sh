#!/bin/bash

#SBATCH -J timemean
#SBATCH -o logs/logg%J.out
#SBATCH -e logs/logg%J.err
#SBATCH --time=00:10:00
#SBATCH --mem=5G

module load cdo/2.3.0

year=$1
month=$2

# create list of weeks

rootFolder="/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/TESTRUN2_run_outs/ZNB_native_mesh"
rootFolder="."

fileType="3DV"
varName="vo"

fileType="2DU"
varName="oce"
fileNameRoot="NEATL36_TESTRUN_1h-m_"

#NEATL36_TESTRUN_1h-m_2DU-oce_20221221-20221221.nc_ZNB
#NEATL36_TESTRUN_1h-m_2DU-oce_202205.nc_ZNB

fileList=${fileNameRoot}${fileType}-${varName}_${year}${month}??-????????.nc_ZNB

for file in ${fileList}; do
	ls -l $file
done

ncra -O ${fileList} ${fileNameRoot}${fileType}-${varName}_${year}${month}.nc_ZNB

rm ${fileList}
