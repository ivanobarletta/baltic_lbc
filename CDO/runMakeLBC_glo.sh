#!/bin/bash


#Purpose:
#	make a loop through cmems files from GLO_ANFC product. Files 3D
#	are separated from 2D (with only ssh). 
#	the folders are organized like:
#
#	ROOT/
#		DAILY/
#			YYYYmm/
#				file3D_YYYYmm01.nc
#				file3D_YYYYmm02.nc
#				file3D_YYYYmm03.nc
#				...	
#			YYYY(mm+1)/	
#				file3D_YYYY(mm+1)01.nc
#				file3D_YYYY(mm+1)02.nc
#				file3D_YYYY(mm+1)03.nc
#				...
#
#		SSH/
#			YYYYmm/
#				file2D_YYYYmm01.nc
#				file2D_YYYYmm02.nc
#				file2D_YYYYmm03.nc
#				... 2    
#			YYYY(mm+1)/ 2    
#				file2D_YYYY(mm+1)01.nc
#				file2D_YYYY(mm+1)02.nc
#				file2D_YYYY(mm+1)03.nc
#				...
#
#
#
#	loops are over years and months
#	The loops build a date (YYYYmm) and checks if the folder exists. 
#
#	If yes, the script does a further loop over the files contained in the
#	folder YYYYmm. If corresponding files  marked with date YYYYmmDD exist the 
#	the script (that takes YYYYmmD, path3D, path2D as arguments) is launched.
#	
# 	the run_job.sh is like:
#	#!/bin/bash
#	
#	#SBATCH -J test
#	#SBATCH -o logg%J.out
#	#SBATCH -e logg%J.err
#	#SBATCH --time=00:02:00
#	#SBATCH --mem=2G
#	
#	#module --force purge 
#	module load cdo/2.3.0
#	#module load nco/4.9.7
#	
#	arg1=$1
#	arg2=$2
#	arg3=$3
#	
#	echo "date:   ${arg1}"
#	echo "file3D: ${arg2}"
#	echo "fileSSH:${arg3}"
#	
#	bash createNEMOLBCwithCDO3.bash ${arg1} ${arg2} ${arg3}
#


ROOT=/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/GLOMFC_PROD

DAILY=DAILY
SSH=SSH

year1=2023
year2=2023

file3DRoot=GLO-NEMO-PHY-DailyMeans-
fileSSHRoot=GLO-NEMO-PHY-DailyMeans-SSH-
count=0

for y in $(seq ${year1} ${year2})	# loop through years
	do
	for m in `seq -f "%02g" 10 12`	# loop through months
		do
			dateYM=${y}${m}
			# check if folder exists
			folder3DPath=${ROOT}/${DAILY}/${dateYM}
			echo ${folder3DPath}
			if [ -d "${folder3DPath}" ]; then
				echo "folder exists"
				for file3D in ${folder3DPath}/*nc;
					do
					# grab day number from file name last but 5 and 4 digits
					echo "file 3D : ${file3D}"
					length=${#file3D}
					day="${file3D:length-5:1}${file3D:length-4:1}"
					#echo "day: $day"
					#fileSSH=${folder3DPath}/../../${SSH}/${dateYM}/${fileSSHRoot}${dateYM}${day}.nc
					fileSSH=${ROOT}/${SSH}/${dateYM}/${fileSSHRoot}${dateYM}${day}.nc
					dateYMD=${dateYM}${day}
					echo "file SSH: $fileSSH"
					if [ -f "${file3D}" ]; then
						echo "file 3D exists" 
						if [ -f "${fileSSH}" ]; then
							echo "both file exist"
							count=$(($count + 1))
							#sbatch run_job_test.job ${dateYMD} ${file3D} ${fileSSH}
							sbatch run_job2.job ${dateYMD} ${file3D} ${fileSSH}
						fi
					fi
					done	

			else
				echo "folder not present"
			fi

		done
	done


echo " "
echo "Processing $count Files"
echo " "
