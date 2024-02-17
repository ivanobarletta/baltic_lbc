#!/bin/bash



ROOT=/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS
ROOT=/mnt/lustre/scratch/nlsas/home/empresa/nrd/NRD/STORE/BALMFC_PRODUCTS/LBC_4NEATL/BALMFC_PROD	# here I have 2 years data

DAILY=DAILY
SSH=SSH_DETIDED

year1=2022
year2=2023

file3DRoot=BAL-NEMO_PHY-DailyMeans-
fileSSHRoot=BAL-NEMO_PHY-detided_ssh-

configFile="config.txt"

echo "ArrayTaskID          dateYYYYMM            		path3D                          	path2D" > $configFile

countID=0

for y in $(seq ${year1} ${year2})	# loop through years
	do
	for m in `seq -f "%02g" 1 12`	# loop through months
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
							#sbatch run_job_test.job ${dateYMD} ${file3D} ${fileSSH}
							#sbatch run_job2.job ${dateYMD} ${file3D} ${fileSSH}
							countID=$(($countID + 1))
							echo "countID: $countID"
							echo "$countID ${dateYMD} ${file3D} ${fileSSH}" >> $configFile
						fi
					fi
					done	

			else
				echo "folder not present"
			fi

		done
	done

	
