#!/bin/bash



Cd=0.001	# drag coefficient (see Curcic - Revised Estimates of Ocean Surface Drag in Strong Winds,2020)
rhoAir=1.2	# air density [kg * m^-3]

# taux = Cd * rhoAir * |U| * Ux [kg * m^-1 * s^-2] = [N * m^-2]
# tauy = Cd * rhoAir * |U| * Uy [kg * m^-1 * s^-2] = [N * m^-2]

cdo --version

if [ $? -ne 0 ]; then
        echo "Error: CDO not loaded"
	exit
fi

fileNameU=$1	#file with u-velocity
fileNameV=$2	#file with v-velocity

# check existence
if [ ! -f "${fileNameU}" ]; then
	echo "Error: file ${fileNameU} not existing"
	exit
fi

if [ ! -f "${fileNameV}" ]; then
	echo "Error: file ${fileNameU} not existing"
	exit
fi

fileNameUV=$(basename "${fileNameU}" | sed 's/BULKU10M/fileUV/')
echo $fileNameUV
fileNameUVmod=$(basename "${fileNameU}" | sed 's/BULKU10M/fileUVmod/')
echo $fileNameUVmod
fileNameUV2=$(basename "${fileNameU}" | sed 's/BULKU10M/fileUV2/')
echo $fileNameUV2
fileNameTau=$(basename "${fileNameU}" | sed 's/BULKU10M/fileTau/')
echo $fileNameTau

cdo -O merge ${fileNameU} ${fileNameV} ${fileNameUV}

# calc |U| 
cdo -O -expr,'uvmod=sqrt(sowinu10*sowinu10 + sowinv10*sowinv10)' ${fileNameUV} ${fileNameUVmod}
# merge files
cdo -O merge ${fileNameUV} ${fileNameUVmod} ${fileNameUV2}
# calc stress components
cdo -O -expr,"utau=${Cd}*${rhoAir}*uvmod*sowinu10;vtau=${Cd}*${rhoAir}*uvmod*sowinv10;qtot=0*uvmod;qsr=0*uvmod;emp=0*uvmod" ${fileNameUV2} ${fileNameTau}

ncatted -a units,utau,a,c,"N * m^-2" ${fileNameTau}
ncatted -a units,vtau,a,c,"N * m^-2" ${fileNameTau}

# remove temporary files
#rm ${fileNameUV}
#rm ${fileNameUV2}
#rm ${fileNameUVmod}

# create symbolic link with NEMO convention name (y0000m00d00)

DATE=$(echo "${fileNameTau}" | cut -d"_" -f 4)
year="${DATE:0:4}"  # Extract year (first 4 characters)
month="${DATE:4:2}"  # Extract month (characters 5 and 6)
day="${DATE:6:2}"  # Extract day (characters 7 and 8)

# Concatenate the strings with the additional characters
FORMATTED_DATE="y${year}m${month}d${day}"

fileNameTau2=filetau_G1280_${FORMATTED_DATE}.nc

# create link with right filename
echo "Filename Tau2: ${fileNameTau2}"
ln -s ${fileNameTau} ${fileNameTau2}



