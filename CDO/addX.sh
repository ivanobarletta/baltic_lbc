#!/bin/bash

rm prova2.nc
#ncap2 -s 'X[$X]=int(1)' -s 'X(:)={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}' prova.nc prova2.nc


#ncap2 -s 'X[$X]=int(1)' prova.nc prova2.nc

ncap2 -A -s 'X[$X]=int(1)' -s 'X[$X]=array(1,1,$X)' prova.nc prova.nc 


#ncap2 -s 'X[$X]=int(1)' -s 'X(:)={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}' prova.nc prova2.nc
