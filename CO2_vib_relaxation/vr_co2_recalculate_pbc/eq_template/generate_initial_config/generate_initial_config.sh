#!/bin/bash

while getopts n:f:w: option
do
case "${option}"
in
n) NMOLECULE=${OPTARG};;
f) LAMMPS_DATA_FILE=${OPTARG};;
w) SIMU_CELL_LENGTH=${OPTARG};;
esac
done

#NMOLECULE=343
#SIMU_CELL_LENGTH=53.557

cp create_config_packmol.inp job.inp

# change the initial configuration 
sed -i "s/343/$NMOLECULE/" job.inp
sed -i "s/53.557/$SIMU_CELL_LENGTH/g"  job.inp

packmol < job.inp

#rm job.inp
