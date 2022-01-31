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

# Define these important variables
#NMOLECULE=343
#LAMMP_DATA_FILE=data_343co2.lmp

# copy necessary file from outside
cp ../generate_lammps_inputs/$LAMMPS_DATA_FILE .
cp ../generate_lammps_inputs/init_raw_au.xyz .

#modify lammps input file
cp in.lmp in_opt.lmp
sed -i "s/data_343co2.lmp/$LAMMPS_DATA_FILE/" in_opt.lmp

# perform the i-pi minimization
cp input_optimize.xml input_opt_this.xml
sed -i "s/53.557/$SIMU_CELL_LENGTH/g" input_opt_this.xml
i-pi input_opt_this.xml &
# please change lmp_mpi to other lmp names (e.g., lmp)
sleep 30s
lmp_mpi < in_opt.lmp 


wait

# capture the optimized geometry
NLINE=$(($NMOLECULE*3+2))

tail -n $NLINE minimize.xc.xyz > init.xyz

sed -i.bak -e '1d' init.xyz
sed -i "1i $NLINE" init.xyz

# append two lines of cavity photons
echo "       L  0.00000e+00  0.00000e+00  0.00000e+00" >> init.xyz
echo "       L  0.00000e+00  0.00000e+00  0.00000e+00" >> init.xyz

# finally, remove irrelevent files
rm mini* RESTART log.* 
rm init.xyz.bak

# copy useful data to ../
cp $LAMMPS_DATA_FILE ../
cp in_opt.lmp ../
cp init.xyz ../

