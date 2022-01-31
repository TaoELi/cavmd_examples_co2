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
#LAMMPS_DATA_FILE=data_343co2.lmp
#SIMU_CELL_LENGTH=53.557

NANGLE=$NMOLECULE
NBOND=$(($NMOLECULE*2))
NATOM=$(($NMOLECULE*3))

# 1. prepare files
cp ../generate_initial_config/init_raw_au.xyz .
cp bak.lt job.lt
sed -i "s/53.557/$SIMU_CELL_LENGTH/g" job.lt
sed -i "s/343/$NMOLECULE/g" job.lt

# 2. create a lmp data file
moltemplate.sh -xyz init_raw_au.xyz job.lt

# 3. clean not important files
rm -rf output_ttree/
rm job.in*
mv job.data $LAMMPS_DATA_FILE

# 4. modify the new lmp data file

# 4.1. clean the original data header
sed -i.bak -e '3,13d' $LAMMPS_DATA_FILE

# 4.2. add the new data header
INFO_ATOMS=" $NATOM  atoms\n $NBOND  bonds\n $NANGLE  angles\n\n 2  atom types\n 1  bond types\n 1  angle types\n"
INFO_FF="Masses\n\n  1 12.0107\n  2 15.9994\n\nBond Coeffs\n\n 1    2.196    0.45106005 -0.67297933 0.5857145\n\nAngle Coeffs\n\n  1    0.0861  180.000000\n"

sed -i.bak "3i  $INFO_ATOMS" $LAMMPS_DATA_FILE
sed -i.bak "16i  $INFO_FF" $LAMMPS_DATA_FILE

# 4.3. clean the backup file
rm $LAMMPS_DATA_FILE.bak 
