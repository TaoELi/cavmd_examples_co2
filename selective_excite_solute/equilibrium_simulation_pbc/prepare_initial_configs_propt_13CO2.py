#!/bin/bash

CavFreq=2320.0

for Nsub in 432 864 1728 3456
do

    # create a dictionary
    DIR=eq_N_$Nsub\_many13CO2
    mkdir $DIR
    cd $DIR
    
    # 1. copy input files from 216 co2
    cp ../template_eq/data_$Nsub.lmp data.lmp
    cp ../template_eq/in.lmp .
    cp ../template_eq/init_$Nsub.xyz init.xyz
    cp ../template_eq/input_equlibrate.xml .
    cp ../template_eq/input_traj.xml.bak .
    cp ../template_eq/create_photon_params.py .
    cp ../template_eq/submit_diabatical_incav.sh .
    
    # 2. change parameters for create_photon_params.sh create photon_params.json
    sed -i "s/216/$Nsub/" create_photon_params.py
    E0=$(echo "" | awk "END {print 4e-4 * (216.0 / $Nsub)^(0.5)}")
    echo "coupling is $E0"
    sed -i "s/0.0005/$E0/" create_photon_params.py
    python create_photon_params.py
    mv test.json photon_params.json

    # 3. Modify the last CO2 molecule to 13^CO2
    N13CO2=$(($Nsub/216 * 10))
    echo "Need to have $N13CO2 13CO2"
    NMOLECULE_START=$(($Nsub-$N13CO2))
    echo "$NMOLECULE_START molecules are CO2, rest are 13^CO2"
    
    NLINE_START=$(($NMOLECULE_START * 3 + 3))
    NLINE_END=$(($Nsub * 3 + 3))
    echo "editing from No.$NLINE_START line"
    sed -i "$NLINE_START,$(($NLINE_END))s/C/C13/g" init.xyz
    

    # 4. Change the simulation cell length for pbc calculations
    SIMU_CELL_LENGTH=$(echo "" | awk "END {print ($Nsub / 216.0)^(1.0/3.0) * 24.292}")
    sed -i "s/24.292/$SIMU_CELL_LENGTH/g" input_equlibrate.xml  

    # 5. change the file of submit jobs
    sed -i "s/Nsub=216/Nsub=$Nsub/" submit_diabatical_incav.sh 
    sed -i "s/one13CO2/many13CO2/" submit_diabatical_incav.sh 
    
    # 6. submit jobs to diabatical
    qsub submit_diabatical_incav.sh

    cd ..

done
