#!/bin/bash

for CONCENTRATION in 0.2 0.4 0.6 0.8
do

    # create a dictionary
    DIR=eq_c_$CONCENTRATION
    mkdir $DIR
    cd $DIR
    
    # 1. copy input files from 216 co2
    cp ../../co2_eq/data.lmp .
    cp ../../co2_eq/in.lmp .
    cp ../../co2_eq/init.xyz .
    cp ../../co2_eq/input_equlibrate.xml .
    cp ../../co2_eq/input_traj.xml.bak .
    cp ../../co2_eq/photon_params.json .
    cp ../../co2_eq/submit_diabatical_incav.sh .
    cp ../../co2_eq/submit_diabatical_outcav.sh .
    
    # 2. change init.xyz so that some molecules are transferred to 14C_18O2
    tmp=$(echo "$CONCENTRATION * 216" | bc)
    NMOLECULE_START=${tmp%.*}
    NLINE_START=$(($NMOLECULE_START * 3 + 3))
    
    echo "editing from No.$NLINE_START line"
    
    sed -i "$NLINE_START,650s/C/C14/g" init.xyz
    sed -i "$NLINE_START,650s/O/O18/g" init.xyz
    
    
    # 3. change the file of submit jobs
    sed -i "s/co2_eq/vr_co2_concentration\/eq_c_$CONCENTRATION/" submit_diabatical_incav.sh 
    sed -i "s/RUN=eq/RUN=eq_c_$CONCENTRATION/" submit_diabatical_incav.sh 
    sed -i "s/ipi_co2.eq/ipi_co2.eq_c_$CONCENTRATION/" submit_diabatical_incav.sh 
    
    sed -i "s/co2_eq/vr_co2_concentration\/eq_c_$CONCENTRATION/" submit_diabatical_outcav.sh 
    sed -i "s/RUN=eq/RUN=eq_c_$CONCENTRATION/" submit_diabatical_outcav.sh 
    sed -i "s/ipi_co2.eq/ipi_co2.eq_c_$CONCENTRATION/" submit_diabatical_outcav.sh 
    
    cd ..

done
