#!/bin/bash

NODE=n10

CavFreq=2320.0

for CONCENTRATION in 1.0 0.954
do

    # create a dictionary
    DIR=eq_c_$CONCENTRATION\_freq_$CavFreq
    mkdir $DIR
    cd $DIR
    
    # 1. copy input files from 216 co2
    cp ../co2_eq/data.lmp .
    cp ../co2_eq/in.lmp .
    cp ../co2_eq/init.xyz .
    cp ../co2_eq/input_equlibrate.xml .
    cp ../co2_eq/input_traj.xml.bak .
    cp ../co2_eq/photon_params.json .
    SUBMIT_SCRIPT=submit_diabatical_incav_4e-4.sh
    cp ../co2_eq/submit_diabatical_incav.sh $SUBMIT_SCRIPT
    
    # 2. change init.xyz so that some molecules are transferred to 13^CO2
    tmp=$(echo "$CONCENTRATION * 216" | bc)
    NMOLECULE_START=${tmp%.*}
    echo "$NMOLECULE_START molecules are CO2, rest are 13^CO2"
    NLINE_START=$(($NMOLECULE_START * 3 + 3))
    
    echo "editing from No.$NLINE_START line"
    
    sed -i "$NLINE_START,650s/C/C13/g" init.xyz
    
    # 3. change the file of submit jobs  
    sed -i "s/2e-4/4e-4/" $SUBMIT_SCRIPT
    sed -i "s/co2_eq/polariton_relaxation_biased\/equilibrium_simulation\/eq_c_$CONCENTRATION\_freq_$CavFreq/" $SUBMIT_SCRIPT 
    sed -i "s/RUN=eq/RUN=eq_c_$CONCENTRATION\_freq_$CavFreq/" $SUBMIT_SCRIPT
    sed -i "s/ipi_co2.eq/ipi_co2.eq_c_$CONCENTRATION\_freq_$CavFreq/" $SUBMIT_SCRIPT
    sed -i "s/n4/$NODE/" $SUBMIT_SCRIPT 
    
    # 4. change photon frequency 
    sed -i "s/2320.0/$CavFreq/" photon_params.json
    
    # 5. submit jobs to diabatical
    qsub $SUBMIT_SCRIPT

    cd ..

done
