#!/bin/bash

NODE=n10

CavFreq=2320.0

for CONCENTRATION in 0.0 0.954 1.0
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
    cp ../co2_eq/submit_diabatical_incav.sh .
    cp ../co2_eq/submit_diabatical_outcav.sh .
    
    # 2. change init.xyz so that some molecules are transferred to 13^CO2
    tmp=$(echo "$CONCENTRATION * 216" | bc)
    NMOLECULE_START=${tmp%.*}
    echo "$NMOLECULE_START molecules are CO2, rest are 13^CO2"
    NLINE_START=$(($NMOLECULE_START * 3 + 3))
    
    echo "editing from No.$NLINE_START line"
    
    sed -i "$NLINE_START,650s/C/C13/g" init.xyz
    #sed -i "$NLINE_START,650s/O/O18/g" init.xyz
    
    
    # 3. change the file of submit jobs
    sed -i "s/co2_eq/polariton_relaxation_biased\/equilibrium_simulation\/eq_c_$CONCENTRATION\_freq_$CavFreq/" submit_diabatical_incav.sh 
    sed -i "s/RUN=eq/RUN=eq_c_$CONCENTRATION\_freq_$CavFreq/" submit_diabatical_incav.sh 
    sed -i "s/ipi_co2.eq/ipi_co2.eq_c_$CONCENTRATION\_freq_$CavFreq/" submit_diabatical_incav.sh 
    sed -i "s/n4/$NODE/" submit_diabatical_incav.sh 
    
    sed -i "s/co2_eq/polariton_relaxation_biased\/equilibrium_simulation\/eq_c_$CONCENTRATION\_freq_$CavFreq/" submit_diabatical_outcav.sh 
    sed -i "s/RUN=eq/RUN=eq_c_$CONCENTRATION\_freq_$CavFreq/" submit_diabatical_outcav.sh 
    sed -i "s/ipi_co2.eq/ipi_co2.eq_c_$CONCENTRATION\_freq_$CavFreq/" submit_diabatical_outcav.sh 
    sed -i "s/n4/$NODE/" submit_diabatical_outcav.sh 
   
    # 4. change photon frequency 
    sed -i "s/2320.0/$CavFreq/" photon_params.json
    
    # 5. submit jobs to diabatical
    qsub submit_diabatical_incav.sh
    qsub submit_diabatical_outcav.sh

    cd ..

done
