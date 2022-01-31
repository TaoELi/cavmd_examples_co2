#!/bin/bash

CavFreqArray=(2320.0)
CW_FreqArray=(2241.0)

C13OFreq=2262.0

SUBMIT_SCRIPT_INCAV=submit_diabatical_incav.sh
SUBMIT_SCRIPT_OUTCAV=submit_diabatical_outcav.sh

for k in ${!CavFreqArray[@]}
do
    CavFreq=${CavFreqArray[$k]}
    CW_Freq=${CW_FreqArray[$k]}

    for CONCENTRATION in 1.0 0.954
    do
    
        # 1. create a dictionary
        DIR=pump_c_$CONCENTRATION\_freq_$CavFreq
        mkdir $DIR
        cd $DIR
        
        # 2. create nonequilibrium scripts for simulation
        cp ../template_pump/* .
    
        # 3. Modify input information
        sed -i "s/Concentration=0.0/Concentration=$CONCENTRATION/" $SUBMIT_SCRIPT_INCAV
        sed -i "s/CavFreq=2200.0/CavFreq=$CavFreq/" $SUBMIT_SCRIPT_INCAV
        sed -i "s/CW_Freq=2241.0/CW_Freq=$CW_Freq/" $SUBMIT_SCRIPT_INCAV
        
        sed -i "s/Concentration=0.0/Concentration=$CONCENTRATION/" $SUBMIT_SCRIPT_OUTCAV
        sed -i "s/CavFreq=2200.0/CavFreq=$CavFreq/" $SUBMIT_SCRIPT_OUTCAV
        sed -i "s/CW_Freq=2241.0/CW_Freq=$C13OFreq/" $SUBMIT_SCRIPT_OUTCAV
        
        # 4. submit jobs
        
        qsub $SUBMIT_SCRIPT_INCAV
        qsub $SUBMIT_SCRIPT_OUTCAV
        cd ..
    
    done
done
