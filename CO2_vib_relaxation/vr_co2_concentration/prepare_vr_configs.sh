#!/bin/bash

for CONCENTRATION in 0.2 0.4 0.6 0.8
do

    # create a dictionary
    DIR=c_$CONCENTRATION
    mkdir $DIR
    cd $DIR
    
    # 1. copy input files from equilibrium simulations
    mkdir equilibrium_configs
    mkdir modified_configs_exc_10
    cp ../eq_c_$CONCENTRATION/E0_0e-4/* equilibrium_configs/
    
    # 2. modify the initial conditions
    cp ../change_initial_condition_exc_10.py .
    python change_initial_condition_exc_10.py equilibrium_configs/init_*

    # 3. copy nonequilibrium scripts for simulation
    cp ../../try_VER/data.lmp .
    cp ../../try_VER/in.lmp .
    cp ../../try_VER/photon_params.json .
    cp ../../try_VER/input_traj.xml.bak .
    cp ../../try_VER/submit_diabatical.sh .

    # 4. change submission job details
    #sed -i "s/try_VER\/N_10/vr_co2_concentration\/c_$CONCENTRATION/" submit_diabatical.sh
    sed -i "s/vib_relax/vr_c_$CONCENTRATION/g" submit_diabatical.sh
    sed -i "s/modified_configs/modified_configs_exc_10/g" submit_diabatical.sh
    cd ..

done
