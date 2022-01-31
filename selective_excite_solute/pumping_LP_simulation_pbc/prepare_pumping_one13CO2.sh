#!/bin/bash

for Nsub in 432 864 1728 3456
do

    # 1. create a dictionary
    DIR=pump_N_$Nsub\_one13CO2
    mkdir $DIR
    cd $DIR
    
    # 2. create nonequilibrium scripts for simulation
    cp ../template_pump/input_traj.xml.bak .
    cp ../template_pump/usr_defined_function.py .
    cp ../template_pump/create_photon_params.py .
    cp ../template_pump/submit_diabatical_incav.sh .
    cp ../../equilibrium_simulation_pbc/eq_N_$Nsub\_one13CO2/data.lmp .
    cp ../../equilibrium_simulation_pbc/eq_N_$Nsub\_one13CO2/in.lmp .
    
    # 2. change parameters for create_photon_params.sh create photon_params.json
    sed -i "s/216/$Nsub/" create_photon_params.py
    E0=$(echo "" | awk "END {print 4e-4 * (216.0 / $Nsub)^(0.5)}")
    echo "coupling is $E0"
    sed -i "s/4e-4/$E0/" create_photon_params.py
    python create_photon_params.py
    mv test.json photon_params.json
   
    # 3. change usr_defined_function
    sed -i "s/Nsub=216/Nsub=$Nsub/" usr_defined_function.py 
    #sed -i "s/Nimpurity=10/Nimpurity=10/" usr_defined_function.py 
    
    # 4. change the file of submit jobs
    sed -i "s/Nsub=216/Nsub=$Nsub/" submit_diabatical_incav.sh 
    #sed -i "s/one13CO2/one13CO2/" submit_diabatical_incav.sh 
    
    # 5. change input file
    Ph1=$(($Nsub*3))
    Ph2=$(($Nsub*3+1))
    sed -i "s/648/$Ph1/" input_traj.xml.bak 
    sed -i "s/649/$Ph2/" input_traj.xml.bak 
    
    # 6. submit jobs
    qsub submit_diabatical_incav.sh
    
    cd ..
   
done
