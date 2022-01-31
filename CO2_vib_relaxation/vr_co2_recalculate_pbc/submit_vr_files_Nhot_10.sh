#!/bin/bash

for Nsub in 216 302 324 432 648 864 1080 1296 1728 2160 2592 3456 4320 5184 6048 6912 8000 12960
do
    echo "Let Nsub = $Nsub"
    DIR=vr_$Nsub
    #mkdir $DIR
    #cp -r vr_template/* $DIR
    cd $DIR
    echo "working in" 
    pwd
   
    # prepare photon_params.json file
    #cp ../eq_$Nsub/photon_params.json .
    #cp ../eq_$Nsub/in.lmp .
    #cp ../eq_$Nsub/data.lmp .

    E0_incav=$(echo "" | awk "END {print 2e-4 * (216.0 / $Nsub)^(0.5)}") 
    E0_incav=$(printf "%.4e\n" ${E0_incav})
    
    # change the output of photons 
    #Ph1=$(($Nsub*3))
    #Ph2=$(($Nsub*3+1))
    #sed -i "s/648/$Ph1/" input_traj.xml.bak
    #sed -i "s/649/$Ph2/" input_traj.xml.bak
    
    for Nhot in 10
    do
	# First, create initial configurations under different Nhot
	#python  change_initial_condition.py ../eq_$Nsub/eq_$Nsub/E0_0e-4/ $Nhot
	
        # second, prepare submit jobs inside cavity
        #cp submit_diabatical_exc_Nhot_incav.sh submit_diabatical_exc_$Nhot\_incav.sh	
	#sed -i "s/Nhot=10/Nhot=$Nhot/" submit_diabatical_exc_$Nhot\_incav.sh
	#sed -i "s/216/$Nsub/" submit_diabatical_exc_$Nhot\_incav.sh
	# also change the E0 inside the cavity to make sure the same Rabi splitting
	#sed -i "s/E0_incav=2e-4/E0_incav=$E0_incav/" submit_diabatical_exc_$Nhot\_incav.sh
        
        qsub submit_diabatical_exc_$Nhot\_incav.sh

        # second, prepare submit jobs outside cavity
        #cp submit_diabatical_exc_Nhot_outcav.sh submit_diabatical_exc_$Nhot\_outcav.sh	
	#sed -i "s/Nhot=10/Nhot=$Nhot/" submit_diabatical_exc_$Nhot\_outcav.sh
	#sed -i "s/216/$Nsub/" submit_diabatical_exc_$Nhot\_outcav.sh
	# also change the E0 inside the cavity to make sure the same Rabi splitting
	#sed -i "s/E0_incav=2e-4/E0_incav=$E0_incav/" submit_diabatical_exc_$Nhot\_outcav.sh
	
	qsub submit_diabatical_exc_$Nhot\_outcav.sh
    done
    cd .. 
done
