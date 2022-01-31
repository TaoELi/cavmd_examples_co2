#!/bin/bash

# a line beginning with # is a comment
# a line beginning with #PBS is a PBS directive
# PBS directives must come first
# any directives after the first executable statement are ignored
# search for $ and replace them with proper strings

# job name
#PBS -N debug

# resource requested
#PBS -l nodes=n10:ppn=32
#PBS -l walltime=3000:00:00
#PBS -l mem=48GB

# working directory
# PBSs -d $DIR

# output file
# PBSdas -o $DIR/x.out

# error file
# PBsdsaS -e $DIR/x.err


OMP_NUM_THREADS=32
MKL_NUM_THREADS=32


cd $PBS_O_WORKDIR

module load arpack/3.6.2 cmake/3.12.2 fftw/3.3.8 gcc/10.2 gsl/2.5 hdf5/1.10.5 lapack/3.8.0 mkl/2018.3.222 mpich/3.2.1 openblas/0.2.20  openmpi/4.0.1 openssl python/3.7 scalapack/2.0.2 zlib/1.2.11 utilities/2020
export PYTHONPATH=$PYTHONPATH:~/lib/python3.7/site-packages/
export PATH=~/bin:$PATH
source ~/.bashrc


#CURRENT_PATH=$(pwd)
CURRENT_PATH=$PBS_O_WORKDIR
Global_path=/scratch/taoli/job/polariton_lifetime_vsc/pumping_detuning/
RUN=detuning

PATTERN='<step>80000</step>'

E0=2e-4
Freq_lst=(2100 2150 2200 2250 2300 2325 2350 2400 2450 2500 2550)
LP_lst=(  2087 2129 2167 2202 2232 2244 2255 2272 2284 2294 2299)
UP_lst=(  2372 2379 2387 2401 2419 2431 2444 2476 2512 2552 2596)

for Polariton in LP UP
do
for Amp in 6e-4 6e-3
do
    for idx in ${!Freq_lst[@]}
    do

	Freq=${Freq_lst[$idx]}
	if [[ $Polariton == LP ]]
	then
    	    PolaritonFreq=${LP_lst[$idx]}
	else
	    PolaritonFreq=${UP_lst[$idx]}
	fi
	
	DIR=$Global_path/Freq_$Freq\_Exc_$Polariton\_Amp_$Amp
	
	echo "Working with $DIR"
	echo "Freq = $Freq, PolaritonFreq = $PolaritonFreq"

        cd $PBS_O_WORKDIR
        mkdir -p $DIR
        cp data.lmp photon_params.json ../eq_detuning/Freq_$Freq/init_* input_traj.xml.bak  in.lmp $DIR
        cd $DIR
        echo "move in $DIR"
        
        # modify parameter control of photon_params.json 
        sed -i "s/2320/$Freq/" photon_params.json
        sed -i "s/0.0005/$E0/" photon_params.json
        sed -i "s/2241/$PolaritonFreq/" photon_params.json
        sed -i "s/0.05/$Amp/" photon_params.json
        # Here, we prepare the input files for each trajectory
        for traj in {1..40}
        do
        (
        	echo "Dealing with $traj slice"
            cp input_traj.xml.bak input_traj_$traj.xml
            cp in.lmp in_$traj.lmp
            # change the i-pi job name for each traj
            sed -i "s/mesitylene-pimd.1/$RUN.Freq_$Freq.$Polariton.Amp_$Amp.traj_$traj/" input_traj_$traj.xml
            sed -i "s/mesitylene-pimd.1/$RUN.Freq_$Freq.$Polariton.Amp_$Amp.traj_$traj/" in_$traj.lmp
            # change the input file for each traj 	
            sed -i "s/RESTART/init_$(($traj-1)).checkpoint/" input_traj_$traj.xml
        	sed -i "s/'simu'/'simu_$traj'/g" input_traj_$traj.xml
    
       
        	CHECKPOINT=simu_$traj.checkpoint
        	if grep -q $PATTERN $CHECKPOINT; then
        	    echo "found checkpoint finished"
        	    echo "Skip $traj-th sequential job"
        	else
        	    echo "not found checkpoint finished"
        	    a=$(wc -c < $CHECKPOINT)
        	    if [ ! -f "$CHECKPOINT" ] || [ $a -le 1 ]; then
        	       echo "Performing $traj-th simulation 40 ps"
        	       i-pi input_traj_$traj.xml &> log_$traj &
        	       sleep 30s
        	       lmp < in_$traj.lmp &> info_lmp_$traj.log
        	       sleep 30s
        	    else
        	       echo "Continuing $traj-th simulation 40 ps"
        	       i-pi $CHECKPOINT &> log_$traj &
        	       sleep 30s
        	       lmp < in_$traj.lmp &> info_lmp_$traj.log
        	       sleep 30s
        	    fi
        	fi
        )&
        done
        wait 
    done
done
done
