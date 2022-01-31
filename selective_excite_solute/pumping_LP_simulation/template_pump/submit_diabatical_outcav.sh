#!/bin/bash

# a line beginning with # is a comment
# a line beginning with #PBS is a PBS directive
# PBS directives must come first
# any directives after the first executable statement are ignored
# search for $ and replace them with proper strings

# job name
#PBS -N debug2

# resource requested
#PBS -l nodes=n10:ppn=1
#PBS -l walltime=3000:00:00
#PBS -l mem=4GB

# working directory
# PBSs -d $DIR

# output file
# PBSdas -o $DIR/x.out

# error file
# PBsdsaS -e $DIR/x.err


#OMP_NUM_THREADS=32
#MKL_NUM_THREADS=32


cd $PBS_O_WORKDIR

module load arpack/3.6.2 cmake/3.12.2 fftw/3.3.8 gcc/10.2 gsl/2.5 hdf5/1.10.5 lapack/3.8.0 mkl/2018.3.222 mpich/3.2.1 openblas/0.2.20  openmpi/4.0.1 openssl python/3.7 scalapack/2.0.2 zlib/1.2.11 utilities/2020
export PYTHONPATH=$PYTHONPATH:~/lib/python3.7/site-packages/
export PATH=~/bin:$PATH
source ~/.bashrc


Concentration=0.0
CavFreq=2200.0
E0=0e-4
CW_Freq=2241.0

RUN=pump_c_$Concentration\_freq_$CavFreq\_$E0\_$CW_Freq

CURRENT_PATH=$PBS_O_WORKDIR
Global_path=/scratch/taoli/job/polariton_relaxation_biased/pumping_LP_simulation/pump_c_$Concentration\_freq_$CavFreq/

cd $PBS_O_WORKDIR

for Amp in 6e-3
do
    DIR=$Global_path/Exc_$CW_Freq\_Amp_$Amp\_E0_$E0
    
    mkdir -p $DIR
    
    cp data.lmp input_traj.xml.bak create_photon_params.py in.lmp ../../equilibrium_simulation/eq_c_$Concentration\_freq_$CavFreq/E0_$E0/init*checkpoint $DIR
    cd $DIR
    echo "move in $DIR"
    
    # create a proper photon_params.json file 
    sed -i "s/0.0005/$E0/" create_photon_params.py
    sed -i "s/0.05,/$Amp,/" create_photon_params.py
    sed -i "s/2241.0/$CW_Freq/" create_photon_params.py
    sed -i "s/2200.0/$CavFreq/" create_photon_params.py
    python create_photon_params.py
    mv test.json photon_params.json
    
    PATTERN='<step>40000</step>'
    
    # Here, we prepare the input files for each trajectory
    for traj in {1..40}
    do
    	echo "Dealing with $traj slice"
        cp input_traj.xml.bak input_traj_$traj.xml
        cp in.lmp in_$traj.lmp
        # change the i-pi job name for each traj
        sed -i "s/mesitylene-pimd.1/co2.$RUN.traj_$traj/" input_traj_$traj.xml
        sed -i "s/mesitylene-pimd.1/co2.$RUN.traj_$traj/" in_$traj.lmp
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
    done
    wait 
done
