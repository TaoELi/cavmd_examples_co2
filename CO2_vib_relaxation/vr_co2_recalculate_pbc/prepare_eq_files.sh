#!/bin/bash

for Nsub in 216 302 324 432 648 864 1080 1296 1728 2160 2592 3456 4320 5184 6048 6912 8000 12960
do
    echo "Let Nsub = $Nsub"
    DIR=eq_$Nsub
    cp -r eq_template $DIR
    cd $DIR
    echo "working in" 
    pwd
    sed -i "s/NMOLECULE=216/NMOLECULE=$Nsub/" create_inputs.sh
    
    ./create_inputs.sh

    cd .. 
done
