#!/bin/bash

for Nsub in 216 302 324 432 648 864 1080 1296 1728 2160 2592 3456 4320 5184 6048 6912 8000 12960
do
    echo "Let Nsub = $Nsub"
    DIR=eq_$Nsub
    cd $DIR
    echo "working in" 
    pwd
    
    qsub submit_diabatical_outcav.sh
    
    cd .. 
done
