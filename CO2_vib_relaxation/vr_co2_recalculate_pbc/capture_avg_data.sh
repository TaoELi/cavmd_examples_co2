#!/bin/bash

for Nsub in 216 302 324 432 648 864 1080 1296 1728 2160 2592 3456 4320 5184 6048 6912 8000 12960
do
    echo "Let Nsub = $Nsub"
    DIR=vr_$Nsub

    E0_incav=$(echo "" | awk "END {print 2e-4 * (216.0 / $Nsub)^(0.5)}")
    E0_incav=$(printf "%.4e\n" ${E0_incav})

    Nhot=$(echo "" | awk "END {print $Nsub * (10.0/216.0)}")
    Nhot=$(printf "%1.f" $Nhot)
    echo "Nhot = $Nhot"
    for Nhot in $Nhot
    do
        python obtain_avg.py $DIR/Exc_$Nhot\_E0_0e-4
        python obtain_avg.py $DIR/Exc_$Nhot\_E0_$E0_incav
        python obtain_avg.py $DIR/Exc_10\_E0_0e-4
        python obtain_avg.py $DIR/Exc_10\_E0_$E0_incav
        python obtain_ph_spectrum.py $DIR/Exc_10\_E0_$E0_incav
        python obtain_ph_spectrum.py $DIR/Exc_$Nhot\_E0_$E0_incav
    done
done

Nsub=2160
DIR=vr_$Nsub
E0_incav=$(echo "" | awk "END {print 2e-4 * (216.0 / $Nsub)^(0.5)}")
E0_incav=$(printf "%.4e\n" ${E0_incav})
for Nhot in 20 30 40 50 60 70 80 90
do
        python obtain_avg.py $DIR/Exc_$Nhot\_E0_0e-4
        python obtain_avg.py $DIR/Exc_$Nhot\_E0_$E0_incav
        python obtain_ph_spectrum.py $DIR/Exc_$Nhot\_E0_$E0_incav
done

