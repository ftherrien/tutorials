#!/bin/bash

mkdir -p POSCARs

for fname in $1/*.mol
do
    echo $fname
    if test -f POSCARs/POSCAR_$(basename $fname .mol)
    then
    echo "Already exists."
    else
    obabel -imol $fname -ovasp -O POSCARs/POSCAR_$(basename $fname .mol)
    python center.py POSCARs/POSCAR_$(basename $fname .mol)
    fi
done
