#!/bin/bash
# He I 2^3S cases, current binary (here) and the 12:40 binary (old/), on lart3.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/hei_metastable
until grep -q "ALL DONE" ../secondary_ionization/driver.log 2>/dev/null; do sleep 20; done
run() { echo "$1/$2 start $(date +%T)" >> $(pwd)/driver.log; (cd $1 && mpirun -np 12 $3 $2.in > $2.log 2>&1); echo "$1/$2 exit=$? $(date +%T)" >> $(pwd)/driver.log; }
for b in hei_lowd hei_highd hei_fuv; do run . $b ./MoCHII_new.x & run old $b ./MoCHII_old.x & done
wait; echo "ALL DONE $(date +%T)" >> driver.log
