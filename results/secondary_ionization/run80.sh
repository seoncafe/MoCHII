#!/bin/bash
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/secondary_ionization
until grep -q "ALL DONE" ../hei_metastable/driver.log 2>/dev/null; do sleep 20; done
run() { echo "$1 start $(date +%T)" >> driver80.log; mpirun -np 24 ./$2 sec_$1.in > sec_$1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver80.log; }
run on_new MoCHII_new.x & run on_old MoCHII_old.x & run off_new MoCHII_new.x &
wait; echo "ALL DONE $(date +%T)" >> driver80.log
