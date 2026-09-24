#!/bin/bash
# Final production runs for the paper (2026-09-24 atomic-data update), on lart3.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m2_atomic
run() { echo "$1 start $(date +%T)" >> driver.log; mpirun -np 32 ./MoCHII_run.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver.log; }
run hii40_atomic &
run hii20_atomic &
wait
echo "ALL DONE $(date +%T)" >> driver.log
