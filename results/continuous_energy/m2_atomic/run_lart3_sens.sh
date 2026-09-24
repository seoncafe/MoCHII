#!/bin/bash
# Sensitivity runs quoted in the paper: zero-density cooling limit, case-B on the spot.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m2_atomic
until grep -q "ALL DONE" driver_nohecx.log 2>/dev/null; do sleep 20; done
run() { echo "$1 start $(date +%T)" >> driver_sens.log; mpirun -np 18 ./MoCHII_run.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver_sens.log; }
run hii40_atomic_lowden & run hii20_atomic_lowden & run hii40_atomic_ots & run hii20_atomic_ots &
wait
echo "ALL DONE $(date +%T)" >> driver_sens.log
