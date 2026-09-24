#!/bin/bash
# Figure runs for the paper, after the two benchmark runs, on lart3.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m2_atomic
until grep -q "ALL DONE" driver.log 2>/dev/null; do sleep 20; done
run() { echo "$1 start $(date +%T)" >> driver_fig.log; mpirun -np 24 ./MoCHII_run.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver_fig.log; }
run pdr_L5_atomic &
run hei_high_density_atomic &
run g2a_thermal_atomic &
wait
echo "ALL DONE $(date +%T)" >> driver_fig.log
