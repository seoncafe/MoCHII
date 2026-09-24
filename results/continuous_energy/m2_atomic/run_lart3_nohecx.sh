#!/bin/bash
# Same production inputs with the He-H charge exchange off, to isolate it.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m2_atomic
until grep -q "ALL DONE" driver_fig.log 2>/dev/null; do sleep 20; done
run() { echo "$1 start $(date +%T)" >> driver_nohecx.log; mpirun -np 32 ./MoCHII_run.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver_nohecx.log; }
run hii40_atomic_nohecx &
run hii20_atomic_nohecx &
wait
echo "ALL DONE $(date +%T)" >> driver_nohecx.log
