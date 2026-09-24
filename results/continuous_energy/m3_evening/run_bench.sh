#!/bin/bash
# The two Lexington production runs of m2_atomic/ rerun with the current
# binary: measures the He + H+ (Kimura) rate that replaced the detailed-balance
# reverse; use_sec_ion stays at its default (off), and the 2^3S update is a
# diagnostic with no feedback.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m3_evening
until grep -q "ALL DONE" driver.log 2>/dev/null; do sleep 20; done
run() { echo "$1 start $(date +%T)" >> driver_bench.log; mpirun -np 32 ./MoCHII_run.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver_bench.log; }
run hii40_evening & run hii20_evening &
wait; echo "ALL DONE $(date +%T)" >> driver_bench.log
