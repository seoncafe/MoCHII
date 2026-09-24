#!/bin/bash
# The two Lexington production runs of m2_atomic/ rerun after the evening
# charge-exchange corrections of 2026-09-24: N0/N+ + H from Kingdon & Ferland
# (1996) in place of the Huang et al. (2023) polynomials, the form-7
# Boltzmann factor at the actual T, and He + H+ back to detailed balance.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m4_nitrogen
run() { echo "$1 start $(date +%T)" >> driver.log; mpirun -np 32 ./MoCHII_run.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver.log; }
run hii40_nitrogen & run hii20_nitrogen &
wait; echo "ALL DONE $(date +%T)" >> driver.log
