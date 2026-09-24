#!/bin/bash
# The paper's He I 10833 figure run (m2_atomic/hei_high_density_atomic.in)
# rerun with the 2^3S update of the evening of 2026-09-24 (exact collision
# constant, electron-impact ionization of 2^3S, He + H+ from Kimura et al.).
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m3_evening
until grep -q "ALL DONE" ../../secondary_ionization/driver80.log 2>/dev/null; do sleep 20; done
echo "start $(date +%T)" >> driver.log
mpirun -np 24 ./MoCHII_run.x hei_high_density_evening.in > hei_high_density_evening.log 2>&1
echo "exit=$? $(date +%T)" >> driver.log; echo "ALL DONE $(date +%T)" >> driver.log
