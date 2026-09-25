#!/bin/bash
# Convergence of the gas iteration: the HII40 production input run for a
# fixed 80 iterations (tolerances out of reach, require_convergence off),
# with the fixed-point distance of x_HII and He0-He+ printed every iteration.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/convergence
echo "hii40_niter80 start $(date +%T)" >> driver.log
mpirun -np 64 ./MoCHII_run.x hii40_niter80.in > hii40_niter80.log 2>&1
echo "hii40_niter80 exit=$? $(date +%T)" >> driver.log; echo "ALL DONE" >> driver.log
