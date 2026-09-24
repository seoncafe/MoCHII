#!/bin/bash
# Convergence check of the HII40 He front: the production input run for a
# fixed 60 iterations (tolerances out of reach), because the stopping test
# looks at x(H+) and Te only.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m4_nitrogen
until grep -q "ALL DONE" driver_rest.log 2>/dev/null; do sleep 30; done
echo "start $(date +%T)" >> driver_niter60.log
mpirun -np 64 ./MoCHII_run.x hii40_niter60_nitrogen.in > hii40_niter60_nitrogen.log 2>&1
echo "exit=$? $(date +%T)" >> driver_niter60.log; echo "ALL DONE $(date +%T)" >> driver_niter60.log
