#!/bin/bash
# G2a with the transported recombination continuum (lart2), against the
# on-the-spot m9 run: where does the He III cascade heat land when its
# photons are followed instead of absorbed in place?  Restarted from the
# g2a_current state (on the spot, 2026-09-26 before the hydrogenic alpha_B).
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/g2a_transported
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
echo "start $(hostname) $(date +%T)" > driver.log
mpirun -np 40 ./MoCHII_g2at.x g2a_thermal_transported.in > g2a_thermal_transported.log 2>&1
echo "exit=$? $(date +%T)" >> driver.log
