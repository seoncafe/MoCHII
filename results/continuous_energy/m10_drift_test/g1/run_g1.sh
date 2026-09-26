#!/bin/bash
# The three G1 Stromgren runs of the paper figure with the binary of the m9
# set (hydrogenic alpha_B of 2026-09-26), on lart2 (made in
# ../../m9_hydrogenic_alpha/g1/ and moved here unchanged the same day); the inputs are those of
# tests/g1_stromgren with the grid path made absolute.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m10_drift_test/g1
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
one() { echo "$1 start $(hostname) $(date +%T)" >> driver.log; mpirun -np 10 ./MoCHII_g1.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver.log; }
export -f one
printf '%s\n' g1_stromgren g1_stromgren_ref g1_stromgren_front | xargs -P 3 -I{} bash -c 'one {}'
echo "ALL DONE $(date +%T)" >> driver.log
