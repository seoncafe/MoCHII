#!/bin/bash
# The level-7 G1 Stromgren run and the two G4 re-refinement runs with the
# current binary (hydrogenic alpha_B of 2026-09-26), for the paper's G1/G4
# numbers; the inputs are tests/g1_stromgren/g1_stromgren_l7.in and
# tests/g4_refine/g4_{refine,uniform7}.in with the grid and data paths made
# absolute.  One job per node: lart1, lart2 (72 ranks), mocafe (88).
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m10_drift_test/g1_g4
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
one() { echo "$1 start $2 $(date '+%F %T')" >> driver.log; mpirun -hosts $2 -n $3 ./MoCHII_m10.x $1.in > $1.log 2>&1; echo "$1 exit=$? $2 $(date '+%F %T')" >> driver.log; }
one g1_stromgren_l7 lart1 72 &
one g4_uniform7 lart2 72 &
one g4_refine mocafe 88 &
wait
echo "ALL DONE $(date '+%F %T')" >> driver.log
