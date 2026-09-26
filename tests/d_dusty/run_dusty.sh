#!/bin/bash
# The three dusty Stromgren runs of check_dusty.py with the current binary,
# one per node: lart1, lart2 (72 ranks), mocafe (88).  Run from this directory.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/d_dusty
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
one() { echo "$1 start $2 $(date '+%F %T')" >> run_dusty.log; mpirun -hosts $2 -n $3 ../../MoCHII.x $1.in > $1.log 2>&1; echo "$1 exit=$? $2 $(date '+%F %T')" >> run_dusty.log; }
one d_dust_full lart1 72 &
one d_dust_l09 lart2 72 &
one d_dust_scat mocafe 88 &
wait
echo "ALL DONE $(date '+%F %T')" >> run_dusty.log
