#!/bin/bash
# Final solve check (2026-10-12): the Phase-5 attribution input it2/d
# (1e5 K blackbody, 2 iterations from the initial guess) on the binary with
# the final gas solve after the final pass; lart3, 16 ranks as it2/d.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/final_solve_check
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
mpirun -hosts lart3 -n 16 ./MoCHII_final_solve.x bb1e5_it2.in > bb1e5_it2.log 2>&1
echo "exit $?" >> bb1e5_it2.log
