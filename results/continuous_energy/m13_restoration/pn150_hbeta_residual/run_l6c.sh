#!/bin/bash
# PN150 level-6 controls continued to convergence (2026-10-17): the 60-iteration
# l6_metals / l6_nometals states did not converge, so they are restarted with
# the drift test on until the stopping test passes (cap 200).
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/pn150_hbeta_residual
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
mpirun -hosts lart1 -n 36 ./MoCHII_l6.x l6c_metals.in > l6c_metals.log 2>&1 &
mpirun -hosts lart2 -n 36 ./MoCHII_l6.x l6c_nometals.in > l6c_nometals.log 2>&1 &
wait
echo done >> run_l6c.done
