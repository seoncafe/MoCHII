#!/bin/bash
# PN150 L(Hbeta) residual: the level-6 grid with metals (lart2) and without
# (lart3), 2^24 packets, 60 iterations from a fresh start, so the H
# recombinations per stellar photon can be set against Cloudy's 1.118 / 1.103.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/pn150_hbeta_residual
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
mpirun -hosts lart2 -n 36 ./MoCHII_l6.x l6_metals.in > l6_metals.log 2>&1 &
mpirun -hosts lart3 -n 36 ./MoCHII_l6.x l6_nometals.in > l6_nometals.log 2>&1 &
wait
echo done >> run_l6.done
