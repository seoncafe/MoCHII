#!/bin/bash
# Collisional H I Balmer component (par%h_coll_effects) on the saved L6 states:
# one transport pass, the restored gas state held (gas_niter = 0).
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/pn150_hbeta_residual
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
for c in nometals metals; do
  mpirun -hosts mocafe -n 44 ./MoCHII_l6.x l6hc_${c}.in > l6hc_${c}.log 2>&1
done
