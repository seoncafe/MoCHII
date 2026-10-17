#!/bin/bash
# Absorption ledger by leaf (2026-10-17): one pass on Cloudy's structure and
# on MoCHII's own converged l6vf_metals state; and the ledger-off check.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/pn150_hbeta_residual
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
mpirun -hosts lart1 -n 36 ./MoCHII_ledger.x fixed_cloudy.in > fixed_cloudy.log 2>&1 &
mpirun -hosts lart2 -n 36 ./MoCHII_ledger.x ledger_l6vf_metals.in > ledger_l6vf_metals.log 2>&1 &
(cd vf_check && mpirun -hosts lart3 -n 36 ../MoCHII_ledger.x chk_ledgeroff.in > chk_ledgeroff.log 2>&1) &
wait
