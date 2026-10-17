#!/bin/bash
# Absorption ledger (par%absorption_ledger) on the converged PN150 L6 states:
# one transport pass with the restored state held (gas_niter = 0); and the
# ledger-off binary against the binary before on the 2-iteration check.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/pn150_hbeta_residual
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
mpirun -hosts lart1 -n 36 ./MoCHII_ledger.x ledger_l6c_metals.in > ledger_l6c_metals.log 2>&1 &
mpirun -hosts lart2 -n 36 ./MoCHII_ledger.x ledger_l6vf_metals.in > ledger_l6vf_metals.log 2>&1 &
mpirun -hosts lart3 -n 36 ./MoCHII_ledger.x ledger_l6c_nometals.in > ledger_l6c_nometals.log 2>&1 &
wait
cd vf_check
mpirun -hosts lart1 -n 36 ../MoCHII_ledger.x chk_ledgeroff.in > chk_ledgeroff.log 2>&1
