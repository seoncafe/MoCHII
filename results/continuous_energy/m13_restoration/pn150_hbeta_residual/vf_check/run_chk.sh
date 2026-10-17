#!/bin/bash
# Channel-8 check (2026-10-17): 'none' with the new binary against the
# binary before (bit identity), and 'valence_fractions' for the channel-8
# luminosity; 2 iterations of 2^20 packets from the converged l6c_metals state.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/pn150_hbeta_residual/vf_check
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
sed "s/chk_none.h5/chk_old.h5/" chk_none.in > chk_old.in
mpirun -hosts lart1 -n 36 ../MoCHII_l6.x chk_old.in > chk_old.log 2>&1 &
mpirun -hosts lart2 -n 36 ../MoCHII_vf.x chk_none.in > chk_none.log 2>&1 &
mpirun -hosts lart3 -n 36 ../MoCHII_vf.x chk_vf.in > chk_vf.log 2>&1 &
wait
