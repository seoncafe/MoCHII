#!/bin/bash
# Lexington/Meudon PN150 under the restored atomic physics of 2026-10-11
# (md/ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md, Phase 5 step 2): level-7
# grid tests/g2_hii/pn150_L7.fits, no restart, lart3 with 40 ranks.  At 56
# ranks the kernel OOM killer stopped the run in the diffuse transport of
# iteration 2 (4.8 GB resident per rank on 2097152 leaves; lart3 has
# 251 GB); that log is kept as pn150_converged_oom56.log.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
n=pn150_converged
echo "${n} start lart3 $(date '+%F %T')" >> driver.log
mpirun -hosts lart3 -n 40 ./MoCHII_m13.x ${n}.in > ${n}.log 2>&1
echo "${n} exit=$? lart3 $(date '+%F %T')" >> driver.log
