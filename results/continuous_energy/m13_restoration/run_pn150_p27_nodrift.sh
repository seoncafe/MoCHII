#!/bin/bash
# PN150 restart (2026-10-13) from pn150_p27 (2^27 packets) after 36 + 32
# iterations (pn150_p27_it32_rates.h5), with the drift test off
# (par%conv_drift_test = .false.): there the x_HII step sat at 1.6e-3 with a
# net outward drift of ~1e-5 of the ionized volume per iteration (t 3.3-3.9),
# every other measure passing.  Tree binary of 2026-10-13 (final gas solve,
# the option) copied as MoCHII_pn150_nodrift.x; lart2 + lart3, 36 ranks each;
# the iterate written every iteration.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
n=pn150_p27_nodrift
echo "${n} start lart2+lart3 36+36 $(date '+%F %T')" >> driver.log
mpirun -hosts lart2,lart3 -ppn 36 -n 72 ./MoCHII_pn150_nodrift.x ${n}.in > ${n}.log 2>&1
echo "${n} exit=$? $(date '+%F %T')" >> driver.log
