#!/bin/bash
# PN150 restart (2026-10-13) from the 2^25-packet run pn150_2node after 36
# iterations (pn150_2node_it36_rates.h5), with 2^27 packets; the earlier header:
# PN150 rerun (2026-10-12): tree binary after the node-shared tables and the
# iterate output (copied as MoCHII_pn150.x), lart2 + lart3 with 36 ranks each,
# the iterate written every 2 iterations to pn150_2node_rates_iterate.h5
# (pn150_p27.in: every iteration, to pn150_p27_rates_iterate.h5).
# Input as pn150_converged.in (2^25 packets) except the output name and
# rates_write_interval.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
n=pn150_p27
echo "${n} start lart2+lart3 36+36 $(date '+%F %T')" >> driver.log
mpirun -hosts lart2,lart3 -ppn 36 -n 72 ./MoCHII_pn150.x ${n}.in > ${n}.log 2>&1
echo "${n} exit=$? $(date '+%F %T')" >> driver.log
