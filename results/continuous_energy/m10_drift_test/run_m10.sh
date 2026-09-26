#!/bin/bash
# The m9 paper set rerun to convergence under the stopping test of
# 2026-09-26 (fixed_point_distance_mod: rho fitted over 8 iterations, and the
# Monte Carlo floor accepted only when the signed change of the ionized
# volume has no resolvable mean).  Under the test before it, 8 of these 10
# runs stopped after 4-5 iterations with their fronts still moving.  Each
# input is its m9 counterpart with the restart taken from the m9 state
# (hii40 from hii40_long, 16 iterations past hii40_converged); physics and
# tolerances are unchanged.  g2a is not rerun: it uses conv_crit = 'cell'.
# One job per node at a time, pulled from a shared queue: lart1 and lart2
# (72 ranks each) and mocafe (88).
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m10_drift_test
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
printf '%s\n' hii40 hii20 hii40_lowden hii20_lowden hii40_nohecx hii20_nohecx \
  hii40_ots hii20_ots pdr_L5 hei_high_density > queue.txt
worker() {
  host=$1; np=$2
  while true; do
    n=$(flock queue.lock bash -c 'n=$(head -1 queue.txt); [ -n "$n" ] && sed -i 1d queue.txt; echo $n')
    [ -z "$n" ] && break
    echo "${n} start ${host} np=${np} $(date '+%F %T')" >> driver.log
    mpirun -hosts ${host} -n ${np} ./MoCHII_m10.x ${n}_converged.in > ${n}_converged.log 2>&1
    echo "${n} exit=$? ${host} $(date '+%F %T')" >> driver.log
  done
}
worker lart1 72 &
worker lart2 72 &
worker mocafe 88 &
wait
echo "ALL DONE $(date '+%F %T')" >> driver.log
