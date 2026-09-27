#!/bin/bash
# Rerun of the G1 set, the G4 pair and the dusty trio with gas_tol = 5e-3:
# under the drift test (conv_crit = 'vol') with gas_tol = 2e-3 they reached
# the Monte Carlo floor (no resolvable drift) but the He0-He+ floor there,
# 2.5-3.4e-3 at 4e7 packets, lies above 2e-3, so the test could never pass
# (run_m12.sh logs, 2026-09-27).  Binary ../../../MoCHII_m12run.x; lart1 and
# lart2 (72 ranks) from a shared queue.
R=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00
cd $R/results/continuous_energy/m12_vol_stop
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
printf '%s\n' tests/g1_stromgren:g1_stromgren_l7 tests/g4_refine:g4_uniform7 tests/g4_refine:g4_refine \
  tests/g1_stromgren:g1_stromgren tests/g1_stromgren:g1_stromgren_ref tests/g1_stromgren:g1_stromgren_front \
  tests/d_dusty:d_dust_full tests/d_dusty:d_dust_l09 tests/d_dusty:d_dust_scat > queue_b.txt
worker() {
  host=$1; np=$2
  while true; do
    e=$(flock queue_b.lock bash -c 'n=$(head -1 queue_b.txt); [ -n "$n" ] && sed -i 1d queue_b.txt; echo $n')
    [ -z "$e" ] && break
    d=${e%%:*}; n=${e##*:}
    echo "${n} start ${host} $(date '+%F %T')" >> driver_b.log
    (cd $R/$d && mpirun -hosts ${host} -n ${np} $R/MoCHII_m12run.x ${n}.in > ${n}.log 2>&1)
    echo "${n} exit=$? ${host} $(date '+%F %T')" >> driver_b.log
  done
}
worker lart1 72 &
worker lart2 72 &
wait
rm -f queue_b.txt queue_b.lock
echo "ALL DONE $(date '+%F %T')" >> driver_b.log
