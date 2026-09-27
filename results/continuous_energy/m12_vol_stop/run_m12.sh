#!/bin/bash
# conv_crit = 'vol' for the runs that stopped on the cell test and always hit
# gas_niter: the G1 Stromgren set and the G4 re-refinement pair (in their
# tests/ directories), the dusty Stromgren trio and the two dust-emission
# smokes (tests/d_dusty; the smokes rerun unchanged with the current binary
# for the paper's dust maps), and G2a restarted from m11 here.  Binary:
# ../../../MoCHII_m12run.x (a copy of MoCHII.x of 2026-09-27, kept at the
# repository root so SEDust/data resolves).  One job per node from a shared
# queue: lart1, lart2 (72 ranks), mocafe (88).
R=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00
cd $R/results/continuous_energy/m12_vol_stop
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
printf '%s\n' tests/g1_stromgren:g1_stromgren_l7 tests/g4_refine:g4_uniform7 \
  results/continuous_energy/m12_vol_stop:g2a_thermal tests/g4_refine:g4_refine \
  tests/g1_stromgren:g1_stromgren tests/g1_stromgren:g1_stromgren_ref tests/g1_stromgren:g1_stromgren_front \
  tests/d_dusty:d_dust_full tests/d_dusty:d_dust_l09 tests/d_dusty:d_dust_scat \
  tests/d_dusty:pahlive_smoke tests/d_dusty:dustemis_smoke > queue.txt
worker() {
  host=$1; np=$2
  while true; do
    e=$(flock queue.lock bash -c 'n=$(head -1 queue.txt); [ -n "$n" ] && sed -i 1d queue.txt; echo $n')
    [ -z "$e" ] && break
    d=${e%%:*}; n=${e##*:}
    echo "${n} start ${host} $(date '+%F %T')" >> driver.log
    (cd $R/$d && mpirun -hosts ${host} -n ${np} $R/MoCHII_m12run.x ${n}.in > ${n}.log 2>&1)
    echo "${n} exit=$? ${host} $(date '+%F %T')" >> driver.log
  done
}
worker lart1 72 &
worker lart2 72 &
worker mocafe 88 &
wait
rm -f queue.txt queue.lock
echo "ALL DONE $(date '+%F %T')" >> driver.log
