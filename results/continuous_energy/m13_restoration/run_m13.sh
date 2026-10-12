#!/bin/bash
# Lexington HII40 / HII20 under the restored atomic physics of 2026-10-11
# (md/ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md, Phase 5 step 2): stages
# I-VII of every metal (77 stages / 66 transitions), full CHIANTI 11 atoms,
# phfit2 subshell absorption, metal free-bound and free-free diffuse
# channels, Ne III from Nahar (2019), HEAD-default switches on.  Restarted
# from the m11 converged states; the added metal transitions start at zero.
# Inputs are the m11 inputs with the restart path changed and the excitation
# model stated.  One job per host: lart1 (72 ranks), mocafe (88 ranks).
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
run() {
  n=$1; host=$2; np=$3
  echo "${n} start ${host} $(date '+%F %T')" >> driver.log
  mpirun -hosts ${host} -n ${np} ./MoCHII_m13.x ${n}.in > ${n}.log 2>&1
  echo "${n} exit=$? ${host} $(date '+%F %T')" >> driver.log
}
run hii40_converged lart1 72 &
run hii20_converged mocafe 88 &
wait
echo "ALL DONE $(date '+%F %T')" >> driver.log
