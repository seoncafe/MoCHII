#!/bin/bash
# The m10 paper set restarted from its converged states with the binary of
# 2026-09-27: He+ + H non-radiative charge exchange from Loreau et al. (2014)
# by reciprocity (was KF96), the radiative channel from Stancil, Lepp &
# Dalgarno (1998) (was GJ07 R26), the hydrogen balance solved for
# u = ln(n_HII/n_HI), x_HII and x_HeIII stored as their own state, and the
# n_e guard lowered to 1e-30 n_H.  Inputs otherwise unchanged.  g2a_thermal
# restarts from the m10 g2a run (cell test, 40 iterations).
# One job per node from a shared queue: lart1 and lart2 (72 ranks).
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m11_hecx_solver
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
printf '%s\n' g2a_thermal hii40_converged hii20_converged hii40_lowden_converged hii40_nohecx_converged \
  hii20_lowden_converged hii20_nohecx_converged hii40_ots_converged hii20_ots_converged \
  pdr_L5_converged hei_high_density_converged > queue.txt
worker() {
  host=$1; np=$2
  while true; do
    n=$(flock queue.lock bash -c 'n=$(head -1 queue.txt); [ -n "$n" ] && sed -i 1d queue.txt; echo $n')
    [ -z "$n" ] && break
    echo "${n} start ${host} $(date '+%F %T')" >> driver.log
    mpirun -hosts ${host} -n ${np} ./MoCHII_m11.x ${n}.in > ${n}.log 2>&1
    echo "${n} exit=$? ${host} $(date '+%F %T')" >> driver.log
  done
}
worker lart1 72 &
worker lart2 72 &
wait
rm -f queue.txt queue.lock
echo "ALL DONE $(date '+%F %T')" >> driver.log
