#!/bin/bash
# The m8 paper set restarted from its converged states with the hydrogenic
# recombination of 2026-09-26: alpha_B(H II) and alpha_B(He III) from the
# nl-resolved Storey & Hummer cascade (hydrogenic_recombination_mod) instead
# of the Badnell RR fit minus the Milne alpha_1 (1.13% lower at 1e4 K for
# H II), the 2s share and the He II Lyman exits from the same cascade, and
# the nebular-continuum table gaps filled.  g2a restarts from g2a_current.
# Three 24-rank runs at a time.  g2a_thermal_m9 and the g1/ runs were moved
# to ../m10_drift_test/ (as g2a_thermal.* and g1/) on 2026-09-26, where the
# paper reads them.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m9_hydrogenic_alpha
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
one() { echo "$1 start $(date +%T)" >> driver.log; mpirun -np 24 ./MoCHII_m9.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver.log; }
export -f one
printf '%s\n' hii40_converged hii20_converged g2a_thermal_m9 \
  hii40_ots_converged hii20_ots_converged pdr_L5_converged \
  hei_high_density_converged hii40_lowden_converged hii20_lowden_converged \
  hii40_nohecx_converged hii20_nohecx_converged \
  | xargs -P 3 -I{} bash -c 'one {}'
echo "ALL DONE $(date +%T)" >> driver.log
