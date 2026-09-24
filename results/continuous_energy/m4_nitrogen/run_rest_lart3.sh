#!/bin/bash
# The rest of the paper's m2_atomic set (figure and sensitivity runs) rerun on
# the binary of m4_nitrogen, three 24-rank runs at a time on lart3.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m4_nitrogen
one() { echo "$1 start $(date +%T)" >> driver_rest.log; mpirun -np 24 ./MoCHII_run.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver_rest.log; }
export -f one
printf '%s\n' hii40_lowden_nitrogen hii40_ots_nitrogen hii40_nohecx_nitrogen \
  hii20_lowden_nitrogen hii20_ots_nitrogen hii20_nohecx_nitrogen \
  g2a_thermal_nitrogen pdr_L5_nitrogen hei_high_density_nitrogen \
  | xargs -P 3 -I{} bash -c 'one {}'
echo "ALL DONE $(date +%T)" >> driver_rest.log
