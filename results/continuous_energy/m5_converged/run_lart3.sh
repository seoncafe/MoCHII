#!/bin/bash
# The paper set rerun with the drift stopping test (conv_crit = 'vol'), the
# vacuum stage-fraction fix of the energy ledger, and the HII40 gas extended
# to 5.2 pc (tests/g2_hii/hii40_L6_rout5p2.fits; the MOCASSIN benchmark's
# Rout = 1.46e19 cm sits on the converged H front).  g2a_thermal uses the
# cell test and is not rerun (m4_nitrogen).  Three 24-rank runs at a time.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m5_converged
one() { echo "$1 start $(date +%T)" >> driver.log; mpirun -np 24 ./MoCHII_run.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver.log; }
export -f one
printf '%s\n' hii40_converged hii40_lowden_converged hii40_ots_converged \
  hii40_nohecx_converged hii20_converged hii20_lowden_converged \
  hii20_ots_converged hii20_nohecx_converged pdr_L5_converged hei_high_density_converged \
  | xargs -P 3 -I{} bash -c 'one {}'
echo "ALL DONE $(date +%T)" >> driver.log
