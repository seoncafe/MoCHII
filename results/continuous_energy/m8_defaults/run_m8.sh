#!/bin/bash
# The paper set of m5_converged (and g2a_thermal of m4_nitrogen) with every
# default of 2026-09-25 in force (charge-exchange heating, He I 2^3S coupling
# with the stratified FUV band from 4.78 eV, collider-resolved two-photon
# quenching, He I/He II on-the-spot photons, He II diffuse channel with the
# transported continuum, Mg + H+ branching, diffuse packet normalization),
# restarted from the m5 (m4) converged states.  Three 24-rank runs at a time.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m8_defaults
one() { echo "$1 start $(date +%T)" >> driver.log; mpirun -np 24 ./MoCHII_m8.x $1.in > $1.log 2>&1; echo "$1 exit=$? $(date +%T)" >> driver.log; }
export -f one
printf '%s\n' hii40_converged hii20_converged hei_high_density_converged \
  hii40_ots_converged hii20_ots_converged pdr_L5_converged \
  hii40_lowden_converged hii20_lowden_converged hii40_nohecx_converged \
  hii20_nohecx_converged g2a_thermal_defaults \
  | xargs -P 3 -I{} bash -c 'one {}'
echo "ALL DONE $(date +%T)" >> driver.log
