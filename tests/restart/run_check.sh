#!/bin/bash
#--- Gate: par%restart_file (previous_state_mod).  A short HII20 run writes a
#--- state; a restart on the same grid with gas_niter = 0 must return it bit
#--- for bit, and a restart on a 32^3 grid must give every gas leaf the values
#--- of the previous leaf containing its centre (check_restart.py).
#--- Usage: tests/restart/run_check.sh [nranks]   (default 8)
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
NP=${1:-8}
W=$(mktemp -d "${TMPDIR:-/tmp}/mochii_restart.XXXXXX")
cd "$W"
python3 "$R/tests/g2_hii/make_hii_grid.py" --rin 0.97223 --rout 3.07874 --level 5 --out coarse.fits > /dev/null
sed -e "s#'../../../#'$R/#g" -e "s#'../../#'$R/#g" \
    -e "s/par%no_photons *= *[0-9]*/par%no_photons = 262144/" \
    -e "s/par%require_convergence *= *.true./par%require_convergence = .false./" \
    "$R/tests/restart/hii20_restart.in" > base.in
sed -e "s/OUTNAME/fresh/" -e "s/NITER/2/" base.in > fresh.in
sed -e "s/OUTNAME/same/"  -e "s/NITER/0/" base.in > same.in
sed -e "s/OUTNAME/coarse/" -e "s/NITER/0/" -e "s#$R/tests/g2_hii/hii20_L6.fits#$W/coarse.fits#" base.in > coarse.in
for f in same coarse; do echo " par%restart_file = '$W/fresh_rates.h5'" > tail.txt
  sed -i "/par%out_file/r tail.txt" $f.in; done
for f in fresh same coarse; do mpirun -np $NP "$R/MoCHII.x" $f.in > $f.log 2>&1; done
python3 "$R/tests/restart/check_restart.py" "$W"
