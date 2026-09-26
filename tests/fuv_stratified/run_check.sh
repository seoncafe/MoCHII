#!/bin/bash
#--- Gate: band-stratified launch of the FUV extension (par%fuv_sampling).
#--- Converges one state (fuv_stratified.in), then transports it with
#--- gas_niter = 0 in 'proportional', 'separate' and add_fuv = .false., four
#--- seeds each, for a 40 kK and a 20 kK blackbody, plus three Sobol runs;
#--- check_fuv_stratified.py evaluates gates (a) and (b).
#--- Usage: tests/fuv_stratified/run_check.sh [nranks]   (default 8)
#--- WORKDIR selects where the runs go (default: a fresh temporary directory).
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
NP=${1:-8}
W=${WORKDIR:-$(mktemp -d "${TMPDIR:-/tmp}/mochii_fuv_stratified.XXXXXX")}
mkdir -p "$W"
python3 "$R/tests/fuv_stratified/check_fuv_stratified.py" make \
        "$R/tests/fuv_stratified/fuv_stratified.in" "$W" "$R"
cd "$W"
mpirun -np "$NP" "$R/MoCHII.x" base.in > base.log 2>&1 < /dev/null
for name in $(cat runs.txt); do
   mpirun -np "$NP" "$R/MoCHII.x" "$name.in" > "$name.log" 2>&1 < /dev/null
done
python3 "$R/tests/fuv_stratified/check_fuv_stratified.py" check "$W"
