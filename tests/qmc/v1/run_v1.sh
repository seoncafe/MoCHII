#!/usr/bin/env bash
# V1 regression + diffuse-fix impact.  Runs the HEAD (pre-QMC) and NEW binaries
# (both built -fp-model precise) on the same inputs and leaves _head/_new rates
# files for check_v1.py.  (a) multi-source and (b) external at np=1 -> expect
# bit-identical random path; (c) diffuse LOG grid at np=8 -> Task-0 bin fix
# inert, expect bit-identical; (d) diffuse ALIGNED grid -> the fix changes
# results (the point), quantified by the checker.
set -e
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
cd "$ROOT"
HEAD=/tmp/claude-1000/MoCHII_head_precise.x
NEW=/tmp/claude-1000/MoCHII_new_precise.x

run () {  # $1=input(no ext) $2=np    (rates written as <basename>_rates.h5 in cwd)
  local in="tests/qmc/v1/$1.in"
  mpirun -np "$2" "$HEAD" "$in" > "tests/qmc/v1/$1_head.log" 2>&1
  mv "$1_rates.h5" "tests/qmc/v1/$1_head_rates.h5"
  mpirun -np "$2" "$NEW"  "$in" > "tests/qmc/v1/$1_new.log" 2>&1
  mv "$1_rates.h5" "tests/qmc/v1/$1_new_rates.h5"
  echo "  done $1 (np=$2)"
}

echo "V1(a) multi-source, random, np=1";  run v1a_multi 1
echo "V1(b) external,     random, np=1";  run v1b_ext   1
echo "V1(c) diffuse LOG,  random, np=8";  run v1c_diff_logedge 8
echo "V1(d) diffuse ALIGNED, random, np=8"; run v1d_diff_aligned 8
echo "V1 runs complete."
