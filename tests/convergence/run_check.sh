#!/bin/bash
#--- Build and run the stopping-test gate against the compiled module.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/fixed_point_distance_mod.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
$FC -O2 -ipo -module src -Isrc -o tests/convergence/check_fixed_point.x \
    tests/convergence/check_fixed_point.f90 src/define.o src/fixed_point_distance_mod.o
./tests/convergence/check_fixed_point.x
