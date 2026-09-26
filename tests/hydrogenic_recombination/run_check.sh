#!/bin/bash
#--- Build and run the hydrogenic recombination gate against the compiled
#--- MoCHII objects.  Exits nonzero when the gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/hydrogenic_recombination_mod.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
#--- -ipo: the production objects are Intel IR, so the gate must be compiled
#--- and linked with the same interprocedural setting.
$FC -cpp -DMPI -ipo -O2 -module src -Isrc \
    -o tests/hydrogenic_recombination/check_hydrogenic_recombination.x \
    tests/hydrogenic_recombination/check_hydrogenic_recombination.f90 \
    src/define.o src/utility.o src/hydrogenic_recombination_mod.o src/recomb_mod.o
cd tests/hydrogenic_recombination && ./check_hydrogenic_recombination.x
