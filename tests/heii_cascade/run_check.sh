#!/bin/bash
#--- Build and run the He III case-B cascade gate (diffuse channel 5) against
#--- the compiled MoCHII objects.  Exits nonzero when the gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/define.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
#--- -ipo: the production objects are Intel IR, so the gate must be compiled
#--- and linked with the same interprocedural setting.
$FC -cpp -DMPI -ipo -O2 -module src -Isrc \
    -o tests/heii_cascade/check_heii_cascade.x \
    tests/heii_cascade/check_heii_cascade.f90 \
    src/define.o src/utility.o src/photo_xsec.o src/recomb_mod.o \
    src/gaunt.o src/gaunt_vh14_mod.o src/recomb_continuum_mod.o \
    src/hydrogenic_twophoton_mod.o src/heii_cascade_mod.o
cd tests/heii_cascade && ./check_heii_cascade.x
