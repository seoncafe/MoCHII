#!/bin/bash
#--- Build and run the He I cascade gates (2^1S two-photon spectrum, and the
#--- branch fractions) against the compiled MoCHII objects.  Exits nonzero when
#--- a gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/define.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
#--- -ipo: the production objects are Intel IR, so the gate must be compiled
#--- and linked with the same interprocedural setting.
$FC -cpp -DMPI -ipo -O2 -module src -Isrc -o tests/hei_cascade/check_hei_2ph.x \
    tests/hei_cascade/check_hei_2ph.f90 \
    src/define.o src/utility.o src/hei_twophoton_mod.o
$FC -cpp -DMPI -ipo -O2 -module src -Isrc -o tests/hei_cascade/check_hei_branches.x \
    tests/hei_cascade/check_hei_branches.f90 \
    src/define.o src/utility.o src/photo_xsec.o src/recomb_mod.o \
    src/hei_cascade_mod.o
$FC -cpp -DMPI -ipo -O2 -module src -Isrc -o tests/hei_cascade/check_hei_584.x \
    tests/hei_cascade/check_hei_584.f90 \
    src/define.o src/utility.o src/photo_xsec.o src/recomb_mod.o \
    src/hei_twophoton_mod.o src/hei_cascade_mod.o
cd tests/hei_cascade && ./check_hei_2ph.x && ./check_hei_branches.x && ./check_hei_584.x
