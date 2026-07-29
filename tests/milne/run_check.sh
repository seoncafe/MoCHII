#!/bin/bash
#--- Build and run the Milne ground-continuum sampler gate against the
#--- already-compiled MoCHII objects.  Exits nonzero when the gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/define.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
#--- -ipo: the production objects are Intel IR, so the gate must be
#--- compiled and linked with the same interprocedural setting.
$FC -cpp -DMPI -ipo -O2 -module src -Isrc -o tests/milne/check_milne.x \
    tests/milne/check_milne.f90 \
    src/define.o src/photo_xsec.o src/recomb_mod.o \
    src/milne_recomb_spectrum_mod.o
./tests/milne/check_milne.x
