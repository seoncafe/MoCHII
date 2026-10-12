#!/bin/bash
#--- Build and run the case-B on-the-spot He I gate against the compiled MoCHII
#--- objects (it calls solve_ion_cell, so it
#--- links all of them but main.o).  Exits nonzero when the gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/define.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
HDF5_PREFIX=${HDF5_PREFIX:-/data/opt/hdf5_intel}
SEDUST_LIBDIR=${SEDUST_LIBDIR:-SEDust/sed/lib}
#--- -ipo: the production objects are Intel IR, so the gate must be compiled
#--- and linked with the same interprocedural setting.
OBJS=$(ls src/*.o | grep -v '/main\.o$')
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module src -Isrc -I"$SEDUST_LIBDIR" \
    -I"$HDF5_PREFIX/include" \
    -o tests/hei_on_the_spot/check_hei_on_the_spot.x \
    tests/hei_on_the_spot/check_hei_on_the_spot.f90 $OBJS \
    -lcfitsio -L/usr/local/lib \
    "$HDF5_PREFIX/lib/libhdf5_fortran.a" "$HDF5_PREFIX/lib/libhdf5.a" \
    -lsz -ldl -lz -lm "$SEDUST_LIBDIR/libsedust.a" -qopenmp -qmkl=sequential
cd tests/hei_on_the_spot && ./check_hei_on_the_spot.x
