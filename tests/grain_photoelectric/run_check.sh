#!/bin/bash
#--- Build and run the grain photoelectric / gas-grain gate against the
#--- compiled MoCHII objects (make MoCHII.x first) and SEDust.  Exits nonzero
#--- when the gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/grain_photoelectric_mod.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
HDF5_PREFIX=${HDF5_PREFIX:-/data/opt/hdf5_intel}
SEDUST_LIBDIR=${SEDUST_LIBDIR:-SEDust/sed/lib}
OBJS=$(ls src/*.o | grep -v '/main\.o$')
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module tests/grain_photoelectric -Isrc \
    -I"$SEDUST_LIBDIR" -I"$HDF5_PREFIX/include" \
    -o tests/grain_photoelectric/check_grain_photoelectric.x \
    tests/grain_photoelectric/check_grain_photoelectric.f90 $OBJS \
    -lcfitsio -L/usr/local/lib "$SEDUST_LIBDIR/libsedust.a" \
    "$HDF5_PREFIX/lib/libhdf5_fortran.a" "$HDF5_PREFIX/lib/libhdf5.a" \
    -lsz -ldl -lz -lm -qopenmp -qmkl=sequential
./tests/grain_photoelectric/check_grain_photoelectric.x "$R/SEDust/data"
