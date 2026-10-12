#!/bin/bash
#--- Build and run the Milne gates against the already-compiled MoCHII
#--- objects: the ground-continuum sampler, and the ground-level alpha_1 the
#--- same relation fixes.  Exits nonzero when a gate fails.
#--- All objects are linked: milne_recomb_spectrum_mod also carries the
#--- metal ground continua and so uses species_mod (2026-10-11).
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/define.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
HDF5_PREFIX=${HDF5_PREFIX:-/data/opt/hdf5_intel}
SEDUST_LIBDIR=${SEDUST_LIBDIR:-SEDust/sed/lib}
OBJS=$(ls src/*.o | grep -v '/main\.o$')
LIBS="-lcfitsio -L/usr/local/lib $HDF5_PREFIX/lib/libhdf5_fortran.a $HDF5_PREFIX/lib/libhdf5.a -lsz -ldl -lz -lm $SEDUST_LIBDIR/libsedust.a -qopenmp -qmkl=sequential"
#--- -ipo: the production objects are Intel IR, so the gate must be
#--- compiled and linked with the same interprocedural setting.
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module src -Isrc -I"$SEDUST_LIBDIR" \
    -o tests/milne/check_milne.x tests/milne/check_milne.f90 $OBJS $LIBS
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module src -Isrc -I"$SEDUST_LIBDIR" \
    -o tests/milne/check_alpha1.x tests/milne/check_alpha1.f90 $OBJS $LIBS
./tests/milne/check_milne.x
./tests/milne/check_alpha1.x
