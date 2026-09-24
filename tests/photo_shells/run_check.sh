#!/bin/bash
#--- Build and run the metal photoionization gate against the compiled MoCHII
#--- objects and Verner's phfit2.f (the source of every PHOTO/PHSHELL row).
#--- Exits nonzero when the gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/define.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
HDF5_PREFIX=${HDF5_PREFIX:-/data/opt/hdf5_intel}
SEDUST_LIBDIR=${SEDUST_LIBDIR:-SEDust/sed/lib}
VERNER_DIR=${VERNER_DIR:-$HOME/RT_Codes/Verner_Fortran_sub}
if [ ! -f "$VERNER_DIR/phfit2.f" ]; then
   echo "ERROR: $VERNER_DIR/phfit2.f not found (set VERNER_DIR)" >&2
   exit 1
fi
OBJS=$(ls src/*.o | grep -v '/main\.o$')
ifort -O0 -c "$VERNER_DIR/phfit2.f" -o tests/photo_shells/phfit2.o
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module src -Isrc -I"$SEDUST_LIBDIR" \
    -I"$HDF5_PREFIX/include" \
    -o tests/photo_shells/check_photo_shells.x \
    tests/photo_shells/check_photo_shells.f90 tests/photo_shells/phfit2.o $OBJS \
    -lcfitsio -L/usr/local/lib \
    "$HDF5_PREFIX/lib/libhdf5_fortran.a" "$HDF5_PREFIX/lib/libhdf5.a" \
    -lsz -ldl -lz -lm "$SEDUST_LIBDIR/libsedust.a" -qopenmp
./tests/photo_shells/check_photo_shells.x
