#!/bin/bash
#--- Build and run the metal ground-continuum gate (diffuse channel 6)
#--- against the compiled MoCHII objects and Verner's phfit2.f.
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
D=tests/metal_recombination_continuum
OBJS=$(ls src/*.o | grep -v '/main\.o$')
ifort -O0 -c "$VERNER_DIR/phfit2.f" -o $D/phfit2.o
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module src -Isrc -I"$SEDUST_LIBDIR" \
    -I"$HDF5_PREFIX/include" \
    -o $D/check_metal_recombination_continuum.x \
    $D/check_metal_recombination_continuum.f90 $D/phfit2.o $OBJS \
    -lcfitsio -L/usr/local/lib \
    "$HDF5_PREFIX/lib/libhdf5_fortran.a" "$HDF5_PREFIX/lib/libhdf5.a" \
    -lsz -ldl -lz -lm "$SEDUST_LIBDIR/libsedust.a" -qopenmp -qmkl=sequential
./$D/check_metal_recombination_continuum.x
