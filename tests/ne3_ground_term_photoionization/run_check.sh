#!/bin/bash
#--- Build and run the Ne III tabulated ground-term cross-section gate
#--- against the compiled MoCHII objects.  reference.txt is written by
#--- make_reference.py from data/atomic/norad/ne3.px.gd.txt (rerun it when
#--- the table is regenerated).  Exits nonzero when the gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/define.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
HDF5_PREFIX=${HDF5_PREFIX:-/data/opt/hdf5_intel}
SEDUST_LIBDIR=${SEDUST_LIBDIR:-SEDust/sed/lib}
D=tests/ne3_ground_term_photoionization
OBJS=$(ls src/*.o | grep -v '/main\.o$')
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module src -Isrc -I"$SEDUST_LIBDIR" \
    -I"$HDF5_PREFIX/include" \
    -o $D/check_ne3_ground_term_photoionization.x \
    $D/check_ne3_ground_term_photoionization.f90 $OBJS \
    -lcfitsio -L/usr/local/lib \
    "$HDF5_PREFIX/lib/libhdf5_fortran.a" "$HDF5_PREFIX/lib/libhdf5.a" \
    -lsz -ldl -lz -lm "$SEDUST_LIBDIR/libsedust.a" -qopenmp -qmkl=sequential
./$D/check_ne3_ground_term_photoionization.x
