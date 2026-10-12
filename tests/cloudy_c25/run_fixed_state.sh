#!/bin/bash
#--- Build tests/cloudy_c25/fixed_state_nlevel.x against the compiled MoCHII
#--- objects (src/*.o except main.o; run `make` at the repository root
#--- first) and run the fixed-state cooling comparison against Cloudy c25.00.
#--- Usage: tests/cloudy_c25/run_fixed_state.sh [options of fixed_state_cooling.py]
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
if [ ! -f src/define.o ]; then
   echo "src/*.o not found: run make at $R first" >&2
   exit 1
fi
FC=${FC:-mpiifort}
HDF5_PREFIX=${HDF5_PREFIX:-/data/opt/hdf5_intel}
SEDUST_LIBDIR=${SEDUST_LIBDIR:-SEDust/sed/lib}
OBJS=$(ls src/*.o | grep -v '/main\.o$')
BLD=$(mktemp -d)
trap 'rm -rf "$BLD"' EXIT
#--- -module "$BLD" keeps the driver's own .mod out of src/
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module "$BLD" -Isrc -I"$SEDUST_LIBDIR" \
    -I"$HDF5_PREFIX/include" \
    -o tests/cloudy_c25/fixed_state_nlevel.x \
    tests/cloudy_c25/fixed_state_nlevel.f90 $OBJS \
    -lcfitsio -L/usr/local/lib \
    "$HDF5_PREFIX/lib/libhdf5_fortran.a" "$HDF5_PREFIX/lib/libhdf5.a" \
    -lsz -ldl -lz -lm "$SEDUST_LIBDIR/libsedust.a" -qopenmp -qmkl=sequential
python3 tests/cloudy_c25/fixed_state_cooling.py "$@"
