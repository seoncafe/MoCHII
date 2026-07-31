#!/bin/bash
#--- Build and run the two-way charge-exchange gate against the already-compiled
#--- MoCHII objects: every charge-exchange reaction must carry the same
#--- volumetric rate in hydrogen's balance and in the metal cascade.  Exits
#--- nonzero when the gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$R"
[ -f src/define.o ] || make MoCHII.x > /dev/null
FC=${FC:-mpiifort}
HDF5_PREFIX=${HDF5_PREFIX:-/data/opt/hdf5_intel}
SEDUST_LIBDIR=${SEDUST_LIBDIR:-SEDust/sed/lib}
#--- the gate calls solve_ion_cell, so it links the compiled objects (all but
#--- main.o).  -ipo: those objects are Intel IR, so the gate must be compiled
#--- and linked with the same interprocedural setting.
OBJS=$(ls src/*.o | grep -v '/main\.o$')
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module src -Isrc -I"$SEDUST_LIBDIR" \
    -I"$HDF5_PREFIX/include" \
    -o tests/charge_exchange/check_charge_exchange.x \
    tests/charge_exchange/check_charge_exchange.f90 $OBJS \
    -lcfitsio -L/usr/local/lib \
    "$HDF5_PREFIX/lib/libhdf5_fortran.a" "$HDF5_PREFIX/lib/libhdf5.a" \
    -lsz -ldl -lz -lm "$SEDUST_LIBDIR/libsedust.a" -qopenmp
./tests/charge_exchange/check_charge_exchange.x

#--- G-M1f: the cell solve must converge inside one call in partially neutral
#--- gas, where the lagged substitution map is marginally stable.
$FC -cpp -DMPI -DHDF5 -ipo -O2 -module src -Isrc -I"$SEDUST_LIBDIR" \
    -I"$HDF5_PREFIX/include" \
    -o tests/charge_exchange/check_cell_convergence.x \
    tests/charge_exchange/check_cell_convergence.f90 $OBJS \
    -lcfitsio -L/usr/local/lib \
    "$HDF5_PREFIX/lib/libhdf5_fortran.a" "$HDF5_PREFIX/lib/libhdf5.a" \
    -lsz -ldl -lz -lm "$SEDUST_LIBDIR/libsedust.a" -qopenmp
./tests/charge_exchange/check_cell_convergence.x
