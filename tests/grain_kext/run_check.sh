#!/bin/bash
#--- Build and run the grain extinction gate: the size-integrated extinction
#--- dust_extinction serves MoCHII from the model's own product (/kext of
#--- data/<model>/sedust_<model>.h5) must be the curve
#--- size_integrated_extinction computes for the same model, on the wavelength
#--- grids grain_model_mod actually builds.  Exits nonzero when the gate fails.
set -e
R="$(cd "$(dirname "$0")/../.." && pwd)"
FC=${FC:-ifort}
SEDUST_LIBDIR=${SEDUST_LIBDIR:-$R/SEDust/sed/lib}
#--- the gate needs only SEDust, not the MoCHII objects: it exercises the
#--- library API grain_model_mod calls, with grain_model_mod's arguments.
if [ ! -f "$SEDUST_LIBDIR/libsedust.a" ]; then
   echo "ERROR: $SEDUST_LIBDIR/libsedust.a is missing; run SEDust/sed/build_lib.sh" >&2
   exit 1
fi
#--- libsedust.a reads the optics products as HDF5, so the archive must come
#--- BEFORE the HDF5 libraries (a static archive resolves left to right), the
#--- same order the MoCHII Makefile links in.  A text-only archive
#--- (HDF5=0 ./build_lib.sh) needs none of them: HDF5_PREFIX= drops them.
HDF5_PREFIX=${HDF5_PREFIX-/data/opt/hdf5_intel}
HDF5_LIBS=""
if [ -n "$HDF5_PREFIX" ]; then
   HDF5_LIBS="$HDF5_PREFIX/lib/libhdf5_fortran.a $HDF5_PREFIX/lib/libhdf5.a -lsz -ldl -lz -lm"
fi
$FC -O2 -qopenmp -I"$SEDUST_LIBDIR" \
    -o "$R/tests/grain_kext/check_kext_table.x" \
    "$R/tests/grain_kext/check_kext_table.f90" "$SEDUST_LIBDIR/libsedust.a" $HDF5_LIBS
#--- build_dust resolves every file the model is made of -- optics product,
#--- size distribution, dielectric functions, ZDA config -- inside the data
#--- directory it is given, so the gate hands it the same one a run gets from
#--- par%sed_data_dir and needs no particular working directory.
"$R/tests/grain_kext/check_kext_table.x" "$R/SEDust/data"
