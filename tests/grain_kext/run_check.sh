#!/bin/bash
#--- Build and run the grain extinction-table gate: the size-integrated
#--- extinction dust_extinction serves MoCHII from data/kext_*.dat must be the
#--- curve size_integrated_extinction computes for the same model, on the
#--- wavelength grids grain_model_mod actually builds.  Exits nonzero when the
#--- gate fails.
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
$FC -O2 -qopenmp -I"$SEDUST_LIBDIR" \
    -o "$R/tests/grain_kext/check_kext_table.x" \
    "$R/tests/grain_kext/check_kext_table.f90" "$SEDUST_LIBDIR/libsedust.a"
#--- SEDust names its Q tables, dielectric functions and kext tables by paths
#--- relative to sed/, and so does grain_model_mod (it chdir's to
#--- par%sed_workdir before building), so the gate runs there too.
cd "$R/SEDust/sed"
"$R/tests/grain_kext/check_kext_table.x"
