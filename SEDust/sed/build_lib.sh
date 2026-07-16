#!/bin/bash
# Build the SEDust dust-emission library for MoCafe v2.00.
#
# This is a self-contained, Intel-built copy of the SEDust `dust_lib` layer
# (source: /home/kiseon/MoCafe/Grain/SEDust on the author's machine).  The module originally
# named `mathlib` is renamed here to `sed_mathlib` to avoid a name clash with
# MoCafe's own `mathlib` module.  Outputs sed/lib/libsedust.a and the .mod
# files that MoCafe links against (Makefile variable SEDUST_INTEL).
#
# Re-run after updating the sources in sed/src/.
set -e
cd "$(dirname "$0")"                 # sed/
FC=${FC:-ifort}
FLAGS="-O3 -xHost -qopenmp -w -module lib"
SRCS="constants mathlib enthalpy_v2 size_dist q_table enthalpy_astrodust mie \
      q_graphite q_graphite_d16 q_graphite_d16_sphere qpah radfield p_sub \
      stoch_qm pah_ioniz grain_dist q_silicate pah_ld01 dust_model_mod \
      zubko_io sed_astrodust dust_lib"
mkdir -p lib
rm -f lib/*.o lib/*.mod lib/libsedust.a
for f in $SRCS; do
   echo "  FC $f"
   $FC $FLAGS -c src/$f.f90 -o lib/$f.o
done
ar rcs lib/libsedust.a lib/*.o
echo "Built: sed/lib/libsedust.a"
