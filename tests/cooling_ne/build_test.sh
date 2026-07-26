#!/bin/bash
# Build the standalone (T, n_e) cooling-table gate (G-a interpolation,
# G-c low-density identity).  Links only define + nlevel_mod +
# nlevel_cooling_mod, which use nothing outside iso_fortran_env, so plain
# ifort links (no MPI, no HDF5).
set -e
cd "$(dirname "$0")"
SRC=../../src
FC=${FC:-ifort}
$FC -O2 -module . -o test_cool_table \
    "$SRC/define.f90" "$SRC/nlevel_mod.f90" "$SRC/nlevel_cooling_mod.f90" \
    test_cool_table.f90
echo "built ./test_cool_table"
