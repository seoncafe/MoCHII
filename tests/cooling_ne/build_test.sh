#!/bin/bash
# Build the standalone (T, n_e) cooling-table gate (G-a interpolation,
# G-c low-density identity, G-e runtime combiner). Links only define,
# nlevel_mod, and nlevel_cooling_mod (no MPI or HDF5). Keep all build products
# in an explicit isolated directory so an existing build is not overwritten.
set -euo pipefail
if [[ $# -ne 1 ]]; then
    echo "usage: $0 OUTPUT_DIRECTORY" >&2
    exit 2
fi
TEST_DIR=$(cd "$(dirname "$0")" && pwd)
SRC_DIR=$(cd "$TEST_DIR/../../src" && pwd)
mkdir -p "$1"
BUILD_DIR=$(cd "$1" && pwd)
FC=${FC:-ifort}
cd "$BUILD_DIR"
"$FC" -cpp -O2 -no-ipo -fp-model precise -check bounds -fpe0 -traceback \
    -module "$BUILD_DIR" -I"$BUILD_DIR" -o test_cool_table \
    "$SRC_DIR/define.f90" "$SRC_DIR/nlevel_mod.f90" \
    "$SRC_DIR/nlevel_cooling_mod.f90" "$TEST_DIR/test_cool_table.f90" -qmkl=sequential
echo "built $BUILD_DIR/test_cool_table"
