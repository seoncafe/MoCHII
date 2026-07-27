#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/../.." && pwd)
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/mochii-io-test.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

hdf5_prefix=${HDF5_PREFIX:-/data/opt/hdf5_intel}
hdf5_fc=${HDF5_FC:-mpiifort}
plain_fc=${PLAIN_FC:-mpif90}

"$plain_fc" -cpp -J"$work_dir" -c "$repo_dir/src/hdf5io_mod.f90" \
  -o "$work_dir/hdf5io_nohdf5.o"

"$hdf5_fc" -cpp -DHDF5 -I"$hdf5_prefix/include" \
  "$repo_dir/src/fitsio_mod.f90" "$repo_dir/src/hdf5io_mod.f90" \
  "$repo_dir/src/iofile_mod.f90" "$repo_dir/tests/io/hdf5_status_test.f90" \
  "$hdf5_prefix/lib/libhdf5_fortran.a" "$hdf5_prefix/lib/libhdf5.a" \
  -lcfitsio -lsz -ldl -lz -lm -o "$work_dir/hdf5_status_test.x"

"$work_dir/hdf5_status_test.x" \
  "$work_dir/status_test.h5" "$work_dir/missing/status_test.h5"
