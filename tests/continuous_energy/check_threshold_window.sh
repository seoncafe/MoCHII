#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/../.." && pwd)
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/mochii-threshold-test.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

fc=${FC:-ifort}

"$fc" -O2 -cpp -module "$work_dir" -I"$work_dir" \
  "$repo_dir/src/define.f90" \
  "$repo_dir/src/photo_xsec.f90" \
  "$repo_dir/tests/continuous_energy/threshold_window_test.f90" \
  -o "$work_dir/threshold_window_test.x"

"$work_dir/threshold_window_test.x" "$repo_dir/data/atomic/element_c.txt"
