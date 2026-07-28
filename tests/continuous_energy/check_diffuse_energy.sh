#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")/../.." && pwd)
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/mochii-diffuse-energy.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT
fc=${FC:-ifort}

"$fc" -O2 -module "$work_dir" -I"$work_dir" \
  "$repo_dir/src/ion_energy_policy_mod.f90" \
  "$repo_dir/tests/continuous_energy/diffuse_energy_policy_test.f90" \
  -o "$work_dir/diffuse_energy_policy_test.x"
"$work_dir/diffuse_energy_policy_test.x"
