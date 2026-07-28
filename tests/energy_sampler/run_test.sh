#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"

FC=${FC:-mpiifort}
"$FC" -ipo -O2 -Isrc -module src tests/energy_sampler/energy_sampler_dump.f90 \
  src/energy_sampler_mod.o src/qmc_mod.o -o tests/energy_sampler/energy_sampler_dump.x
"$FC" -ipo -O2 -Isrc -module src tests/energy_sampler/energy_sampler_mpi.f90 \
  src/energy_sampler_mod.o src/qmc_mod.o -o tests/energy_sampler/energy_sampler_mpi.x
python3 tests/energy_sampler/check_energy_sampler.py
