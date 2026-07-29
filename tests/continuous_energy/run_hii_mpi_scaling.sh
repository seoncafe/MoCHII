#!/usr/bin/env bash
# MPI-rank agreement and lightweight scaling smoke for continuous HII.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo"
exe=${MOCHII_EXE:-./MoCHII.x}
ranks=(1 2 3 5 8 68)
cases=("$@")
if [ ${#cases[@]} -eq 0 ]; then
    cases=(hii20 hii40)
fi

for tag in "${cases[@]}"; do
    input="tests/continuous_energy/${tag}_continuous.in"
    for rank in "${ranks[@]}"; do
        work=$(mktemp "/tmp/mochii_${tag}_mpi_${rank}_XXXXXX.in")
        trap 'rm -f "$work"' EXIT
        sed \
            -e 's/par%no_photons      = 33554432/par%no_photons      = 100000/' \
            -e 's/par%no_print        = 33554432/par%no_print        = 100000/' \
            -e 's/par%gas_niter       = 100/par%gas_niter       = 3/' \
            -e 's/par%require_convergence = .true./par%require_convergence = .false./' \
            -e "s/${tag}_continuous.h5/${tag}_mpi_np${rank}.h5/" \
            "$input" > "$work"
        echo "=== ${tag}: MPI ranks=${rank} ==="
        /usr/bin/time -f 'elapsed_s=%e max_rss_kb=%M' \
            -o "${tag}_mpi_np${rank}.time" \
            mpirun -np "$rank" "$exe" "$work" > "${tag}_mpi_np${rank}.log" 2>&1
        tail -n 4 "${tag}_mpi_np${rank}.log"
        cat "${tag}_mpi_np${rank}.time"
        rm -f "$work"
        trap - EXIT
    done
done
