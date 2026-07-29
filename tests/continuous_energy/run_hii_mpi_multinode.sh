#!/usr/bin/env bash
# Three-node continuous HII MPI smoke: lart2/lart3/lart4, 68 ranks each.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo"
results="$repo/results/continuous_energy/mpi_multinode"
mkdir -p "$results"
exe=${MOCHII_EXE:-$repo/MoCHII.x}
hosts=${MOCHII_HOSTS:-lart2,lart3,lart4}
ppn=${MOCHII_PPN:-68}
total=$((3 * ppn))
cases=("$@")
if [ ${#cases[@]} -eq 0 ]; then
    cases=(hii20 hii40)
fi

for tag in "${cases[@]}"; do
    input="tests/continuous_energy/${tag}_continuous.in"
    # The temporary input must live on the shared repository filesystem so
    # every remote MPI rank can open it.
    work=$(mktemp "$repo/tests/continuous_energy/.mochii_${tag}_multinode_XXXXXX.in")
    trap 'rm -f "$work"' EXIT
    sed \
        -e 's/par%no_photons      = 33554432/par%no_photons      = 100000/' \
        -e 's/par%no_print        = 33554432/par%no_print        = 100000/' \
        -e 's/par%gas_niter       = 100/par%gas_niter       = 3/' \
        -e 's/par%require_convergence = .true./par%require_convergence = .false./' \
        -e "s#${tag}_continuous.h5#${results}/${tag}_mpi_multinode.h5#" \
        "$input" > "$work"
    echo "=== ${tag}: hosts=${hosts}, ppn=${ppn}, total=${total} ==="
    I_MPI_HYDRA_BOOTSTRAP=ssh /usr/bin/time -f 'elapsed_s=%e max_rss_kb=%M' \
        -o "$results/${tag}_mpi_multinode.time" \
        mpirun -hosts "$hosts" -ppn "$ppn" -np "$total" "$exe" "$work" \
        > "$results/${tag}_mpi_multinode.log" 2>&1
    tail -n 4 "$results/${tag}_mpi_multinode.log"
    cat "$results/${tag}_mpi_multinode.time"
    rm -f "$work"
    trap - EXIT
done
