#!/usr/bin/env bash
# Independent scrambled-Sobol production ensemble for HII20/HII40.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo"
exe=${MOCHII_EXE:-./MoCHII.x}
ranks=${MOCHII_RANKS:-68}
seeds=(101 202 303 404)
cases=("$@")
if [ ${#cases[@]} -eq 0 ]; then
    cases=(hii20 hii40)
fi

for tag in "${cases[@]}"; do
    input="tests/continuous_energy/${tag}_continuous.in"
    if [ ! -f "$input" ]; then
        echo "unknown production case: $tag" >&2
        exit 2
    fi
    for seed in "${seeds[@]}"; do
        work=$(mktemp "/tmp/mochii_${tag}_rqmc_${seed}_XXXXXX.in")
        trap 'rm -f "$work"' EXIT
        sed \
            -e "s/par%qmc_seed        = 20260728/par%qmc_seed        = ${seed}/" \
            -e "s/${tag}_continuous.h5/${tag}_rqmc_${seed}.h5/" \
            "$input" > "$work"
        echo "=== ${tag}: qmc_seed=${seed}, ranks=${ranks} ==="
        mpirun -np "$ranks" "$exe" "$work" | tee "${tag}_rqmc_${seed}.log"
        rm -f "$work"
        trap - EXIT
    done
done
