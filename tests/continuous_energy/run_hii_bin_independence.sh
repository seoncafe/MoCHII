#!/usr/bin/env bash
# Full-production diagnostic-bin gate for continuous HII20/HII40.
#
# Run from any directory.  Each case uses one fixed 2^25-packet Sobol history;
# only nnu_ion changes, so the continuous solver/state arrays must be bitwise
# identical.  The default leaves four of the 72 hardware threads available.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo"

exe=${MOCHII_EXE:-./MoCHII.x}
ranks=${MOCHII_RANKS:-68}
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
    for nbin in 8 16 32 64 128; do
        work=$(mktemp "/tmp/mochii_${tag}_n${nbin}_XXXXXX.in")
        trap 'rm -f "$work"' EXIT
        sed \
            -e "s/par%nnu_ion         = 32/par%nnu_ion         = ${nbin}/" \
            -e "s/${tag}_continuous.h5/${tag}_continuous_n${nbin}.h5/" \
            "$input" > "$work"
        echo "=== ${tag}: nnu_ion=${nbin}, ranks=${ranks} ==="
        mpirun -np "$ranks" "$exe" "$work" | tee "${tag}_continuous_n${nbin}.log"
        rm -f "$work"
        trap - EXIT
    done
done
