#!/usr/bin/env bash
# Paired HII40 experiments for the residual C III / He I-front comparison.
#
# Each case starts from the accepted continuous-energy production input and
# changes exactly one physical or numerical ingredient.  Outputs are ignored
# generated data under results/continuous_energy/residual_investigation/.
#
# Usage:
#   MOCHII_RANKS=68 ./tests/continuous_energy/run_hii40_residual_experiments.sh \
#       diffuse_off no_hei l7 mocrec threshold milne
#
# `mocrec` needs MOCREC_ATOMIC_DIR to point to the alternate atomic-data tree.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo"

base="tests/continuous_energy/hii40_continuous.in"
results="$repo/results/continuous_energy/residual_investigation"
exe=${MOCHII_EXE:-./MoCHII.x}
ranks=${MOCHII_RANKS:-68}
cases=("$@")

if [ ${#cases[@]} -eq 0 ]; then
    cases=(diffuse_off no_hei l7 mocrec threshold milne)
fi
mkdir -p "$results"

for case_name in "${cases[@]}"; do
    work=$(mktemp "/tmp/mochii_hii40_${case_name}_XXXXXX.in")
    trap 'rm -f "$work"' EXIT
    output="$results/hii40_${case_name}.h5"
    sed -e "s#hii40_continuous.h5#${output}#" "$base" > "$work"

    case "$case_name" in
    diffuse_off)
        sed -i \
            -e "s/par%case_ab         = 'A'/par%case_ab         = 'B'/" \
            -e 's/par%diffuse_field   = .true./par%diffuse_field   = .false./' \
            -e 's/par%hei_diffuse     = .true./par%hei_diffuse     = .false./' \
            "$work"
        ;;
    no_hei)
        sed -i 's/par%hei_diffuse     = .true./par%hei_diffuse     = .false./' "$work"
        ;;
    l7)
        sed -i 's#tests/g2_hii/hii40_L6.fits#tests/g2_hii/hii40_L7.fits#' "$work"
        ;;
    mocrec)
        : "${MOCREC_ATOMIC_DIR:?set MOCREC_ATOMIC_DIR for the mocrec case}"
        sed -i "s#par%atomic_dir      = 'data/atomic'#par%atomic_dir      = '${MOCREC_ATOMIC_DIR}'#" "$work"
        ;;
    threshold)
        # Sensitivity bracket for the ground-recombination continuum energy.
        sed -i "/par%diffuse_field/a par%diffuse_energy_model = 'threshold'" "$work"
        ;;
    milne)
        # Detailed-balance free-bound ground continuum: the physical model the
        # exponential default approximates.
        sed -i "/par%diffuse_field/a par%diffuse_energy_model = 'milne'" "$work"
        ;;
    *)
        echo "unknown residual experiment: $case_name" >&2
        exit 2
        ;;
    esac

    cp "$work" "$results/hii40_${case_name}.in"
    echo "=== HII40 residual experiment: ${case_name}, ranks=${ranks} ==="
    mpirun -np "$ranks" "$exe" "$work" | tee "$results/hii40_${case_name}.log"
    # The legacy output layer retains only basename(par%out_file), so collect
    # its three root-level products after a successful run.
    for product in rates.h5 lines.txt nebcont.txt; do
        generated="hii40_${case_name}_${product}"
        if [ -f "$generated" ]; then
            mv "$generated" "$results/"
        fi
    done
    rm -f "$work"
    trap - EXIT
done
