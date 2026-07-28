#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 SERIAL_RATES.h5 PARALLEL_RATES.h5" >&2
    exit 2
fi

command -v h5diff >/dev/null || {
    echo "ERROR: h5diff is required" >&2
    exit 2
}

serial_file=$1
parallel_file=$2

for field in x_HI x_HeI x_HeII n_e T_e; do
    h5diff "$serial_file" "$parallel_file" \
        "/${field}/data" "/${field}/data"
done

echo "PASS: thermal gas state is exactly identical"
