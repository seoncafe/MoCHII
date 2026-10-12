#!/usr/bin/env bash
set -euo pipefail

if [[ "${NEI_PRODUCTS_APPROVED:-}" != "yes" ]]; then
    echo 'Explicit Ne I product validation approval is required.' >&2
    exit 2
fi
if [[ -z "${NEI_CORRECTED_CANDIDATE:-}" || -z "${NEI_VALIDATION_MANIFEST:-}" ]]; then
    echo 'Corrected candidate and post-correction validation manifest are required.' >&2
    exit 2
fi

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "${here}/prepare_plan.py" --candidate "${NEI_CORRECTED_CANDIDATE}" \
    --validation "${NEI_VALIDATION_MANIFEST}"
tmux new-session -d -s mochii_nei_cooling_products \
    "python3 ${HOME}/.codex/scripts/calculation_workflow.py --plan ${here}/plan.json --state-dir ${here}/controller"
tmux list-sessions | rg mochii_nei_cooling_products
