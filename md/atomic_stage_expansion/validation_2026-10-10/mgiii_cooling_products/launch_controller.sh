#!/usr/bin/env bash
set -euo pipefail

if [[ "${MGIII_PRODUCTS_APPROVED:-}" != "yes" ]]; then
    echo 'Explicit Mg III product validation approval is required.' >&2
    exit 2
fi
if [[ -z "${MGIII_CORRECTED_CANDIDATE:-}" || -z "${MGIII_VALIDATION_MANIFEST:-}" ]]; then
    echo 'Corrected candidate and post-correction validation manifest are required.' >&2
    exit 2
fi

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "${here}/prepare_plan.py" --candidate "${MGIII_CORRECTED_CANDIDATE}" \
    --validation "${MGIII_VALIDATION_MANIFEST}"
tmux new-session -d -s mochii_mgiii_cooling_products \
    "python3 ${HOME}/.codex/scripts/calculation_workflow.py --plan ${here}/plan.json --state-dir ${here}/controller"
tmux list-sessions | rg mochii_mgiii_cooling_products
