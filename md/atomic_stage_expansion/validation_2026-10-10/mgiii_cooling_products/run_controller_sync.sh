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
python3 "${HOME}/.codex/scripts/calculation_workflow.py" \
    --plan "${here}/plan.json" --state-dir "${here}/controller"
python3 - "${here}/candidate/result.json" "${here}/controller/state.json" <<'PY'
import json
import sys
from pathlib import Path
product = json.loads(Path(sys.argv[1]).read_text())
state = json.loads(Path(sys.argv[2]).read_text())
if product.get('passed') is not True or state.get('status') != 'completed':
    raise SystemExit('Mg III products or controller did not pass')
if any(step.get('status') != 'completed' for step in state['steps']):
    raise SystemExit('Mg III product sequence is incomplete')
print(sys.argv[1])
PY
