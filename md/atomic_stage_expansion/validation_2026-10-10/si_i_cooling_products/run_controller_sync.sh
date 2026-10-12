#!/usr/bin/env bash
set -euo pipefail

if [[ "${SII_PRODUCTS_APPROVED:-}" != "yes" ]]; then
    echo 'Reviewed Si I numerical-product approval is required.' >&2
    exit 2
fi
if [[ -z "${SII_CORRECTED_CANDIDATE:-}" || -z "${SII_VALIDATION_MANIFEST:-}" ]]; then
    echo 'Corrected source atom and complete numerical manifest are required.' >&2
    exit 2
fi

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "${here}/prepare_plan.py" --candidate "${SII_CORRECTED_CANDIDATE}" \
    --validation "${SII_VALIDATION_MANIFEST}"
python3 "${HOME}/.codex/scripts/calculation_workflow.py" \
    --plan "${here}/plan.json" --state-dir "${here}/controller"
python3 - "${here}/candidate/result.json" "${here}/controller/state.json" <<'PY'
import json
import sys
from pathlib import Path

product = json.loads(Path(sys.argv[1]).read_text())
state = json.loads(Path(sys.argv[2]).read_text())
if product.get('passed') is not True or product.get('ion') != 'SiI':
    raise SystemExit('Si I numerical products did not pass')
if state.get('status') != 'completed' or any(
    step.get('status') != 'completed' for step in state['steps']
):
    raise SystemExit('Si I product sequence is incomplete')
print(sys.argv[1])
PY
