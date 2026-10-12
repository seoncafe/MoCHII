#!/usr/bin/env bash
set -euo pipefail

if [[ "${OIV_CONDITIONAL_DIAGNOSTIC_APPROVED:-}" != "yes" ]]; then
    echo 'Root review of the conditional O IV diagnostic is required.' >&2
    exit 2
fi
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
session='mochii_oiv_full_20261010'
if tmux -L mochii-workflows has-session -t "${session}" 2>/dev/null; then
    echo "Existing ${session} tmux session must be inspected." >&2
    exit 2
fi
if [[ -e "${here}/controller" || -e "${here}/candidate" ]]; then
    echo 'Existing O IV controller or product must be preserved and inspected.' >&2
    exit 2
fi
if [[ ! -e "${here}/plan.json" ]]; then
    python3 "${here}/prepare_plan.py"
fi
tmux -L mochii-workflows new-session -d -s "${session}" \
    "python3 '${HOME}/.codex/scripts/calculation_workflow.py' --plan '${here}/plan.json' --state-dir '${here}/controller' > '${here}/controller_launcher.log' 2>&1"
tmux -L mochii-workflows list-panes -t "${session}" -F '#{pane_pid} #{pane_start_command}' > "${here}/tmux_pane.txt"
printf '%s\n' "${session}" > "${here}/tmux_session.txt"
printf 'Started %s; inspect controller/state.json and controller/progress.md.\n' "${session}"
