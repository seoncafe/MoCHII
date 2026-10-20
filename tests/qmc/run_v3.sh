#!/usr/bin/env bash
# Continuous-energy fixed-state MC/RQMC study. Preserve historical v3 outputs.
# The controller validates each run and automatically analyzes the ensemble.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"
python tests/qmc/run_continuous.py prepare
OUT="$ROOT/results/qmc_continuous_2026-10-18"
tmux new-session -d -s qmc_continuous_20261018 \
  python "$HOME/.codex/scripts/calculation_workflow.py" \
  --plan "$OUT/plan.json" --state-dir "$OUT/controller"
echo "Started qmc_continuous_20261018; progress: $OUT/controller/progress.md"
