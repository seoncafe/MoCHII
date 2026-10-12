#!/usr/bin/env bash
set -euo pipefail
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00
exec /opt/miniconda3/bin/python /nfs/mocafe/kiseon/.codex/scripts/calculation_workflow.py \
    --plan /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/md/atomic_stage_expansion/validation_2026-10-10/si_i_radiative_completion/plan.json \
    --state-dir /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/md/atomic_stage_expansion/validation_2026-10-10/si_i_radiative_completion/controller \
    >> /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/md/atomic_stage_expansion/validation_2026-10-10/si_i_radiative_completion/controller_main.log 2>&1
