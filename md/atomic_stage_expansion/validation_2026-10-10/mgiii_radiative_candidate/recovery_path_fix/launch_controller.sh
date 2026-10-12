#!/usr/bin/env bash
set -euo pipefail
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00
exec /opt/miniconda3/bin/python /nfs/mocafe/kiseon/.codex/scripts/calculation_workflow.py --plan /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/md/atomic_stage_expansion/validation_2026-10-10/mgiii_radiative_candidate/recovery_path_fix/plan.json --state-dir /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/md/atomic_stage_expansion/validation_2026-10-10/mgiii_radiative_candidate/recovery_path_fix/controller >> /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/md/atomic_stage_expansion/validation_2026-10-10/mgiii_radiative_candidate/recovery_path_fix/controller_main.log 2>&1
