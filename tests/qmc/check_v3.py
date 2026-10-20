#!/usr/bin/env python3
"""Analyze the current continuous-energy MC/RQMC ensemble.

Historical bin-center results remain in v3/. They require their historical
reference and are not analyzed by this command.
"""
from pathlib import Path
import subprocess
import sys

here = Path(__file__).resolve().parent
output = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else (
    here.parents[1] / "results/qmc_continuous_2026-10-18")
subprocess.run([sys.executable, str(here / "check_continuous.py"), str(output)], check=True)
