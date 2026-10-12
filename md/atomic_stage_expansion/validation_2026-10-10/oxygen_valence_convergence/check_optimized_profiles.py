"""Compare optimized installed-atom solves against debug validated profiles."""

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = HERE.parents[3]
manifest = json.loads((HERE / "installation_manifest.json").read_text())
grid = HERE / "native_temperature_grid.txt"
executable = HERE / "optimized_build" / "source_check.x"
report = {}
for stage in (2, 3):
    installed = Path(manifest[f"O{stage}"]["installed_path"])
    if hashlib.sha256(installed.read_bytes()).hexdigest() != manifest[f"O{stage}"]["new_sha256"]:
        raise ValueError("Installed oxygen atom changed after validation")
    output = HERE / f"optimized_O{stage}_profiles.txt"
    subprocess.run([str(executable), str(installed), str(grid), str(output)],
                   check=True, timeout=1800)
    measured = np.loadtxt(output)
    debug = np.loadtxt(HERE / f"native_O{stage}_base.txt")
    if measured.shape != debug.shape or not np.array_equal(measured[:, :2], debug[:, :2]):
        raise ValueError("Optimized and debug population grids differ")
    cooling = float(np.max(abs(measured[:, 2] / debug[:, 2] - 1)))
    population = float(np.max(abs(measured[:, 3:] - debug[:, 3:])))
    if cooling > 1e-9 or population > 1e-10:
        raise ValueError(f"O {stage} optimized/debug mismatch: {cooling}, {population}")
    report[f"O{stage}"] = {"rows": int(len(measured)),
                            "maximum_cooling_relative_error": cooling,
                            "maximum_population_absolute_error": population,
                            "installed_sha256": manifest[f"O{stage}"]["new_sha256"]}
result = {"passed": True, "ions": report}
(HERE / "optimized_profile_result.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
