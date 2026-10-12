"""Check one installed oxygen ion through the optimized registered table."""

import argparse
import hashlib
import json
import math
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).parent
ROOT = HERE.parents[3]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("stage", type=int, choices=(2, 3))
args = parser.parse_args()
stage = args.stage
manifest = json.loads((HERE / "installation_manifest.json").read_text())
installed = Path(manifest[f"O{stage}"]["installed_path"])
if hashlib.sha256(installed.read_bytes()).hexdigest() != manifest[f"O{stage}"]["new_sha256"]:
    raise ValueError("Installed oxygen atom changed before registration")
profile = json.loads((HERE / "optimized_profile_result.json").read_text())
if not profile.get("passed") or profile["ions"][f"O{stage}"]["installed_sha256"] != manifest[f"O{stage}"]["new_sha256"]:
    raise ValueError("Optimized direct profile was not validated")
driver = HERE / "optimized_build" / "registered_check.x"
table = HERE / f"registered_O{stage}.txt"
if table.exists():
    raise FileExistsError("Preserve existing registered result")
started = time.monotonic()
subprocess.run([str(driver), str(ROOT / "data/atomic"), "o", str(stage), str(table)],
               check=True, timeout=8600)
elapsed = time.monotonic() - started
lines = table.read_text().splitlines()
if not lines or lines[-1] != "# PASS":
    raise ValueError("Registered result lacks final PASS")
values = {}
for line in lines:
    if line.startswith("# ") and len(line.split()) >= 3:
        parts = line.split(maxsplit=2)
        values[parts[1]] = parts[2]
rows = [line for line in lines if line.startswith("grid ") or line.startswith("middle ")]
expected_levels, expected_collisions = ((196, 19110) if stage == 2 else (177, 15576))
if (len(rows) != 40 or int(values["rows"]) != 40
    or int(values["nlev"]) != expected_levels
    or int(values["nups"]) != expected_collisions
    or values["density_ratio_cooling"] != "T"
    or int(values["failures"]) != 0):
    raise ValueError("Registered structural or physical flag check failed")
for field, limit in (("worst_grid_relative_error", 1e-10),
                     ("worst_middle_relative_error", .005),
                     ("worst_apply_relative_error", .005)):
    if not math.isfinite(float(values[field])) or float(values[field]) > limit:
        raise ValueError(f"Registered table threshold failed: {field}")
result = {"passed": True, "stage": stage, "rows": len(rows),
          "installed_sha256": manifest[f"O{stage}"]["new_sha256"],
          "driver": "isolated Intel -O2 build from current source",
          "elapsed_wall_seconds": elapsed,
          "registration_wall_seconds": float(values["registration_wall_seconds"]),
          "worst_grid_relative_error": float(values["worst_grid_relative_error"]),
          "worst_middle_relative_error": float(values["worst_middle_relative_error"]),
          "worst_apply_relative_error": float(values["worst_apply_relative_error"]),
          "physical_source_temperature_range_K": [1000.0, 100000.0],
          "grid_samples_outside_source_range_are_numerical_only": True}
(HERE / f"registered_O{stage}_result.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
