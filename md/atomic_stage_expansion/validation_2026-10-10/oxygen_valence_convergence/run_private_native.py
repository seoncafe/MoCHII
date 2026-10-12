"""Run scoped native oxygen populations, cooling, and line-balance tests."""

import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
MANIFEST = json.loads((HERE / "candidate_manifest.json").read_text())
if not MANIFEST.get("passed"):
    raise ValueError("Candidate construction did not complete")
DRIVER = HERE / "native_source_check.x"
TEMPERATURE = np.geomspace(1000.0, 100000.0, 401)
DENSITY = (1e-6, 1e2, 1e6, 1e10)
grid = HERE / "native_temperature_grid.txt"
grid.write_text("0 401\n" + " ".join(format(t, ".17g") for t in TEMPERATURE) + "\n")
profiles = {"temperature_K": TEMPERATURE, "density_cm-3": np.asarray(DENSITY)}
results = {}

fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
for axis, stage in zip(axes, (2, 3)):
    values = {}
    checks = {}
    for name in ("base", "M1", "E2"):
        atom = MANIFEST["ions"][f"O{stage}"]["artifacts"][name]["path"]
        output = HERE / f"native_O{stage}_{name}.txt"
        subprocess.run([str(DRIVER), atom, str(grid), str(output)],
                       check=True, timeout=1800)
        arr = np.loadtxt(output)
        nlev = MANIFEST["ions"][f"O{stage}"]["artifacts"][name]["nlev"]
        if arr.shape != (1604, nlev + 3):
            raise ValueError(f"O {stage} {name} incomplete native profile: {arr.shape}")
        arr = arr.reshape(401, 4, nlev + 3)
        if not np.all(np.isfinite(arr)) or np.any(arr[:, :, 2] <= 0):
            raise ValueError(f"O {stage} {name} invalid cooling or populations")
        if not np.allclose(arr[:, :, 0], TEMPERATURE[:, None], rtol=1e-14):
            raise ValueError("Native temperature grid mismatch")
        if not np.allclose(arr[:, :, 1], DENSITY, rtol=0, atol=0):
            raise ValueError("Native density grid mismatch")
        pop = arr[:, :, 3:]
        norm = np.max(abs(pop.sum(axis=2) - 1))
        if np.min(pop) < -1e-13 or norm > 1e-12:
            raise ValueError(f"O {stage} {name} population normalization")
        values[name] = arr[:, :, 2]
        profiles[f"O{stage}_{name}_cooling"] = values[name]
        checks[name] = {"rows": 1604, "max_population_normalization_error": float(norm),
                        "minimum_population": float(pop.min()),
                        "cooling_min": float(values[name].min()),
                        "cooling_max": float(values[name].max())}
    for name, style in (("M1", "-"), ("E2", "--")):
        relative = 100 * (values[name] / values["base"] - 1)
        checks[name]["max_absolute_cooling_change_percent"] = float(np.max(abs(relative)))
        for i, density in enumerate(DENSITY):
            axis.semilogx(TEMPERATURE, relative[:, i], style,
                          label=f"{name}, n_e={density:.0e}")
    axis.set_ylabel(f"O {stage}: cooling change (percent)")
    axis.grid(True, alpha=.25)
    axis.legend(fontsize=7, ncol=2)
    results[f"O{stage}"] = checks
axes[-1].set_xlabel("Temperature (K)")
fig.tight_layout()
fig.savefig(HERE / "same_parity_native_sensitivity.pdf")
np.savez_compressed(HERE / "native_full_profiles.npz", **profiles)
report = {"passed": True, "temperature_count": 401, "density_count": 4,
          "model_sensitivity_only": True, "ions": results}
(HERE / "native_result.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
