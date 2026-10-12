"""Measure individual radiative-line sensitivity of private oxygen candidates."""

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
MANIFEST = json.loads((HERE / "candidate_manifest.json").read_text())
RESULT = {}
HC = 1.9864458571489287e-16


def radiative_rates(path):
    lines = Path(path).read_text().splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith("NLEV "))
    count = int(lines[index].split()[1])
    energies = {int(line.split()[0]): float(line.split()[2])
                for line in lines[index + 1:index + 1 + count]}
    index = next(i for i, line in enumerate(lines) if line.startswith("NRAD "))
    count = int(lines[index].split()[1])
    rates = {}
    for line in lines[index + 1:index + 1 + count]:
        lower, upper, value = line.split()
        rates[int(lower), int(upper)] = float(value)
    return energies, rates


for stage in (2, 3):
    files = MANIFEST["ions"][f"O{stage}"]["artifacts"]
    energies, base_rates = radiative_rates(files["base"]["path"])
    base = np.loadtxt(HERE / f"native_O{stage}_base.txt")
    temperature, density = base[:, 0], base[:, 1]
    base_cooling, base_pop = base[:, 2], base[:, 3:]
    measured = {}
    for name in ("M1", "E2"):
        current = np.loadtxt(HERE / f"native_O{stage}_{name}.txt")
        if not np.array_equal(current[:, :2], base[:, :2]):
            raise ValueError("Candidate line profiles use different grids")
        _, rates = radiative_rates(files[name]["path"])
        if set(rates) != set(base_rates):
            raise ValueError("Candidate radiative line identities differ")
        pop = current[:, 3:]
        strongest = {"relative_difference": 0.0}
        protected = 0.0
        included = 0
        for (lower, upper), source_a in base_rates.items():
            photon = (energies[upper] - energies[lower]) * HC
            base_power = base_pop[:, upper - 1] * source_a * photon / density
            meaningful = base_power > 1e-6 * base_cooling
            if not np.any(meaningful):
                continue
            included += 1
            ratio = pop[meaningful, upper - 1] * rates[lower, upper]
            ratio /= base_pop[meaningful, upper - 1] * source_a
            point = int(np.argmax(abs(ratio - 1)))
            difference = float(ratio[point] - 1)
            if abs(difference) > abs(strongest["relative_difference"]):
                original_index = np.flatnonzero(meaningful)[point]
                strongest = {"line_lower": lower, "line_upper": upper,
                             "relative_difference": difference,
                             "temperature_K": float(temperature[original_index]),
                             "density_cm-3": float(density[original_index]),
                             "base_line_fraction_of_cooling": float(base_power[original_index] / base_cooling[original_index])}
            if stage == 3 and (lower, upper) in ((2, 4), (3, 4)):
                protected = max(protected, float(np.max(abs(ratio - 1))))
        measured[name] = {"line_count_above_fraction_threshold": included,
                          "fraction_threshold": 1e-6,
                          "largest_line_change": strongest,
                          "maximum_OIII_ground_optical_line_change": protected if stage == 3 else None}
    RESULT[f"O{stage}"] = measured
(HERE / "line_sensitivity_result.json").write_text(json.dumps(RESULT, indent=2) + "\n")
print(json.dumps(RESULT, indent=2))
