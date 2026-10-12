"""Assess radiative ordering and ground-excitation cooling of O II/O III."""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

repo = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(repo / "tools" / "fitting"))
from atomic_excitation_energies import selected_energies  # noqa: E402
from chianti_cooling import (  # noqa: E402
    COLL_PREF, HC_ERG_CM, K_B_ERG, excitation_levels, read_excitation_atom,
    upsilon,
)

temperature = np.geomspace(1000.0, 100000.0, 801)
prefixes = {2: (5, 16, 31, 32, 64, 196), 3: (6, 16, 32, 64, 89, 90, 177)}
result = {}
profiles = {"temperature_K": temperature}


def ground_paths(ids, radiation, energy):
    paths = {i: set() for i in ids}
    for lower, upper, wavelength, gf, rate in radiation:
        if lower in ids and upper in ids and rate > 0 and energy[upper][1] > energy[lower][1]:
            paths[upper].add(lower)

    def reaches_ground(level, visited):
        if level == 1:
            return True
        if level in visited:
            return False
        return any(reaches_ground(lower, visited | {level}) for lower in paths[level])

    return [i for i in ids if not reaches_ground(i, set())]


for stage in (2, 3):
    raw, radiation, collisions = read_excitation_atom("o", stage)
    energy = selected_energies("o", stage, raw, radiation, collisions)
    selected = {i for i in energy if i <= max(prefixes[stage])}
    values = {}
    for tr in collisions:
        lower, upper = excitation_levels(energy, tr)
        if lower != 1 or upper not in selected:
            continue
        gap = (energy[upper][1] - energy[1][1]) * HC_ERG_CM
        strength = upsilon(tr, temperature)
        if np.any(strength < -1e-14):
            raise ValueError(f"Negative effective collision strength for O {stage} {lower}-{upper}")
        power = COLL_PREF / (energy[1][0] * np.sqrt(temperature))
        power *= strength * np.exp(-gap / (K_B_ERG * temperature)) * gap
        values[upper] = values.get(upper, 0) + power
    totals = {}
    inspection = {}
    for prefix in prefixes[stage]:
        ids = {i for i in selected if i <= prefix}
        reversed_rows = [(lower, upper, rate) for lower, upper, wavelength, gf, rate in radiation
                         if lower in ids and upper in ids and rate > 0
                         and energy[upper][1] <= energy[lower][1]]
        no_ground_path = ground_paths(ids, radiation, energy)
        total = sum((power for i, power in values.items() if i <= prefix), np.zeros_like(temperature))
        totals[prefix] = total
        profiles[f"O{stage}_{prefix}"] = total
        inspection[str(prefix)] = {
            "reversed_A_count": len(reversed_rows),
            "first_reversed_A": reversed_rows[:3],
            "no_valid_radiative_path_to_ground": no_ground_path,
            "ground_excitation_collision_count": sum(i <= prefix for i in values),
        }
    full = totals[max(prefixes[stage])]
    for prefix, total in totals.items():
        difference = np.abs(total / full - 1.0)
        inspection[str(prefix)]["coronal_difference_from_full_max_percent"] = float(100 * difference.max())
        inspection[str(prefix)]["coronal_difference_from_full_at_10000K_percent"] = float(
            100 * difference[400]
        )
    result[f"O{stage}"] = inspection

fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
for axis, stage in zip(axes, (2, 3)):
    full = profiles[f"O{stage}_{max(prefixes[stage])}"]
    for prefix in prefixes[stage][:-1]:
        ratio = profiles[f"O{stage}_{prefix}"] / full
        axis.semilogx(temperature, 100 * (1 - ratio), label=f"{prefix} levels")
    axis.set_ylabel(f"O {stage}: omitted ground excitation (percent)")
    axis.legend(fontsize=8)
    axis.grid(True, alpha=0.25)
axes[-1].set_xlabel("Temperature (K)")
fig.tight_layout()
output = Path(__file__).parent
fig.savefig(output / "ground_excitation_prefixes.pdf")
np.savez_compressed(output / "ground_excitation_prefixes.npz", **profiles)
(output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({ion: {n: x["coronal_difference_from_full_max_percent"]
                        for n, x in report.items()} for ion, report in result.items()}, indent=2))
