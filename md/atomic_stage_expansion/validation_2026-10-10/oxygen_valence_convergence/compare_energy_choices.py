"""Compare oxygen ground-excitation power for two selected energy sets."""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools" / "fitting"))
from atomic_excitation_energies import selected_energies  # noqa: E402
from chianti_cooling import (  # noqa: E402
    COLL_PREF, HC_ERG_CM, K_B_ERG, bound_excitation_atom,
    excitation_levels, read_excitation_atom, upsilon,
)

TEMPERATURE = np.geomspace(1000.0, 100000.0, 801)
OUT = Path(__file__).parent
profiles = {"temperature_K": TEMPERATURE}
summary = {}
fig, ax = plt.subplots(figsize=(7, 3.8))
for stage, observed_prefix in ((2, 5), (3, 6)):
    raw, radiation, collisions = read_excitation_atom("o", stage)
    bound, _, retained = bound_excitation_atom("o", stage)
    hybrid = selected_energies("o", stage, raw, radiation, collisions)
    source_order = {i: (g, observed if i <= observed_prefix and (observed > 0 or i == 1)
                        else theoretical)
                    for i, (g, observed, theoretical) in raw.items()}
    values = {}
    for name, energy in (("hybrid", hybrid), ("source_order", source_order)):
        power = np.zeros_like(TEMPERATURE)
        for transition in retained:
            lower, upper = excitation_levels(energy, transition)
            if lower != 1 or upper not in bound:
                continue
            gap = abs(energy[upper][1] - energy[1][1]) * HC_ERG_CM
            power += (COLL_PREF / (energy[1][0] * np.sqrt(TEMPERATURE))
                      * upsilon(transition, TEMPERATURE)
                      * np.exp(-gap / (K_B_ERG * TEMPERATURE)) * gap)
        values[name] = power
        profiles[f"O{stage}_{name}"] = power
    ratio = values["source_order"] / values["hybrid"]
    summary[f"O{stage}"] = {
        "maximum_absolute_relative_difference_percent": float(100 * np.max(abs(ratio - 1))),
        "difference_at_10000K_percent": float(100 * (ratio[400] - 1)),
        "difference_at_100000K_percent": float(100 * (ratio[-1] - 1)),
    }
    ax.semilogx(TEMPERATURE, 100 * (ratio - 1), label=f"O {stage}")
ax.axhline(0, color="black", lw=.6)
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Coronal cooling change (percent)")
ax.legend()
ax.grid(True, alpha=.25)
fig.tight_layout()
fig.savefig(OUT / "energy_choice_ground_excitation.pdf")
np.savez_compressed(OUT / "energy_choice_ground_excitation.npz", **profiles)
(OUT / "energy_choice_result.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
