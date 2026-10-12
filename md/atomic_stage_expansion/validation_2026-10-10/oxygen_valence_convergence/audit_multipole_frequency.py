"""Quantify frequency sensitivity of same-parity oxygen radiative rates.

The cubic/quintic values are conditional M1/E2 sensitivities, not uncertainties
of the underlying atomic structure or proofs of the source A frequency.
"""

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools" / "fitting"))
from atomic_excitation_energies import selected_energies  # noqa: E402
from chianti_cooling import read_excitation_atom  # noqa: E402

OUT = Path(__file__).parent
result = {}
for stage, prefix in ((2, 5), (3, 6)):
    raw, radiation, collisions = read_excitation_atom("o", stage)
    hybrid = selected_energies("o", stage, raw, radiation, collisions)
    source_order = {i: (g, obs if i <= prefix and (obs > 0 or i == 1) else theory)
                    for i, (g, obs, theory) in raw.items()}
    path = Path.home() / f"RT_Codes/CHIANTI/dbase/o/o_{stage}/o_{stage}.elvlc"
    configs = {}
    for line in path.read_text().splitlines():
        if line.strip() == "-1":
            break
        configs[int(line[:7])] = line[7:37]

    def parity(index):
        parts = re.findall(r"(\d+)([spdfgh])(\d*)", configs[index])
        return sum((int(n) if n else 1) * "spdfgh".index(l) for _, l, n in parts) % 2

    evaluations = {}
    for name, energies in (("hybrid", hybrid), ("source_order", source_order)):
        rows = []
        for lower, upper, wavelength, gf, rate in radiation:
            if rate <= 0 or wavelength == 0 or parity(lower) != parity(upper):
                continue
            ratio = abs(energies[upper][1] - energies[lower][1]) * abs(wavelength) / 1e8
            rows.append({"lower": lower, "upper": upper, "wavelength_A": wavelength,
                         "source_A_s-1": rate, "gf": gf,
                         "selected_to_wavelength_frequency": ratio,
                         "M1_frequency_factor": ratio ** 3,
                         "E2_frequency_factor": ratio ** 5})
        evaluations[name] = {
            "same_parity_positive_rate_count": len(rows),
            "conditional_M1_change_over_10_percent_count": sum(abs(x["M1_frequency_factor"] - 1) > .1 for x in rows),
            "conditional_E2_change_over_10_percent_count": sum(abs(x["E2_frequency_factor"] - 1) > .1 for x in rows),
            "largest_conditional_M1_change_percent": max(100 * abs(x["M1_frequency_factor"] - 1) for x in rows),
            "largest_conditional_E2_change_percent": max(100 * abs(x["E2_frequency_factor"] - 1) for x in rows),
            "rows": rows,
        }
    result[f"O{stage}"] = evaluations
(OUT / "multipole_frequency.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({ion: {choice: {k: v for k, v in data.items() if k != "rows"}
                        for choice, data in options.items()}
                  for ion, options in result.items()}, indent=2))
