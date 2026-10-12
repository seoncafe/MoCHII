"""Audit oxygen radiative rates under the selected observed/theoretical energies.

This diagnostic does not alter the production atomic data. CHIANTI gf is the
weighted oscillator strength for allowed transitions; its Einstein identity
recovers the frequency at which the source A and gf were calculated.
"""

import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools" / "fitting"))
from atomic_excitation_energies import selected_energies  # noqa: E402
from chianti_cooling import read_excitation_atom  # noqa: E402

EINSTEIN = 6.6703e15  # A_ki = EINSTEIN gf/(g_k lambda_A**2)
OUT = Path(__file__).parent
results = {}
for stage in (2, 3):
    levels, rows, collisions = read_excitation_atom("o", stage)
    selected = selected_energies("o", stage, levels, rows, collisions)
    source = Path.home() / f"RT_Codes/CHIANTI/dbase/o/o_{stage}/o_{stage}.elvlc"
    configurations = {}
    for line in source.read_text().splitlines():
        if line.strip() == "-1":
            break
        configurations[int(line[:7])] = line[7:37]

    def parity(index):
        occupations = re.findall(r"(\d+)([spdfgh])(\d*)", configurations[index])
        if not occupations:
            raise ValueError(f"Missing configuration for O {stage} level {index}")
        return sum((int(count) if count else 1) * "spdfgh".index(orbital)
                   for _, orbital, count in occupations) % 2

    records = []
    for lower, upper, wavelength, gf, rate in rows:
        if rate <= 0:
            continue
        delta = selected[upper][1] - selected[lower][1]
        reversed_order = delta < 0
        opposite_parity = parity(lower) != parity(upper)
        if gf > 0:
            source_frequency = 1e8 * np.sqrt(rate * levels[upper][0] / (EINSTEIN * gf))
            theory_frequency = levels[upper][2] - levels[lower][2]
            wavelength_frequency = 1e8 / abs(wavelength)
            theoretical_error = abs(source_frequency / theory_frequency - 1) if theory_frequency > 0 else None
            wavelength_error = abs(source_frequency / wavelength_frequency - 1)
            if theoretical_error is not None and theoretical_error < .02:
                source_frequency_kind = "theoretical"
            elif wavelength_error < .02:
                source_frequency_kind = "wavelength"
            else:
                source_frequency_kind = "unidentified"
            if opposite_parity:
                physical_upper = upper if delta > 0 else lower
                factor = (abs(delta) / source_frequency) ** 3 * levels[upper][0] / levels[physical_upper][0]
                corrected = float(rate * factor)
            else:
                corrected = None
                factor = None
        else:
            source_frequency = None
            theoretical_error = None
            wavelength_error = None
            source_frequency_kind = "non_E1"
            corrected = None
            factor = None
        records.append({"source_lower": lower, "source_upper": upper,
                        "wavelength_A": wavelength, "gf": gf, "source_A_s-1": rate,
                        "selected_gap_cm-1": delta, "reversed": reversed_order,
                        "opposite_parity": opposite_parity,
                        "source_frequency_cm-1": float(source_frequency) if source_frequency is not None else None,
                        "frequency_matches": source_frequency_kind,
                        "theory_relative_error": float(theoretical_error) if theoretical_error is not None else None,
                        "wavelength_relative_error": float(wavelength_error) if wavelength_error is not None else None,
                        "conditional_E1_A_s-1": corrected,
                        "conditional_E1_factor": float(factor) if factor is not None else None})
    e1 = [r for r in records if r["gf"] > 0 and r["opposite_parity"]]
    reversed_e1 = [r for r in e1 if r["reversed"]]
    non_e1_reversed = [r for r in records if r["reversed"] and r["gf"] == 0]
    results[f"O{stage}"] = {
        "source_level_count": len(levels), "source_radiative_count": len(rows),
        "positive_E1_count": len(e1), "reversed_E1_count": len(reversed_e1),
        "reversed_non_E1": non_e1_reversed,
        "unidentified_E1_count": sum(r["frequency_matches"] == "unidentified" for r in e1),
        "reversed_E1_source_frequency_max_relative_theory_error": max(r["theory_relative_error"] for r in reversed_e1),
        "reversed_E1_factor_min": min(r["conditional_E1_factor"] for r in reversed_e1),
        "reversed_E1_factor_max": max(r["conditional_E1_factor"] for r in reversed_e1),
        "E1_frequency_counts": {kind: sum(r["frequency_matches"] == kind for r in e1)
                                for kind in ("theoretical", "wavelength", "unidentified")},
        "rows": records,
    }
(OUT / "e1_reciprocity.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps({ion: {k: v for k, v in report.items() if k != "rows"}
                  for ion, report in results.items()}, indent=2))
