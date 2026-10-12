"""Identify electric-dipole reversals in the bound Mg III source; write no atom."""
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
from chianti_cooling import bound_excitation_atom, read_excitation_atom


def main():
    source = Path.home() / "RT_Codes/CHIANTI/dbase/mg/mg_3"
    original, _, _ = read_excitation_atom("mg", 3)
    selected, radiation, collisions = bound_excitation_atom("mg", 3)
    configurations = {}
    for line in (source / "mg_3.elvlc").read_text().splitlines():
        if line.strip() == "-1":
            break
        configurations[int(line[:7])] = line[7:37]

    def parity(index):
        configuration = re.sub(r"\([^)]*\)", "", configurations[index])
        occupations = re.findall(r"(\d+)([spdfgh])(\d*)", configuration)
        if not occupations:
            raise ValueError(f"Missing configuration for level {index}")
        return sum((int(count) if count else 1) * "spdfgh".index(orbital)
                   for _, orbital, count in occupations) % 2

    records = []
    for lower, upper, wavelength, gf, rate in radiation:
        if rate <= 0:
            continue
        gap = selected[upper][1] - selected[lower][1]
        opposite = parity(lower) != parity(upper)
        jl, ju = ((original[i][0] - 1) / 2 for i in (lower, upper))
        dipole_selection = opposite and abs(ju - jl) <= 1 and ju + jl > 0
        inferred = (1e8 * np.sqrt(rate * original[upper][0] / (6.6703e15 * gf))
                    if gf > 0 else None)
        theory = original[upper][2] - original[lower][2]
        records.append({
            "source_lower": lower, "source_upper": upper,
            "physical_gap_cm-1": gap, "reversed": gap < 0,
            "source_A_s-1": rate, "source_gf": gf,
            "source_wavelength_A": wavelength,
            "opposite_parity": opposite, "dipole_selection": dipole_selection,
            "A_gf_frequency_cm-1": inferred,
            "theoretical_gap_cm-1": theory,
            "frequency_relative_theory_difference":
                abs(inferred / theory - 1) if inferred is not None and theory > 0 else None,
            "frequency_relative_wavelength_difference":
                abs(inferred * abs(wavelength) / 1e8 - 1) if inferred is not None else None,
        })
    reversed_rows = [row for row in records if row["reversed"]]
    result = {
        "scope": "Transition identity and source-frequency diagnostic; no accepted atom",
        "source_sha256": {
            str(source / f"mg_3.{suffix}"):
                hashlib.sha256((source / f"mg_3.{suffix}").read_bytes()).hexdigest()
            for suffix in ("elvlc", "wgfa", "scups")
        },
        "bound_levels": len(selected), "radiative_rows": len(radiation),
        "collision_rows": len(collisions),
        "reversed_positive_A_rows": len(reversed_rows),
        "all_reversed_rows_obey_E1_selection_and_have_positive_gf":
            all(row["dipole_selection"] and row["source_gf"] > 0 for row in reversed_rows),
        "maximum_reversed_frequency_relative_theory_difference":
            max(row["frequency_relative_theory_difference"] for row in reversed_rows),
        "rows": records,
    }
    output = Path(__file__).with_name("result.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in {"rows", "source_sha256"}}, indent=2))


if __name__ == "__main__":
    main()
