#!/usr/bin/env python3
"""Audit oxygen radiative identities against source and selected energies."""

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C

BASE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect(element, stage):
    raw, radiation, collisions = C.read_excitation_atom(element, stage)
    selected, _, _ = C.bound_excitation_atom(element, stage)
    rows = []
    for lower, upper, wavelength, gf, rate in radiation:
        if rate <= 0 or lower not in selected or upper not in selected:
            continue
        selected_gap = selected[upper][1] - selected[lower][1]
        if selected_gap > 0:
            continue
        theoretical_gap = raw[upper][2] - raw[lower][2]
        observed_gap = raw[upper][1] - raw[lower][1]
        source_upper_g = raw[upper][0]
        reverse_upper_g = raw[lower][0]
        implied_upper_g = (gf / (1.49919e-16 * abs(wavelength)**2 * rate)
                           if gf > 0 and wavelength != 0 else None)
        implied_theory_upper_g = (gf / (1.49919e-16 * (1e8 / abs(theoretical_gap))**2 * rate)
                                  if gf > 0 and theoretical_gap != 0 else None)
        rows.append({"source_lower": lower, "source_upper": upper,
                     "source_wavelength_A": wavelength, "source_gf": gf,
                     "source_A_s1": rate,
                     "selected_gap_cm1": selected_gap,
                     "selected_lower_energy_cm1": selected[lower][1],
                     "selected_upper_energy_cm1": selected[upper][1],
                     "observed_gap_cm1": observed_gap,
                     "theoretical_gap_cm1": theoretical_gap,
                     "source_upper_g": source_upper_g,
                     "reverse_upper_g": reverse_upper_g,
                     "gf_A_implied_upper_g": implied_upper_g,
                     "gf_A_theoretical_gap_implied_upper_g": implied_theory_upper_g,
                     "gf_A_theoretical_original_g_rel_error": (
                         abs(implied_theory_upper_g/source_upper_g-1)
                         if implied_theory_upper_g is not None else None),
                     "gf_A_theoretical_reversed_g_rel_error": (
                         abs(implied_theory_upper_g/reverse_upper_g-1)
                         if implied_theory_upper_g is not None else None),
                     "original_upper_g_rel_error": (
                         abs(implied_upper_g/source_upper_g-1)
                         if implied_upper_g is not None else None),
                     "reversed_upper_g_rel_error": (
                         abs(implied_upper_g/reverse_upper_g-1)
                         if implied_upper_g is not None else None),
                     "source_wavelength_to_theory_ratio": (
                         abs(wavelength)*abs(theoretical_gap)/1e8
                         if theoretical_gap != 0 else None),
                     "source_wavelength_to_selected_ratio": (
                         abs(wavelength)*abs(selected_gap)/1e8
                         if selected_gap != 0 else None)})
    array = [row for row in rows if row["gf_A_implied_upper_g"] is not None]
    return {"bound_levels": len(selected), "positive_A_rows": sum(
                rate > 0 and lower in selected and upper in selected
                for lower, upper, wavelength, gf, rate in radiation),
            "reversed_rows": len(rows),
            "exact_zero_selected_gap_rows": sum(
                row["selected_gap_cm1"] == 0 for row in rows),
            "original_theoretical_gap_positive": sum(
                row["theoretical_gap_cm1"] > 0 for row in rows),
            "original_theoretical_gap_negative": sum(
                row["theoretical_gap_cm1"] < 0 for row in rows),
            "positive_gf_reversed_rows": len(array),
            "gf_A_original_upper_g_matches_1_percent": sum(
                row["original_upper_g_rel_error"] < .01 for row in array),
            "gf_A_reversed_upper_g_matches_1_percent": sum(
                row["reversed_upper_g_rel_error"] < .01 for row in array),
            "gf_A_theoretical_original_upper_g_matches_1_percent": sum(
                row["gf_A_theoretical_original_g_rel_error"] < .01
                for row in array),
            "gf_A_theoretical_reversed_upper_g_matches_1_percent": sum(
                row["gf_A_theoretical_reversed_g_rel_error"] < .01
                for row in array),
            "rows": rows}


def main():
    dependencies = [Path(__file__), ROOT / "tools/fitting/chianti_cooling.py",
                    ROOT / "tools/fitting/atomic_excitation_energies.py"]
    for stage in (2, 3):
        stem = Path(C.ion_dir("o", stage)) / f"o_{stage}"
        dependencies.extend(Path(f"{stem}.{extension}")
                            for extension in ("elvlc", "wgfa", "scups"))
    before = {str(path): digest(path) for path in dependencies}
    ions = {f"O{stage}": inspect("o", stage) for stage in (2, 3)}
    after = {str(path): digest(path) for path in dependencies}
    result = {"passed": before == after, "sources_unchanged": before == after,
              "source_sha256": before, "ions": ions}
    (BASE / "oxygen_radiation_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    if before != after:
        raise RuntimeError("Oxygen source changed during audit")
    for name, data in ions.items():
        print(name, data["reversed_rows"], data["positive_gf_reversed_rows"],
              data["gf_A_theoretical_original_upper_g_matches_1_percent"],
              data["gf_A_theoretical_reversed_upper_g_matches_1_percent"])


if __name__ == "__main__":
    main()
