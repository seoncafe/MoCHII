"""Audit physical transition gaps under CHIANTI level-energy choices.

This diagnostic does not alter CHIANTI files or silently discard transitions.
CHIANTI indices identify levels; they do not by themselves establish the
direction of a radiative or collisional transition after energies are chosen.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from chianti_cooling import ion_dir, read_elvlc, read_scups, read_wgfa

RY_WAVENUMBER = 109737.31568539  # cm^-1, CODATA 2022 Rydberg constant


def select_level_energies(levels, policy):
    """Return level wavenumbers without changing CHIANTI level identities.

    A negative entry is missing. The zero energy of level 1 is the ground
    state. Other zero entries are retained as supplied and will fail the audit
    if they produce a nonphysical transition gap.
    """
    if policy not in ("observed", "theoretical", "mixed"):
        raise ValueError(f"unknown energy policy: {policy}")
    result = {}
    for index, (_, observed, theoretical) in levels.items():
        observed_valid = observed >= 0 and (index == 1 or observed > 0)
        theoretical_valid = theoretical >= 0 and (index == 1 or theoretical > 0)
        if policy == "observed":
            energy = observed if observed_valid else None
        elif policy == "theoretical":
            energy = theoretical if theoretical_valid else None
        else:
            energy = observed if observed_valid else (
                theoretical if theoretical_valid else None
            )
        result[index] = energy
    return result


def audit_transitions(levels, radiation, collisions, policy):
    """Classify gaps in the file's directed level pairs for one energy policy."""
    energies = select_level_energies(levels, policy)
    summary = {
        "policy": policy,
        "levels": len(levels),
        "missing_level_energies": sum(value is None for value in energies.values()),
        "radiative_positive": 0,
        "radiative_reversed": [],
        "radiative_equal": [],
        "radiative_missing": [],
        "collision_positive": 0,
        "collision_reversed": [],
        "collision_equal": [],
        "collision_missing": [],
    }
    for lower, upper, wavelength, _, rate in radiation:
        if rate <= 0:
            continue
        pair = [lower, upper]
        if lower not in energies or upper not in energies or energies[lower] is None or energies[upper] is None:
            summary["radiative_missing"].append(pair)
            continue
        gap = energies[upper] - energies[lower]
        if gap > 0:
            summary["radiative_positive"] += 1
        else:
            record = {"pair": pair, "gap_cm1": gap,
                      "wavelength_A": wavelength, "rate_s1": rate}
            summary["radiative_equal" if gap == 0 else "radiative_reversed"].append(record)
    for transition in collisions:
        lower, upper = transition["ll"], transition["ul"]
        pair = [lower, upper]
        if lower not in energies or upper not in energies or energies[lower] is None or energies[upper] is None:
            summary["collision_missing"].append(pair)
            continue
        gap = energies[upper] - energies[lower]
        if gap > 0:
            summary["collision_positive"] += 1
        else:
            record = {"pair": pair, "gap_cm1": gap,
                      "scups_scaling_gap_cm1": transition["de"] * RY_WAVENUMBER}
            summary["collision_equal" if gap == 0 else "collision_reversed"].append(record)
    return summary


def audit_ion(element, stage):
    """Audit CHIANTI raw files for one ion under all three energy choices."""
    stem = Path(ion_dir(element, stage)) / f"{element}_{stage}"
    levels = read_elvlc(str(stem) + ".elvlc")
    radiation = read_wgfa(str(stem) + ".wgfa")
    collisions = read_scups(str(stem) + ".scups")
    return {
        "element": element,
        "stage": stage,
        "files": {extension: str(stem) + extension for extension in
                  (".elvlc", ".wgfa", ".scups")},
        "policies": {policy: audit_transitions(levels, radiation, collisions, policy)
                     for policy in ("observed", "theoretical", "mixed")},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ion", nargs="+", help="CHIANTI ion name, such as n_4")
    parser.add_argument("--output", type=Path, help="JSON output file; default stdout")
    args = parser.parse_args()
    records = []
    for name in args.ion:
        element, separator, stage = name.rpartition("_")
        if not separator or not element or not stage.isdecimal():
            parser.error(f"invalid ion name: {name}")
        records.append(audit_ion(element, int(stage)))
    output = json.dumps({"ions": records}, indent=2) + "\n"
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
