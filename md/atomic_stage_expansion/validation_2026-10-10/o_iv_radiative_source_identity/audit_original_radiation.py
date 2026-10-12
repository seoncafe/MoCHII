"""Compare O IV CHIANTI radiation with the published OPEN-ADAS ADF04 source.

The raw ADF04 file is retained outside the repository under /tmp. No raw
numerical source table is copied into this diagnostic directory.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C

ADF = Path("/tmp/blike_lgy12_o3_adf04.dat")
ADF_SHA = "f49d7a8b38d0b82ccff9f4f21871f4e8c785d875521c5be00b017cb9b61407fb"
PARENT_RESULT = ROOT / "md/atomic_stage_expansion/validation_2026-10-10/o_iv_main_radiative_review/result.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_float(value: str) -> float:
    return float(re.sub(r"(?<=\d)([+-]\d+)$", r"e\1", value))


def main() -> None:
    if sha256(ADF) != ADF_SHA:
        raise ValueError("OPEN-ADAS source hash changed")
    lines = ADF.read_text().splitlines()
    adf_energy = {}
    adf_weight = {}
    for line in lines[1:205]:
        match = re.match(r"^\s*(\d+)\s+.*\(\s*([\d.]+)\)\s+([\d.]+)\s*$", line)
        if match is None:
            raise ValueError(f"Unparsed ADF04 level: {line}")
        index, j, energy = match.groups()
        adf_energy[int(index)] = float(energy)
        adf_weight[int(index)] = 2 * float(j) + 1
    if sorted(adf_energy) != list(range(1, 205)):
        raise ValueError("ADF04 level numbering is incomplete")
    adf_radiation = {}
    for line in lines[207:]:
        match = re.match(r"^\s*(\d+)\s+(\d+)\s+([\d.]+[+-]\d+)\b", line)
        if match is None:
            continue
        upper, lower, rate = match.groups()
        key = (int(lower), int(upper))
        if key in adf_radiation:
            raise ValueError(f"Duplicate ADF04 radiative pair: {key}")
        adf_radiation[key] = compact_float(rate)

    stem = Path(C.ion_dir("o", 4)) / "o_4"
    raw_levels = C.read_elvlc(str(stem) + ".elvlc")
    selected_levels, selected_radiation, _ = C.bound_excitation_atom("o", 4)
    parent_result = json.loads(PARENT_RESULT.read_text())
    mismatches = parent_result["unmatched_higher_E1"]
    checked = []
    missing = []
    different = []
    physical_gap_differences = {"upper_at_most_20": [], "upper_above_20": []}
    for lower, upper, wavelength, gf, rate in selected_radiation:
        key = (lower, upper)
        if key not in adf_radiation:
            missing.append(key)
            continue
        relative_error = abs(rate / adf_radiation[key] - 1) if adf_radiation[key] else (0 if rate == 0 else math.inf)
        checked.append(relative_error)
        if relative_error > 0.001:
            different.append((lower, upper, relative_error))
        original_gap = adf_energy[upper] - adf_energy[lower]
        selected_gap = selected_levels[upper][1] - selected_levels[lower][1]
        if original_gap <= 0 or selected_gap <= 0:
            raise ValueError(f"Nonpositive radiative gap for {key}")
        gap_group = "upper_at_most_20" if upper <= 20 else "upper_above_20"
        physical_gap_differences[gap_group].append(abs(selected_gap / original_gap - 1))

    def relevant(pair):
        lower, upper = pair["lower"], pair["upper"]
        gap = adf_energy[upper] - adf_energy[lower]
        inferred = pair["inferred_Einstein_frequency_cm_inverse"]
        return {
            "pair": [lower, upper],
            "source_gf": pair["source_gf"],
            "chianti_A_s_inverse": pair["source_A_s_inverse"],
            "adas_A_matches_within_0p1_percent": abs(pair["source_A_s_inverse"] / adf_radiation[(lower, upper)] - 1) <= 0.001,
            "adas_gap_cm_inverse": gap,
            "einstein_frequency_cm_inverse": inferred,
            "einstein_to_adas_gap_ratio": inferred / gap if gap > 0 else None,
            "chianti_wavelength_angstrom": pair["wavelength_angstrom"],
        }

    level_differences = []
    weight_differences = []
    for index in range(1, 205):
        chianti_weight, _, theory = raw_levels[index]
        level_differences.append(abs(adf_energy[index] - theory))
        if chianti_weight != adf_weight[index]:
            weight_differences.append(index)
    result = {
        "passed": True,
        "physical_model_accepted": False,
        "adas_level_count": len(adf_energy),
        "adas_radiative_count": len(adf_radiation),
        "selected_radiative_count": len(selected_radiation),
        "selected_pairs_missing_in_adas": missing,
        "selected_A_difference_over_0p1_percent_count": len(different),
        "selected_A_difference_over_0p1_percent_with_lower_at_most_20": sum(lower <= 20 for lower, _, _ in different),
        "selected_A_difference_over_1_percent_with_lower_above_20": sum(lower > 20 and error > 0.01 for lower, _, error in different),
        "largest_A_disagreement_pair_ids": [[lower, upper] for lower, upper, _ in sorted(different, key=lambda row: row[2], reverse=True)[:8]],
        "maximum_selected_relative_A_difference": max(checked),
        "maximum_absolute_level_theory_difference_cm_inverse": max(level_differences),
        "level_weight_mismatches": weight_differences,
        "selected_physical_gap_vs_ADF04_theory": {
            group: {
                "pair_count": len(errors),
                "maximum_relative_difference": max(errors),
                "pairs_above_one_percent": sum(error > 0.01 for error in errors),
            }
            for group, errors in physical_gap_differences.items()
        },
        "parent_343_mismatch_count": len(mismatches),
        "parent_343_with_adas_A_match_within_0p1_percent": sum(
            (q := adf_radiation.get((r["lower"], r["upper"]))) is not None
            and abs(r["source_A_s_inverse"] / q - 1) <= 0.001
            for r in mismatches
        ),
        "selected_examples": [relevant(row) for row in sorted(mismatches, key=lambda row: row["minimum_available_frequency_relative_difference"], reverse=True)[:12]],
        "pair_77_87": relevant(next(row for row in parent_result["largest_conditional_wavelength_scaling_examples"] if (row["lower"], row["upper"]) == (77, 87))),
        "source_sha256": {
            str(ADF): sha256(ADF),
            str(PARENT_RESULT): sha256(PARENT_RESULT),
            **{str(Path(str(stem) + suffix)): sha256(Path(str(stem) + suffix)) for suffix in (".elvlc", ".wgfa", ".scups")},
        },
    }
    target = Path(__file__).with_name("result.json")
    if target.exists():
        raise FileExistsError("Preserve previous audit result")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print({key: value for key, value in result.items() if not isinstance(value, (list, dict))})


if __name__ == "__main__":
    main()
