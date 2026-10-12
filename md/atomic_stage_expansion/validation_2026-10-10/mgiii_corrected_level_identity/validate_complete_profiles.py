"""Gate the corrected Mg III profile and balance products for fitting."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "mgiii_radiative_candidate"
RESULT = HERE / "validation_manifest.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name: str) -> dict:
    report = json.loads((HERE / name).read_text())
    if report.get("passed") is not True:
        raise ValueError(f"Validation result failed: {name}")
    return report


def finite_nonnegative(values, name: str) -> None:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)) or np.any(array < 0):
        raise ValueError(f"Invalid reported numerical result: {name}")


def maximum_condition(values: np.ndarray, temperatures: np.ndarray,
                      densities: np.ndarray) -> dict:
    index = np.unravel_index(int(np.argmax(values)), values.shape)
    return {"value": float(values[index]),
            "temperature_K": float(temperatures[index[0]]),
            "density_cm3": float(densities[index[1]])}


def main() -> None:
    if RESULT.exists():
        raise FileExistsError(f"Preserve existing manifest: {RESULT}")
    candidate = read("candidate_result.json")
    pilot = read("pilot_step.json")
    native = read("native_step.json")
    direct = read("direct_step.json")
    compare = read("comparison_step.json")
    ledger = read("independent_balance_result.json")
    atom = Path(candidate["candidate_path"])
    if digest(atom) != candidate["candidate_sha256"]:
        raise ValueError("Corrected candidate changed")
    if candidate["selected_level_167_cm1"] != 634430.0:
        raise ValueError("Corrected ID 167 energy changed")
    if candidate["source_collision_records_unchanged"] is not True:
        raise ValueError("Collision nodes or scaling changed")
    if native["full_atom_hash"] != candidate["candidate_sha256"]:
        raise ValueError("Native calculation used a different atom")
    if ledger["input_sha256"][str(atom)] != candidate["candidate_sha256"]:
        raise ValueError("Independent ledger used a different atom")
    if (pilot["full_levels"], native["full_levels"], direct["full_levels"],
            candidate["nlev"]) != (167, 167, 167, 167):
        raise ValueError("Full-atom level counts changed")
    if not all(report["temperature_points"] == 801 for report in
               (native, direct, compare)):
        raise ValueError("Profile temperature count changed")
    if ledger["temperature_count"] != 801 or ledger["density_cm3"] != [1e-6, 100.0, 1e6, 1e10]:
        raise ValueError("Ledger grid differs")
    errors = compare["maximum_relative_native_cooling_errors"]["full"]
    for label in ("full", "first27"):
        finite_nonnegative(compare["maximum_relative_native_cooling_errors"][label],
                           f"{label} cooling errors")
        finite_nonnegative(compare["maximum_absolute_native_population_errors"][label],
                           f"{label} population errors")
    finite_nonnegative(compare["first27_full_relative_cooling_errors"],
                       "first-27 truncation errors")
    finite_nonnegative(native["maximum_relative_source_strength_error"],
                       "native source strengths")
    finite_nonnegative(candidate["maximum_relative_serialized_vs_independent_A_error"],
                       "radiative serialization")
    finite_nonnegative(direct["minimum_cooling_erg_cm3_s"], "minimum cooling")
    for key in ("maximum_stationarity_residual_relative_to_total_gross_collision_flow",
                "maximum_local_stationarity_where_level_flow_exceeds_1e-12_total",
                "maximum_collision_minus_radiation_power_relative_to_gross_collision_power",
                "maximum_radiation_vs_native_cooling_relative_error"):
        finite_nonnegative(ledger[key]["value"], key)
    if len(errors) != 4 or max(errors) >= 1e-12:
        raise ValueError("Native and direct cooling profiles differ")
    if max(compare["maximum_absolute_native_population_errors"]["full"]) >= 1e-10:
        raise ValueError("Native and direct populations differ")
    if (ledger["maximum_stationarity_residual_relative_to_total_gross_collision_flow"]["value"] >= 1e-8
        or ledger["maximum_collision_minus_radiation_power_relative_to_gross_collision_power"]["value"] >= 1e-8
        or ledger["maximum_radiation_vs_native_cooling_relative_error"]["value"] >= 1e-8):
        raise ValueError("Statistical or energy balance failed")
    if (digest(HERE / "profile_native_full.txt") != native["full_profile_hash"]
        or digest(HERE / "direct_profiles.npz") != direct["profile_hash"]):
        raise ValueError("Saved population profile changed")
    old = OLD / "direct_profiles.npz"
    if digest(old) != "9e3786c39cca72aa09111ad8ad0727f9a366788afe7959afc002435601615722":
        raise ValueError("Preserved original candidate profile changed")
    with np.load(old) as previous, np.load(HERE / "direct_profiles.npz") as current:
        temperatures = current["temperature_K"]
        densities = current["density_cm3"]
        if not np.array_equal(previous["temperature_K"], temperatures):
            raise ValueError("Old and corrected profile grids differ")
        changes = np.stack([abs(current[f"cooling_full_{index}"] /
                                previous[f"cooling_full_{index}"] - 1)
                            for index in range(4)], axis=1)
        finite_nonnegative(changes, "old-to-corrected conditional difference")
        change = maximum_condition(changes, temperatures, densities)
        for index in range(4):
            if not np.array_equal(previous[f"cooling_first27_{index}"],
                                  current[f"cooling_first27_{index}"]):
                raise ValueError("Reused first-27 direct profile changed")
    result = {
        "passed": True,
        "candidate_path": str(atom),
        "candidate_sha256": candidate["candidate_sha256"],
        "level_167_identity_checked": True,
        "selected_level_167_cm1": 634430.0,
        "temperature_points": 801,
        "densities_cm3": [1e-6, 100.0, 1e6, 1e10],
        "maximum_relative_cooling_errors": errors,
        "maximum_absolute_population_errors":
            compare["maximum_absolute_native_population_errors"]["full"],
        "ledger_numerical_checks_passed": True,
        "maximum_stationarity_relative_to_gross_flow":
            ledger["maximum_stationarity_residual_relative_to_total_gross_collision_flow"],
        "maximum_energy_residual_relative_to_gross_collision_power":
            ledger["maximum_collision_minus_radiation_power_relative_to_gross_collision_power"],
        "first27_full_relative_cooling_errors":
            compare["first27_full_relative_cooling_errors"],
        "old_misassigned_candidate_to_corrected_maximum_relative_change": change,
        "old_candidate_physically_rejected": True,
        "physical_source_accuracy_verified": False,
        "source_hashes": candidate["raw_source_sha256"],
        "native_profile_sha256": native["full_profile_hash"],
        "direct_profile_sha256": direct["profile_hash"],
        "ledger_result_sha256": digest(HERE / "independent_balance_result.json"),
    }
    RESULT.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
