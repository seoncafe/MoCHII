"""Gate complete Ne I source profiles and balance before private cooling fits."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RESULT = HERE / "validation_manifest.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name: str) -> dict:
    value = json.loads((HERE / name).read_text())
    if value.get("passed") is not True:
        raise ValueError(f"Ne I numerical predecessor failed: {name}")
    return value


def finite_nonnegative(values, name: str) -> None:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)) or np.any(array < 0):
        raise ValueError(f"Invalid Ne I numerical result: {name}")


def main() -> None:
    if RESULT.exists():
        raise FileExistsError(f"Preserve existing Ne I manifest: {RESULT}")
    candidate = read("candidate_result.json")
    pilot = read("pilot_step.json")
    native = read("native_step.json")
    direct = read("direct_step.json")
    compare = read("comparison_step.json")
    ledger = read("independent_balance_result.json")
    atom = Path(candidate["candidate_path"])
    if digest(atom) != candidate["candidate_sha256"]:
        raise ValueError("Completed Ne I atom changed")
    if not candidate["stout_added_radiative_pairs_checked"]:
        raise ValueError("Two Stout A values were not checked")
    if candidate["added_pairs_s-1"] != {"3,15": 9.00e5, "5,15": 6.82e7}:
        raise ValueError("Stout A identities changed")
    if (candidate["nlev"], candidate["nrad"], candidate["nups"],
            candidate["nodes"], candidate["bt_tmin_K"]) != (31, 144, 465, 6975, 1000.36):
        raise ValueError("Completed source atom counts or floor changed")
    if not (native["atom_sha256"] == compare["source_atom_sha256"] ==
            candidate["candidate_sha256"] == ledger["input_sha256"][str(atom)]):
        raise ValueError("Different Ne I atoms entered the numerical checks")
    if (native["temperature_points"], direct["temperature_points"],
            compare["temperature_points"], ledger["temperature_count"]) != (401, 401, 401, 401):
        raise ValueError("Ne I profile temperature grid changed")
    densities = [1e-6, 100.0, 1e6, 1e10]
    if any(report["densities_cm3"] != densities for report in
           (pilot, native, direct, compare)) or ledger["density_cm3"] != densities:
        raise ValueError("Ne I electron-density grid changed")
    if (pilot["radiative_rows"], native["nrad"], direct["nrad"],
            ledger["radiative_rows"]) != (144, 144, 144, 144):
        raise ValueError("One Ne I calculation omitted a restored decay")
    cooling_errors = compare["maximum_relative_cooling_errors"]
    population_errors = compare["maximum_absolute_population_errors"]
    strength_error = native["maximum_relative_strength_error"]
    if len(cooling_errors) != 4 or len(population_errors) != 4:
        raise ValueError("Incomplete density-dependent error report")
    finite_nonnegative(cooling_errors, "native cooling errors")
    finite_nonnegative(population_errors, "native population errors")
    finite_nonnegative(strength_error, "native source strength error")
    finite_nonnegative(direct["minimum_cooling_erg_cm3_s"], "minimum source cooling")
    if max(cooling_errors) >= 1e-12 or max(population_errors) >= 1e-10 or strength_error >= 1e-12:
        raise ValueError("Native and independent source profiles differ")
    for key in ("maximum_stationarity_residual_relative_to_total_gross_collision_flow",
                "maximum_local_stationarity_where_level_flow_exceeds_1e-12_total",
                "maximum_collision_minus_radiation_power_relative_to_gross_collision_power",
                "maximum_radiation_vs_native_cooling_relative_error"):
        finite_nonnegative(ledger[key]["value"], key)
    if (ledger["maximum_stationarity_residual_relative_to_total_gross_collision_flow"]["value"] >= 1e-8
        or ledger["maximum_collision_minus_radiation_power_relative_to_gross_collision_power"]["value"] >= 1e-8
        or ledger["maximum_radiation_vs_native_cooling_relative_error"]["value"] >= 1e-8):
        raise ValueError("Ne I level stationarity or energy balance failed")
    if (digest(HERE / "native_profiles.txt") != native["profile_sha256"]
        or digest(HERE / "direct_profiles.npz") != direct["profile_sha256"]):
        raise ValueError("Saved Ne I population profile changed")
    result = {
        "passed": True,
        "candidate_path": str(atom),
        "candidate_sha256": candidate["candidate_sha256"],
        "temperature_points": 401,
        "densities_cm3": densities,
        "maximum_relative_cooling_errors": cooling_errors,
        "maximum_absolute_population_errors": population_errors,
        "maximum_relative_strength_error": strength_error,
        "stout_added_radiative_pairs_checked": True,
        "ledger_numerical_checks_passed": True,
        "maximum_stationarity_relative_to_gross_flow":
            ledger["maximum_stationarity_residual_relative_to_total_gross_collision_flow"],
        "maximum_energy_residual_relative_to_gross_collision_power":
            ledger["maximum_collision_minus_radiation_power_relative_to_gross_collision_power"],
        "source_hashes": candidate["raw_source_sha256"],
        "native_profile_sha256": native["profile_sha256"],
        "direct_profile_sha256": direct["profile_sha256"],
        "ledger_result_sha256": digest(HERE / "independent_balance_result.json"),
        "physical_source_accuracy_verified": False,
        "previous_false_closed_atom_physically_rejected": True,
    }
    RESULT.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
