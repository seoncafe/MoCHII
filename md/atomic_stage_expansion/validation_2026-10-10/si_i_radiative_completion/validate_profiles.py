"""Gate the Si I numerical source comparison before product fitting."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from run_profiles import (DENSITY, EXPECTED_COLLISIONS, EXPECTED_LEVELS,
                          EXPECTED_RADIATION, HERE, candidate, digest, save_json)


def passed(name: str) -> dict:
    value = json.loads((HERE / name).read_text())
    if value.get("passed") is not True:
        raise ValueError(f"Si I numerical predecessor failed: {name}")
    return value


def finite_below(values, bound: float, label: str) -> None:
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)) or np.any(array < 0) or np.any(array >= bound):
        raise ValueError(f"Si I {label} exceeds its numerical tolerance")


def main() -> None:
    candidate_report, atom, _, _, _ = candidate()
    pilot = passed("pilot_step.json")
    native = passed("native_step.json")
    direct = passed("direct_step.json")
    comparison = passed("comparison_step.json")
    ledger = passed("energy_balance_step.json")
    hash_value = candidate_report["candidate_sha256"]
    if any(value["candidate_sha256"] != hash_value
           for value in (pilot, native, direct, comparison, ledger)):
        raise ValueError("Different Si I atoms entered the numerical checks")
    if (pilot["nlev"], pilot["nrad"], pilot["nups"]) != (
            EXPECTED_LEVELS, EXPECTED_RADIATION, EXPECTED_COLLISIONS):
        raise ValueError("Si I pilot source counts changed")
    for value in (native, direct):
        if (value["nlev"], value["nrad"], value["nups"], value["nodes"]) != (
                EXPECTED_LEVELS, EXPECTED_RADIATION, EXPECTED_COLLISIONS,
                candidate_report["nodes"]):
            raise ValueError("Si I full source counts changed")
    for value in (native, direct, comparison, ledger):
        if value["temperature_points"] != 401 or value["densities_cm3"] != DENSITY:
            raise ValueError("Incomplete Si I temperature-density grid")
    if (digest(HERE / "native_profiles.txt") != native["profile_sha256"] or
            digest(HERE / "direct_profiles.npz") != direct["profile_sha256"] or
            ledger["native_profile_sha256"] != native["profile_sha256"] or
            comparison["native_profile_sha256"] != native["profile_sha256"] or
            comparison["direct_profile_sha256"] != direct["profile_sha256"]):
        raise ValueError("Saved Si I profile changed")
    for name, recorded in (("cooling_profiles.pdf", comparison["cooling_figure_sha256"]),
                           ("native_differences.pdf", comparison["difference_figure_sha256"])):
        if digest(HERE / name) != recorded:
            raise ValueError(f"Si I figure changed: {name}")
    finite_below(native["maximum_relative_strength_error"], 1e-10,
                 "source strength error")
    finite_below(comparison["maximum_relative_cooling_errors"], 1e-10,
                 "cooling profile errors")
    finite_below(comparison["maximum_absolute_population_errors"], 1e-10,
                 "population profile errors")
    if len(comparison["maximum_relative_cooling_errors"]) != 4 or len(
            comparison["maximum_absolute_population_errors"]) != 4:
        raise ValueError("Incomplete Si I density comparison")
    for key in ("maximum_global_stationarity", "maximum_significant_local_stationarity",
                "maximum_energy_error", "maximum_cooling_power_error"):
        finite_below(ledger[key]["value"], 1e-8, key)
    save_json(HERE / "validation_manifest.json", {
        "passed": True, "candidate_path": str(atom),
        "candidate_sha256": hash_value, "temperature_points": 401,
        "densities_cm3": DENSITY,
        "nlev": EXPECTED_LEVELS, "nrad": EXPECTED_RADIATION,
        "nups": EXPECTED_COLLISIONS, "nodes": candidate_report["nodes"],
        "bt_tmin_K": candidate_report.get("bt_tmin_K"),
        "maximum_relative_cooling_errors":
            comparison["maximum_relative_cooling_errors"],
        "maximum_absolute_population_errors":
            comparison["maximum_absolute_population_errors"],
        "maximum_relative_strength_error": native["maximum_relative_strength_error"],
        "ledger_numerical_checks_passed": True,
        "maximum_global_stationarity": ledger["maximum_global_stationarity"],
        "maximum_energy_error": ledger["maximum_energy_error"],
        "source_hashes": candidate_report["source_sha256"],
        "native_profile_sha256": native["profile_sha256"],
        "direct_profile_sha256": direct["profile_sha256"],
        "source_accuracy_verified": False,
    })


if __name__ == "__main__":
    main()
