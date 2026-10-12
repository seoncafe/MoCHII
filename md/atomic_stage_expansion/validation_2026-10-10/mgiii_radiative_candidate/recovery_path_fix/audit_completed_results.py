"""Verify preserved Mg III profiles after a controller artifact-name error."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
BASE = HERE.parent


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    state = json.loads((HERE / "controller/state.json").read_text())
    assert state["status"] == "failed"
    assert [(step["id"], step["status"]) for step in state["steps"]] == [
        ("pilot", "completed"),
        ("native", "completed"),
        ("direct", "completed"),
        ("compare", "failed"),
    ]
    assert state["steps"][-1]["return_code"] == 0
    plan = json.loads((HERE / "plan.json").read_text())
    assert plan["steps"][-1]["check"]["path"].endswith("compare_step.json")
    assert not (BASE / "compare_step.json").exists()
    assert (BASE / "comparison_step.json").exists()

    frozen = state["source_sha256"]
    changed = [name for name, expected in frozen.items()
               if digest(Path(name)) != expected]
    if changed:
        raise ValueError(f"Frozen input changed: {changed}")

    reports = {name: json.loads((BASE / f"{name}_step.json").read_text())
               for name in ("pilot", "native", "direct", "comparison")}
    assert all(report["passed"] for report in reports.values())
    native, direct, comparison = (reports[name] for name in
                                  ("native", "direct", "comparison"))
    products = {
        "full_atom": (BASE / "nlevel_mg_3_line_strength_candidate.txt",
                      native["full_atom_hash"]),
        "first27_atom": (BASE / "first27_original_A/nlevel_mg_3.txt",
                         native["first27_atom_hash"]),
        "driver": (BASE / "check_source_excitation_atom.x", native["driver_hash"]),
        "native_object": (Path(__file__).resolve().parents[5] / "src/nlevel_mod.o",
                          native["native_object_hash"]),
        "native_full": (BASE / "profile_native_full.txt", native["full_profile_hash"]),
        "native_first27": (BASE / "profile_native_first27.txt",
                           native["first27_profile_hash"]),
        "direct": (BASE / "direct_profiles.npz", direct["profile_hash"]),
    }
    for name, (path, expected) in products.items():
        if digest(path) != expected:
            raise ValueError(f"Product hash changed: {name}")

    with np.load(BASE / "direct_profiles.npz") as profile:
        assert profile["temperature_K"].shape == (801,)
        assert np.all(np.diff(profile["temperature_K"]) > 0)
        assert np.isclose(profile["temperature_K"][0], 1000)
        assert np.isclose(profile["temperature_K"][-1], 100000)
        assert np.array_equal(profile["density_cm3"], [1e-6, 1e2, 1e6, 1e10])
        maximum_cutoff = []
        for index in range(4):
            full = profile[f"cooling_full_{index}"]
            first27 = profile[f"cooling_first27_{index}"]
            assert full.shape == first27.shape == (801,)
            assert np.all(np.isfinite(full)) and np.all(full >= 0)
            assert np.all(np.isfinite(first27)) and np.all(first27 >= 0)
            for label, count in (("full", 167), ("first27", 27)):
                populations = profile[f"populations_{label}_{index}"]
                assert populations.shape == (801, count)
                assert np.all(np.isfinite(populations))
                assert np.min(populations) >= -1e-13
                assert np.max(abs(populations.sum(axis=1) - 1)) < 1e-12
            positive = full > 64 * np.finfo(float).tiny
            maximum_cutoff.append(float(np.max(abs(first27[positive] / full[positive] - 1))))
        assert np.allclose(maximum_cutoff,
                           comparison["first27_full_relative_cooling_errors"],
                           rtol=1e-12, atol=0)

    for label, count in (("full", 167), ("first27", 27)):
        rows = np.loadtxt(BASE / f"profile_native_{label}.txt")
        assert rows.shape == (3204, count + 3)
        assert np.all(np.isfinite(rows))
        assert np.all(rows[:, 2] >= 0)
    assert comparison["temperature_points"] == 801
    assert comparison["first27_within_0p3_percent"] is False
    for name in ("full_and_first27_cooling_profiles.pdf",
                 "first27_truncation_differences.pdf"):
        assert (BASE / name).stat().st_size > 1000

    result = {
        "passed": True,
        "original_controller_status": "failed",
        "original_error": "Artifact filename compare_step.json differs from comparison_step.json",
        "preserved_numerical_steps": ["pilot", "native", "direct", "comparison"],
        "frozen_input_count": len(frozen),
        "verified_product_hashes": {name: expected for name, (_, expected) in products.items()},
        "temperature_points": 801,
        "density_cm3": [1e-6, 1e2, 1e6, 1e10],
        "maximum_relative_native_cooling_errors":
            comparison["maximum_relative_native_cooling_errors"],
        "maximum_absolute_native_population_errors":
            comparison["maximum_absolute_native_population_errors"],
        "first27_full_relative_cooling_errors": maximum_cutoff,
        "first27_within_0p3_percent": False,
        "physical_source_accuracy_verified": False,
    }
    (HERE / "audit_completed_results.json").write_text(
        json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
