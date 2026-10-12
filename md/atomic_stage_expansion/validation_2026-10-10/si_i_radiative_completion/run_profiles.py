"""Check a private Si I atom against its independent physical source records."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C  # noqa: E402
from validate_excitation_density_profiles import full_excitation_balance  # noqa: E402


TEMPERATURE = np.geomspace(1000.0, 100000.0, 401)
DENSITY = [1e-6, 100.0, 1e6, 1e10]
DRIVER = HERE.parent / "oxygen_valence_convergence/optimized_build/source_check.x"
SOURCE_RESULT = (HERE.parent / "si_i_electric_dipole_source_audit/candidate_v2/"
                 "candidate_result.json")
EXPECTED_LEVELS = 35
EXPECTED_COLLISIONS = 595
EXPECTED_RADIATION = 250  # 229 existing pairs plus 21 new NIST multipole pairs.


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path: Path, payload: dict) -> None:
    if path.exists():
        raise FileExistsError(f"Preserve previous Si I result: {path}")
    path.write_text(json.dumps(payload, indent=2) + "\n")


def candidate():
    report = json.loads(SOURCE_RESULT.read_text())
    if report.get("passed") is not True:
        raise ValueError("Si I candidate has not passed its source checks")
    atom = Path(report["candidate_path"])
    module_path = Path(report["source_module_path"])
    if digest(atom) != report["candidate_sha256"]:
        raise ValueError("Si I candidate atom changed")
    if digest(module_path) != report["source_sha256"][str(module_path)]:
        raise ValueError("Si I source module changed")
    if (report["nlev"], report["nrad"], report["nups"]) != (
            EXPECTED_LEVELS, EXPECTED_RADIATION, EXPECTED_COLLISIONS):
        raise ValueError("Si I candidate is incomplete")
    for name, expected in report["source_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Si I raw source changed: {name}")
    spec = importlib.util.spec_from_file_location("private_si_i_source", module_path)
    if spec is None or spec.loader is None:
        raise ValueError("Si I source module cannot be loaded")
    sys.path.insert(0, str(module_path.parent))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    levels, radiation, collisions = module.source_atom()
    if (sorted(levels) != list(range(1, EXPECTED_LEVELS + 1))
            or len(radiation) != EXPECTED_RADIATION
            or len(collisions) != EXPECTED_COLLISIONS):
        raise ValueError("Independent Si I source has different counts")
    return report, atom, levels, radiation, collisions


def radiative_rows(radiation):
    rows = []
    seen = set()
    for entry in radiation:
        lower, upper = int(entry[0]), int(entry[1])
        rate = float(entry[4] if len(entry) >= 5 else entry[2])
        if (lower, upper) in seen or not np.isfinite(rate) or rate <= 0:
            raise ValueError("Duplicate or invalid Si I spontaneous transition")
        seen.add((lower, upper))
        rows.append((lower, upper, rate))
    return rows


def check_source(levels, radiation, collisions):
    energies = {key: level[1] for key, level in levels.items()}
    if any(not np.isfinite(level[0]) or level[0] <= 0 or
           not np.isfinite(level[1]) for level in levels.values()):
        raise ValueError("Invalid Si I source levels")
    for lower, upper, _ in radiative_rows(radiation):
        if energies[upper] <= energies[lower]:
            raise ValueError("Si I radiative transition has no positive photon gap")
    pairs = set()
    for row in collisions:
        first, second = int(row["ll"]), int(row["ul"])
        pair = tuple(sorted((first, second)))
        if first not in levels or second not in levels or first == second or pair in pairs:
            raise ValueError("Invalid or duplicate Si I electron transition")
        pairs.add(pair)
        strengths = C.upsilon(row, TEMPERATURE)
        if not np.all(np.isfinite(strengths)) or np.any(strengths < 0):
            raise ValueError("Negative or nonfinite Si I source strength")
    if len(pairs) != EXPECTED_COLLISIONS:
        raise ValueError("Si I bound collision pair set is incomplete")


def save_grid(path: Path, mode: int, temperature: np.ndarray) -> None:
    if path.exists():
        raise FileExistsError(f"Preserve native grid: {path}")
    path.write_text(f"{mode} {len(temperature)}\n" +
                    " ".join(f"{value:.17g}" for value in temperature) + "\n")


def check_balance(cooling, populations, ntemperature: int) -> None:
    if (cooling.shape != (ntemperature,) or
            populations.shape != (ntemperature, EXPECTED_LEVELS) or
            not np.all(np.isfinite(cooling)) or np.any(cooling < 0) or
            not np.all(np.isfinite(populations)) or
            np.min(populations) < -1e-13 or
            np.max(abs(populations.sum(axis=1) - 1)) > 1e-12):
        raise ValueError("Invalid Si I cooling or populations")


def pilot() -> None:
    report, atom, levels, radiation, collisions = candidate()
    check_source(levels, radiation, collisions)
    probe = np.asarray([1000.0, 10000.0, 100000.0])
    source = full_excitation_balance(levels, radiation, collisions, probe, DENSITY)
    for density in DENSITY:
        check_balance(*source[density], len(probe))
    save_json(HERE / "pilot_step.json", {
        "passed": True, "candidate_sha256": digest(atom),
        "temperatures_K": probe.tolist(), "densities_cm3": DENSITY,
        "nlev": len(levels), "nrad": len(radiation), "nups": len(collisions),
        "nodes": report["nodes"],
        "adopted_bt_tmin_K": report.get("bt_tmin_K"),
        "source_accuracy_verified": False,
    })


def native() -> None:
    started = time.monotonic()
    report, atom, levels, radiation, collisions = candidate()
    check_source(levels, radiation, collisions)
    probe = np.asarray([1000.0, 10000.0, 100000.0])
    strength_grid = HERE / "native_strength_grid.txt"
    profile_grid = HERE / "native_profile_grid.txt"
    strengths = HERE / "native_strengths.txt"
    profile = HERE / "native_profiles.txt"
    if strengths.exists() or profile.exists():
        raise FileExistsError("Preserve previous Si I native profiles")
    save_grid(strength_grid, 1, probe)
    save_grid(profile_grid, 2, TEMPERATURE)
    subprocess.run([str(DRIVER), str(atom), str(strength_grid), str(strengths)],
                   check=True, timeout=300)
    measured = np.loadtxt(strengths)
    collisions.sort(key=lambda row: tuple(sorted((row["ll"], row["ul"]),
                    key=lambda index: (levels[index][1], index))))
    expected = np.column_stack([C.upsilon(row, probe) for row in collisions])
    if (measured.shape != (len(probe), EXPECTED_COLLISIONS + 1) or
            not np.allclose(measured[:, 0], probe, rtol=1e-13, atol=0) or
            not np.all(np.isfinite(measured)) or
            not np.all(np.isfinite(expected)) or np.any(expected < 0)):
        raise ValueError("Incomplete or nonfinite Si I native strength profile")
    positive = expected > 0
    if np.any(measured[:, 1:][~positive] != 0):
        raise ValueError("Native Si I source zeros changed")
    strength_error = float(np.max(abs(measured[:, 1:][positive] /
                                      expected[positive] - 1)))
    if not np.isfinite(strength_error) or strength_error >= 1e-10:
        raise ValueError("Native Si I source strengths differ")
    subprocess.run([str(DRIVER), str(atom), str(profile_grid), str(profile)],
                   check=True, timeout=1200)
    rows = np.loadtxt(profile)
    if (rows.shape != (1604, EXPECTED_LEVELS + 3) or
            not np.allclose(rows[:, 0], np.repeat(TEMPERATURE, 4), rtol=1e-13, atol=0) or
            not np.array_equal(rows[:, 1], np.tile(DENSITY, 401))):
        raise ValueError("Incomplete Si I native population profile")
    check_balance(rows[:, 2].reshape(401, 4)[:, 0],
                  rows[:, 3:].reshape(401, 4, EXPECTED_LEVELS)[:, 0, :], 401)
    for index in range(1, 4):
        check_balance(rows[:, 2].reshape(401, 4)[:, index],
                      rows[:, 3:].reshape(401, 4, EXPECTED_LEVELS)[:, index, :], 401)
    save_json(HERE / "native_step.json", {
        "passed": True, "candidate_sha256": report["candidate_sha256"],
        "temperature_points": 401, "densities_cm3": DENSITY,
        "nlev": EXPECTED_LEVELS, "nrad": EXPECTED_RADIATION,
        "nups": EXPECTED_COLLISIONS, "nodes": report["nodes"],
        "maximum_relative_strength_error": strength_error,
        "exact_zero_strength_samples": int(np.sum(~positive)),
        "profile_sha256": digest(profile), "driver_sha256": digest(DRIVER),
        "native_object_sha256": digest(DRIVER.parent / "nlevel_mod.o"),
        "native_seconds": time.monotonic() - started,
    })


def direct() -> None:
    started = time.monotonic()
    report, atom, levels, radiation, collisions = candidate()
    check_source(levels, radiation, collisions)
    source = full_excitation_balance(levels, radiation, collisions,
                                     TEMPERATURE, DENSITY)
    payload = {"temperature_K": TEMPERATURE, "density_cm3": np.asarray(DENSITY)}
    for index, density in enumerate(DENSITY):
        cooling, populations = source[density]
        check_balance(cooling, populations, 401)
        payload[f"cooling_{index}"] = cooling
        payload[f"populations_{index}"] = populations
    output = HERE / "direct_profiles.npz"
    if output.exists():
        raise FileExistsError("Preserve previous Si I independent profiles")
    np.savez_compressed(output, **payload)
    save_json(HERE / "direct_step.json", {
        "passed": True, "candidate_sha256": report["candidate_sha256"],
        "temperature_points": 401, "densities_cm3": DENSITY,
        "nlev": EXPECTED_LEVELS, "nrad": EXPECTED_RADIATION,
        "nups": EXPECTED_COLLISIONS, "nodes": report["nodes"],
        "profile_sha256": digest(output),
        "independent_seconds": time.monotonic() - started,
    })


def compare() -> None:
    native_report = json.loads((HERE / "native_step.json").read_text())
    direct_report = json.loads((HERE / "direct_step.json").read_text())
    native_file = HERE / "native_profiles.txt"
    direct_file = HERE / "direct_profiles.npz"
    if (digest(native_file) != native_report["profile_sha256"] or
            digest(direct_file) != direct_report["profile_sha256"]):
        raise ValueError("Saved Si I profile changed")
    rows = np.loadtxt(native_file)
    if rows.shape != (1604, EXPECTED_LEVELS + 3):
        raise ValueError("Incomplete saved Si I native profile")
    cooling_errors, population_errors, unresolved = [], [], []
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True,
                             constrained_layout=True)
    difference, difference_axes = plt.subplots(2, 2, figsize=(9, 7),
                                               sharex=True, constrained_layout=True)
    with np.load(direct_file) as source:
        if (not np.array_equal(source["temperature_K"], TEMPERATURE) or
                not np.array_equal(source["density_cm3"], DENSITY)):
            raise ValueError("Independent Si I profile axes changed")
        for index, density in enumerate(DENSITY):
            cooling = source[f"cooling_{index}"]
            populations = source[f"populations_{index}"]
            measured = rows[index::4, 2]
            measured_pop = rows[index::4, 3:]
            check_balance(cooling, populations, 401)
            check_balance(measured, measured_pop, 401)
            resolved = cooling > 64 * np.finfo(float).tiny
            if np.any(measured[~resolved] != cooling[~resolved]):
                raise ValueError("Unresolved Si I cooling differs")
            if not np.any(resolved):
                raise ValueError("No resolved Si I cooling samples")
            cooling_error = float(np.max(abs(measured[resolved] /
                                             cooling[resolved] - 1)))
            population_error = float(np.max(abs(measured_pop - populations)))
            if not np.isfinite(cooling_error) or not np.isfinite(population_error):
                raise ValueError("Nonfinite Si I comparison error")
            cooling_errors.append(cooling_error)
            population_errors.append(population_error)
            unresolved.append(int(np.sum(~resolved)))
            panel = axes.flat[index]
            panel.loglog(TEMPERATURE, cooling, label="Corrected physical source")
            panel.loglog(TEMPERATURE, measured, "--", label="Compiled source atom")
            panel.set_title(r"$n_e=" + f"{density:g}" + r"\,\mathrm{cm}^{-3}$")
            panel.set_xlabel(r"$T$ [K]")
            panel.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
            panel.legend(fontsize=8)
            panel = difference_axes.flat[index]
            panel.semilogx(TEMPERATURE[resolved],
                           100 * (measured[resolved] / cooling[resolved] - 1))
            panel.set_title(r"$n_e=" + f"{density:g}" + r"\,\mathrm{cm}^{-3}$")
            panel.set_xlabel(r"$T$ [K]")
            panel.set_ylabel(r"Native / source - 1 [\%]")
    if max(cooling_errors) >= 1e-10 or max(population_errors) >= 1e-10:
        raise ValueError("Native and independent Si I profiles differ")
    for path in (HERE / "cooling_profiles.pdf", HERE / "native_differences.pdf"):
        if path.exists():
            raise FileExistsError(f"Preserve existing Si I figure: {path}")
    fig.savefig(HERE / "cooling_profiles.pdf")
    difference.savefig(HERE / "native_differences.pdf")
    plt.close(fig)
    plt.close(difference)
    save_json(HERE / "comparison_step.json", {
        "passed": True, "candidate_sha256": native_report["candidate_sha256"],
        "temperature_points": 401, "densities_cm3": DENSITY,
        "maximum_relative_cooling_errors": cooling_errors,
        "maximum_absolute_population_errors": population_errors,
        "unresolved_cooling_samples": unresolved,
        "native_profile_sha256": digest(native_file),
        "direct_profile_sha256": digest(direct_file),
        "cooling_figure_sha256": digest(HERE / "cooling_profiles.pdf"),
        "difference_figure_sha256": digest(HERE / "native_differences.pdf"),
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", required=True,
                        choices=("pilot", "native", "direct", "compare"))
    args = parser.parse_args()
    globals()[args.step]()


if __name__ == "__main__":
    main()
