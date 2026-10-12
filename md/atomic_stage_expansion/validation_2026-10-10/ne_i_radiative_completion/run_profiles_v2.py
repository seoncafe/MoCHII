"""Validate the completed Ne I radiation model at four electron densities."""

from __future__ import annotations

import argparse
import hashlib
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
from prepare_candidate import OLD, STOUT, records, source_rates  # noqa: E402


SOURCE = Path.home() / "RT_Codes/CHIANTI/dbase/ne/ne_1/ne_1"
ATOM = HERE / "atom/nlevel_ne_1.txt"
DRIVER = HERE.parent / "oxygen_valence_convergence/optimized_build/source_check.x"
TEMPERATURE = np.geomspace(1000.0, 100000.0, 401)
DENSITY = [1e-6, 100.0, 1e6, 1e10]
FLOOR = 1000.36


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path: Path, value: dict) -> None:
    if path.exists():
        raise FileExistsError(f"Preserve existing result: {path}")
    path.write_text(json.dumps(value, indent=2) + "\n")


def stout_level_energy_difference(levels):
    """Check matching source indices and weights with rounded Stout energies."""
    stout_levels = {}
    for line in Path(str(STOUT) + ".nrg").read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0].isdigit() and int(fields[0]) <= 31:
            stout_levels[int(fields[0])] = (float(fields[2]), float(fields[1]))
    if sorted(stout_levels) != list(range(1, 32)):
        raise ValueError("Stout first-31 level count or indices changed")
    if any(stout_levels[index][0] != levels[index][0] for index in levels):
        raise ValueError("Stout and CHIANTI statistical weights differ")
    differences = {index: abs(stout_levels[index][1] - levels[index][1])
                   for index in levels}
    if max(differences.values()) > 0.002:
        raise ValueError("Stout and CHIANTI first-31 energies differ beyond rounding")
    for index in (3, 5, 15):
        if differences[index] != 0:
            raise ValueError("A level supplying the added Stout rates differs")
    return differences


def source_atom():
    levels, old_rates, _, _ = records(OLD)
    raw_levels = C.read_elvlc(str(SOURCE) + ".elvlc")
    if len(levels) != 31 or any(raw_levels[index][0] != levels[index][0]
                                for index in levels):
        raise ValueError("CHIANTI level identities changed")
    stout_level_energy_difference(levels)
    raw = C.read_wgfa(str(SOURCE) + ".wgfa")
    radiation = []
    source_rates_by_pair = {}
    for lower, upper, wavelength, gf, rate in raw:
        if lower not in levels or upper not in levels or rate <= 0:
            continue
        if levels[upper][1] <= levels[lower][1]:
            raise ValueError("Original CHIANTI radiative row has nonpositive photon energy")
        pair = (lower, upper)
        if pair in source_rates_by_pair:
            raise ValueError("Original CHIANTI radiative pair is duplicated")
        source_rates_by_pair[pair] = rate
        radiation.append((lower, upper, wavelength, gf, rate))
    if source_rates_by_pair != old_rates or len(radiation) != 142:
        raise ValueError("Historical CHIANTI radiation differs from the raw source")
    stout = source_rates()
    for lower, upper, expected in ((3, 15, 9.00e5), (5, 15, 6.82e7)):
        if (lower, upper) in source_rates_by_pair or stout.get((lower, upper)) != expected:
            raise ValueError("Added Stout decay is missing or duplicated")
        radiation.append((lower, upper, 0.0, 0.0, expected))
    collisions = C.read_scups(str(SOURCE) + ".scups")
    if len(collisions) != 465:
        raise ValueError("Raw Ne I collision-row count changed")
    for transition in collisions:
        transition["bt_tmin_K"] = FLOOR
        if transition["ll"] not in levels or transition["ul"] not in levels:
            raise ValueError("Raw collision row is outside the 31-level atom")
    return levels, radiation, collisions


def grid_file(path: Path, mode: int, temperatures: np.ndarray) -> None:
    path.write_text(f"{mode} {len(temperatures)}\n" +
                    " ".join(f"{temperature:.17g}" for temperature in temperatures) + "\n")


def pilot() -> None:
    levels, radiation, collisions = source_atom()
    stout_differences = stout_level_energy_difference(levels)
    temperatures = np.asarray([1000.0, FLOOR, 10000.0, 100000.0])
    source = full_excitation_balance(levels, radiation, collisions,
                                     temperatures, DENSITY)
    minimum = []
    for density in DENSITY:
        cooling, populations = source[density]
        if (not np.all(np.isfinite(cooling)) or np.any(cooling < 0)
            or not np.all(np.isfinite(populations))
            or np.min(populations) < -1e-13
            or np.max(abs(populations.sum(axis=1) - 1)) > 1e-12):
            raise ValueError("Ne I pilot equilibrium is invalid")
        minimum.append(float(np.min(cooling)))
    save_json(HERE / "pilot_step.json", {
        "passed": True, "temperatures_K": temperatures.tolist(),
        "densities_cm3": DENSITY, "levels": 31,
        "radiative_rows": len(radiation), "collision_rows": len(collisions),
        "maximum_stout_chianti_energy_difference_cm1": max(stout_differences.values()),
        "stout_chianti_rounded_energy_differences_cm1": {
            str(index): value for index, value in stout_differences.items() if value > 0},
        "minimum_cooling_erg_cm3_s": min(minimum),
    })


def native() -> None:
    started = time.monotonic()
    candidate = json.loads((HERE / "candidate_result.json").read_text())
    if not candidate["passed"] or digest(ATOM) != candidate["candidate_sha256"]:
        raise ValueError("Prepared Ne I atom changed")
    levels, _, collisions = source_atom()
    probe = np.asarray([1000.0, 1000.35, FLOOR, 1000.37, 100000.0])
    strength_grid = HERE / "native_strength_grid.txt"
    profile_grid = HERE / "native_profile_grid.txt"
    strengths = HERE / "native_strengths.txt"
    profile = HERE / "native_profiles.txt"
    if strengths.exists() or profile.exists():
        raise FileExistsError("Preserve existing native profiles")
    grid_file(strength_grid, 1, probe)
    grid_file(profile_grid, 2, TEMPERATURE)
    subprocess.run([str(DRIVER), str(ATOM), str(strength_grid), str(strengths)],
                   check=True, timeout=300)
    measured = np.loadtxt(strengths)
    collisions.sort(key=lambda row: tuple(sorted((row["ll"], row["ul"]),
                    key=lambda index: (levels[index][1], index))))
    expected = np.column_stack([C.upsilon(row, probe) for row in collisions])
    if measured.shape != (5, 466) or not np.allclose(measured[:, 0], probe,
                                                       rtol=1e-13, atol=0):
        raise ValueError("Incomplete native source-strength profile")
    positive = expected > 0
    if np.any(measured[:, 1:][~positive] != 0):
        raise ValueError("Native zero effective strengths changed")
    strength_error = float(np.max(abs(measured[:, 1:][positive] /
                                       expected[positive] - 1)))
    if strength_error > 1e-10:
        raise ValueError("Native collision strengths differ from raw source")
    subprocess.run([str(DRIVER), str(ATOM), str(profile_grid), str(profile)],
                   check=True, timeout=900)
    rows = np.loadtxt(profile)
    if rows.shape != (1604, 34):
        raise ValueError("Incomplete Ne I native population profile")
    if (not np.allclose(rows[:, 0], np.repeat(TEMPERATURE, 4), rtol=1e-13, atol=0)
        or not np.array_equal(rows[:, 1], np.tile(DENSITY, 401))
        or not np.all(np.isfinite(rows)) or np.any(rows[:, 2] < 0)
        or np.min(rows[:, 3:]) < -1e-13
        or np.max(abs(rows[:, 3:].sum(axis=1) - 1)) > 1e-12):
        raise ValueError("Invalid native Ne I cooling or populations")
    save_json(HERE / "native_step.json", {
        "passed": True, "temperature_points": 401, "densities_cm3": DENSITY,
        "nlev": 31, "nrad": 144, "nups": 465,
        "source_strength_samples": int(expected.size),
        "maximum_relative_strength_error": strength_error,
        "exact_zero_strength_samples": int(np.sum(~positive)),
        "native_seconds": time.monotonic() - started,
        "atom_sha256": digest(ATOM), "profile_sha256": digest(profile),
        "driver_sha256": digest(DRIVER),
        "native_object_sha256": digest(DRIVER.parent / "nlevel_mod.o"),
    })


def direct() -> None:
    started = time.monotonic()
    levels, radiation, collisions = source_atom()
    source = full_excitation_balance(levels, radiation, collisions,
                                     TEMPERATURE, DENSITY)
    payload = {"temperature_K": TEMPERATURE, "density_cm3": np.asarray(DENSITY)}
    minimum = []
    for index, density in enumerate(DENSITY):
        cooling, populations = source[density]
        if (cooling.shape != (401,) or populations.shape != (401, 31)
            or not np.all(np.isfinite(cooling)) or np.any(cooling < 0)
            or not np.all(np.isfinite(populations))
            or np.min(populations) < -1e-13
            or np.max(abs(populations.sum(axis=1) - 1)) > 1e-12):
            raise ValueError("Invalid independent Ne I profile")
        payload[f"cooling_{index}"] = cooling
        payload[f"populations_{index}"] = populations
        minimum.append(float(np.min(cooling)))
    output = HERE / "direct_profiles.npz"
    if output.exists():
        raise FileExistsError("Preserve existing independent profile")
    np.savez_compressed(output, **payload)
    save_json(HERE / "direct_step.json", {
        "passed": True, "temperature_points": 401, "densities_cm3": DENSITY,
        "nlev": 31, "nrad": len(radiation), "nups": len(collisions),
        "minimum_cooling_erg_cm3_s": min(minimum),
        "direct_seconds": time.monotonic() - started,
        "profile_sha256": digest(output),
        "independent_radiation": "Raw CHIANTI wgfa plus exact two Stout tp rates",
    })


def compare() -> None:
    with np.load(HERE / "direct_profiles.npz") as source:
        if (not np.array_equal(source["temperature_K"], TEMPERATURE)
            or not np.array_equal(source["density_cm3"], DENSITY)):
            raise ValueError("Independent profile grid changed")
        native_rows = np.loadtxt(HERE / "native_profiles.txt")
        if native_rows.shape != (1604, 34):
            raise ValueError("Native profile count changed")
        cooling_errors = []
        population_errors = []
        fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True,
                                 constrained_layout=True)
        difference, difference_axes = plt.subplots(2, 2, figsize=(9, 7),
                                                    sharex=True, constrained_layout=True)
        for index, density in enumerate(DENSITY):
            cooling = source[f"cooling_{index}"]
            populations = source[f"populations_{index}"]
            measured = native_rows[index::4, 2]
            measured_pop = native_rows[index::4, 3:]
            positive = cooling > 64 * np.finfo(float).tiny
            if np.any(measured[cooling == 0] != 0):
                raise ValueError("Native zero cooling differs")
            cooling_error = float(np.max(abs(measured[positive] /
                                             cooling[positive] - 1)))
            population_error = float(np.max(abs(measured_pop - populations)))
            cooling_errors.append(cooling_error)
            population_errors.append(population_error)
            panel = axes.flat[index]
            panel.loglog(TEMPERATURE, cooling, label="Raw CHIANTI plus Stout A")
            panel.loglog(TEMPERATURE, measured, "--", label="Compiled source atom")
            panel.set_title(r"$n_e=" + f"{density:g}" + r"\,\mathrm{cm}^{-3}$")
            panel.set_xlabel(r"$T$ [K]")
            panel.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
            panel.legend(fontsize=8)
            panel = difference_axes.flat[index]
            panel.semilogx(TEMPERATURE[positive],
                           100 * (measured[positive] / cooling[positive] - 1))
            panel.set_title(r"$n_e=" + f"{density:g}" + r"\,\mathrm{cm}^{-3}$")
            panel.set_xlabel(r"$T$ [K]")
            panel.set_ylabel("Native / direct - 1 [percent]")
        if max(cooling_errors) > 1e-10 or max(population_errors) > 1e-10:
            raise ValueError("Native Ne I source calculation differs")
        fig.savefig(HERE / "cooling_profiles.pdf")
        difference.savefig(HERE / "native_differences.pdf")
        plt.close(fig)
        plt.close(difference)
    save_json(HERE / "comparison_step.json", {
        "passed": True, "temperature_points": 401, "densities_cm3": DENSITY,
        "maximum_relative_cooling_errors": cooling_errors,
        "maximum_absolute_population_errors": population_errors,
        "source_atom_sha256": digest(ATOM),
        "native_profile_sha256": digest(HERE / "native_profiles.txt"),
        "direct_profile_sha256": digest(HERE / "direct_profiles.npz"),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", required=True,
                        choices=("pilot", "native", "direct", "compare"))
    args = parser.parse_args()
    globals()[args.step]()


if __name__ == "__main__":
    main()
