#!/usr/bin/env python3
"""Compare installed five-state atoms with complete bound excitation sources."""

import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C
from validate_excitation_density_profiles import full_excitation_balance
from validate_source_excitation_atom import grid_file

BASE = Path(__file__).resolve().parent
DRIVER = (ROOT / "md/atomic_stage_expansion/validation_2026-10-10/"
          "full_bound_native/optimized_build/check_source_excitation_atom.x")
DENSITY = [1e-6, 100., 1e6, 1e10]
TEMPERATURE = np.geomspace(1e3, 1e5, 801)
IONS = (("CI", "c", 1), ("CII", "c", 2), ("NI", "n", 1))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def maximum_relative(actual, reference):
    normal = reference > 64 * np.finfo(float).tiny
    if not np.any(normal):
        raise ValueError("No normal positive source cooling")
    if np.any((reference > 0) & (actual == 0)):
        raise ValueError("Positive source cooling became zero")
    difference = np.abs(actual[normal] / reference[normal] - 1)
    index = np.argmax(difference)
    return dict(maximum=float(difference[index]),
                temperature_K=float(TEMPERATURE[np.flatnonzero(normal)[index]]),
                normal_points=int(normal.sum()),
                zero_points=int(np.sum(reference == 0)))


def inspect_atom(label, element, stage):
    work = BASE / label
    work.mkdir(parents=True, exist_ok=False)
    levels, radiation, collisions = C.bound_excitation_atom(element, stage)
    selected_ids = sorted(levels)[:5]
    chosen = {i: levels[i] for i in selected_ids}
    chosen_radiation = [r for r in radiation if r[0] in chosen and r[1] in chosen]
    chosen_collisions = [r for r in collisions if r["ll"] in chosen and r["ul"] in chosen]
    nonpositive_rad = [(l, u) for l, u, _, _, a in radiation
                       if a > 0 and levels[u][1] <= levels[l][1]]
    zero_collision = [(tr["ll"], tr["ul"]) for tr in collisions
                      if levels[tr["ll"]][1] == levels[tr["ul"]][1]]
    if nonpositive_rad:
        raise ValueError(f"{label}: nonpositive radiative gaps {nonpositive_rad}")

    full = full_excitation_balance(levels, radiation, collisions,
                                   TEMPERATURE, DENSITY)
    chosen_result = full_excitation_balance(chosen, chosen_radiation,
                                            chosen_collisions, TEMPERATURE, DENSITY)

    atom_file = ROOT / f"data/atomic/nlevel_{element}_{stage}.txt"
    grid = work / "grid.txt"
    native_file = work / "installed_native.txt"
    grid_file(grid, 2, TEMPERATURE)
    subprocess.run([str(DRIVER), str(atom_file), str(grid), str(native_file)],
                   check=True, timeout=300)
    rows = np.loadtxt(native_file)
    if rows.shape != (4 * TEMPERATURE.size, 8):
        raise ValueError(f"{label}: incomplete installed native profile")
    if not np.allclose(rows[:, 0], np.repeat(TEMPERATURE, 4), rtol=1e-13, atol=0):
        raise ValueError(f"{label}: wrong native temperature grid")
    if not np.array_equal(rows[:, 1], np.tile(DENSITY, TEMPERATURE.size)):
        raise ValueError(f"{label}: wrong native density grid")

    saved = {"temperature_K": TEMPERATURE}
    comparison = []
    for j, density in enumerate(DENSITY):
        reference, _ = full[density]
        first_five, _ = chosen_result[density]
        native = rows[j::4, 2]
        if (not np.all(np.isfinite(reference)) or not np.all(np.isfinite(first_five))
                or not np.all(np.isfinite(native)) or min(np.min(reference),
                np.min(first_five), np.min(native)) < 0):
            raise ValueError(f"{label}: invalid cooling")
        comparison.append({"density_cm3": density,
                           "five_source_vs_full": maximum_relative(first_five, reference),
                           "installed_native_vs_full": maximum_relative(native, reference),
                           "installed_native_vs_five_source": maximum_relative(native, first_five),
                           "minimum_positive_full": float(np.min(reference[reference > 0]))})
        saved[f"full_cooling_density_{j}"] = reference
        saved[f"five_cooling_density_{j}"] = first_five
        saved[f"installed_cooling_density_{j}"] = native
    np.savez_compressed(work / "profiles.npz", **saved)
    return {"full_levels": len(levels), "full_radiative_rows": len(radiation),
            "full_collision_rows": len(collisions), "installed_levels": 5,
            "selected_source_ids": selected_ids,
            "five_radiative_rows": len(chosen_radiation),
            "five_collision_rows": len(chosen_collisions),
            "zero_energy_collision_pairs": zero_collision,
            "nonpositive_radiative_gap_pairs": nonpositive_rad,
            "temperature_points": len(TEMPERATURE), "densities_cm3": DENSITY,
            "comparison": comparison, "profiles": str(work / "profiles.npz"),
            "profiles_sha256": digest(work / "profiles.npz")}


def main():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    if not DRIVER.is_file():
        raise FileNotFoundError(DRIVER)
    source_files = [Path(__file__), DRIVER, ROOT / "tools/fitting/chianti_cooling.py",
                    ROOT / "tools/fitting/validate_excitation_density_profiles.py",
                    ROOT / "tools/fitting/atomic_excitation_energies.py",
                    ROOT / "tools/fitting/ionization_potentials.py"]
    from ionization_potentials import SOURCE as ionization_source
    source_files.append(ionization_source)
    for _, element, stage in IONS:
        stem = Path(C.ion_dir(element, stage)) / f"{element}_{stage}"
        source_files.extend(Path(f"{stem}.{extension}")
                            for extension in ("elvlc", "wgfa", "scups"))
        source_files.append(ROOT / f"data/atomic/nlevel_{element}_{stage}.txt")
    before = {str(path): digest(path) for path in source_files}
    ions = {label: inspect_atom(label, element, stage)
            for label, element, stage in IONS}
    after = {str(path): digest(path) for path in source_files}
    result = {"passed": before == after,
              "sources_unchanged": before == after,
              "source_sha256": before, "ions": ions,
              "scope": "Complete bound source versus five-state source and installed native atoms"}
    (BASE / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    if before != after:
        raise RuntimeError("Source changed during the calculation")
    print(json.dumps({"passed": result["passed"], "ions": list(ions)}))


if __name__ == "__main__":
    main()
