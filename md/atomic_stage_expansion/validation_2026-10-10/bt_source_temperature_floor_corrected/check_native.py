#!/usr/bin/env python3
"""Check cold BT strengths, Maxwellian rates, and optional native metadata."""

import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys

os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C
from validate_source_excitation_atom import grid_file

DRIVER = BASE / "build/check_source_excitation_atom.x"
SOURCE_TMIN = 3200.
GAP_CM1 = 1000.
X = np.array([0., .23, .51, .79, 1.])
Y = np.array([.1, .3, .2, .5, .4])
Y6 = np.array([-1., -.8, -1.1, -.5, -.4])
STRENGTH_T = np.array([500., 1000., 1800., 3199.9, 3200., 3200.1,
                       5000., 20000., 100000.])
PROFILE_T = np.array([500., 1000., 1800., 3199.9, 3200., 3200.1,
                      5000., 20000.])
DENSITY = [1e-6, 100., 1e6, 1e10]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_record(kind, minimum=True):
    record = dict(ll=1, ul=2, de=.02, cups=2., ttype=kind,
                  xs=X.copy(), ys=Y6.copy() if kind == 6 else Y.copy())
    if minimum:
        record["bt_tmin_K"] = SOURCE_TMIN
    return record


def atom_text(kind, minimum=True):
    y = Y6 if kind == 6 else Y
    lines = ["NLEVEL_VERSION 2", "NLEV 2", "1 1 0", f"2 3 {GAP_CM1:.17g}",
             "NRAD 1", "1 2 1", f"NODES {len(X)}", "NUPS 1",
             f"1 2 BT_SOURCE {kind} 2 0.02 CUBIC {len(X)}"]
    lines += [f"{x:.17g} {ordinate:.17g}" for x, ordinate in zip(X, y)]
    if minimum:
        lines.append(f"BT_TMIN_K {SOURCE_TMIN:.17g}")
    return "\n".join(lines) + "\n"


def run_atom(name, contents, grid, expect_success=True):
    atom = BASE / f"{name}.txt"
    output = BASE / f"{name}_output.txt"
    atom.write_text(contents)
    command = [str(DRIVER), str(atom), str(grid), str(output)]
    run = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, timeout=30, check=False)
    (BASE / f"{name}.log").write_text(run.stdout)
    if expect_success and run.returncode:
        raise AssertionError(f"{name}: native exit {run.returncode}: {run.stdout[-500:]}")
    if not expect_success and run.returncode == 0:
        raise AssertionError(f"{name}: malformed input accepted")
    if not expect_success and ("ERROR:" not in run.stdout or
                               "BT_TMIN_K" not in run.stdout):
        raise AssertionError(f"{name}: intended reader rejection absent")
    return np.loadtxt(output) if expect_success else run.returncode


def check_rates(rows, record, kind):
    gap = GAP_CM1 * C.HC_ERG_CM
    maximum_cooling = 0.
    maximum_population = 0.
    maximum_balance = 0.
    for index, temperature in enumerate(PROFILE_T):
        strength = C.upsilon(record, [temperature])[0]
        downward = C.COLL_PREF * strength / (3 * np.sqrt(temperature))
        upward = C.COLL_PREF * strength / np.sqrt(temperature) * \
                 np.exp(-gap / (C.K_B_ERG * temperature))
        balance = upward / downward / (3 * np.exp(-gap / (C.K_B_ERG * temperature)))
        maximum_balance = max(maximum_balance, abs(balance - 1))
        for j, density in enumerate(DENSITY):
            row = rows[4 * index + j]
            expected_upper = density * upward / (1 + density * (upward + downward))
            expected_cooling = expected_upper * gap / density
            measured_cooling = row[2]
            measured_upper = row[4]
            maximum_cooling = max(maximum_cooling,
                                  abs(measured_cooling / expected_cooling - 1))
            maximum_population = max(maximum_population,
                                     abs(measured_upper - expected_upper))
            if not np.isclose(row[0], temperature, rtol=1e-13, atol=0):
                raise AssertionError(f"Type {kind}: wrong temperature")
            if row[1] != density:
                raise AssertionError(f"Type {kind}: wrong density")
    if maximum_cooling > 1e-10 or maximum_population > 1e-10 or maximum_balance > 1e-14:
        raise AssertionError(f"Type {kind}: rate balance failed")
    return dict(maximum_relative_cooling_error=maximum_cooling,
                maximum_absolute_population_error=maximum_population,
                maximum_relative_detailed_balance_error=maximum_balance)


def main():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    sources = [Path(__file__), ROOT / "src/define.f90", ROOT / "src/nlevel_mod.f90",
               ROOT / "tools/fitting/chianti_cooling.py",
               ROOT / "tests/ion_stages/check_source_excitation_atom.f90", DRIVER]
    before = {str(path): digest(path) for path in sources}
    strength_grid = BASE / "strength_grid.txt"
    profile_grid = BASE / "profile_grid.txt"
    grid_file(strength_grid, 1, STRENGTH_T)
    grid_file(profile_grid, 2, PROFILE_T)
    results = {}
    type5_cold_strength = None
    for kind in range(1, 7):
        name = f"type{kind}"
        record = source_record(kind)
        measured = run_atom(name, atom_text(kind), strength_grid)
        expected = C.upsilon(record, STRENGTH_T)
        if measured.shape != (len(STRENGTH_T), 2):
            raise AssertionError(f"{name}: incomplete strengths")
        difference = np.abs(measured[:, 1] / expected - 1)
        if np.max(difference) > 1e-11:
            raise AssertionError(f"{name}: Python/native strengths differ")
        if np.any(measured[:5, 1] != measured[4, 1]):
            raise AssertionError(f"{name}: Upsilon is not exactly constant below Tmin")
        if kind == 5:
            type5_cold_strength = measured[0, 1]
        profile = run_atom(name + "_profile", atom_text(kind), profile_grid)
        if profile.shape != (4 * len(PROFILE_T), 5):
            raise AssertionError(f"{name}: incomplete profile")
        rates = check_rates(profile, record, kind)
        results[name] = {"maximum_relative_strength_error": float(np.max(difference)),
                         **rates}

    without_floor = run_atom("no_floor", atom_text(5, False), strength_grid)
    no_floor_expected = C.upsilon(source_record(5, False), STRENGTH_T)
    no_floor_error = float(np.max(np.abs(without_floor[:, 1] / no_floor_expected - 1)))
    if no_floor_error > 1e-11 or without_floor[0, 1] == type5_cold_strength:
        raise AssertionError("Version-2 source without floor changed")

    legacy = ("NLEV 2\n1 1 0\n2 3 1000\nNRAD 1\n1 2 1\nNUPS 1\n"
              "1 2 2 1 0.02 0 1 0 1 0.3\n")
    legacy_values = run_atom("legacy", legacy, strength_grid)
    if not np.allclose(legacy_values[:, 1], .3, rtol=1e-14, atol=0):
        raise AssertionError("Legacy excitation changed")

    base = atom_text(2)
    table = ("NLEVEL_VERSION 2\nNLEV 2\n1 1 0\n2 3 1000\nNRAD 1\n1 2 1\n"
             "NODES 2\nNUPS 1\n1 2 TEMPERATURE_TABLE 0 0 0 LOG_PCHIP_T 2\n"
             "1000 0.1\n10000 0.2\nBT_TMIN_K 3200\n")
    malformed = {
        "duplicate": base + "BT_TMIN_K 3200\n",
        "zero": base.replace("BT_TMIN_K 3200", "BT_TMIN_K 0"),
        "negative": base.replace("BT_TMIN_K 3200", "BT_TMIN_K -1"),
        "nan": base.replace("BT_TMIN_K 3200", "BT_TMIN_K NaN"),
        "inf": base.replace("BT_TMIN_K 3200", "BT_TMIN_K Infinity"),
        "nonnumeric": base.replace("BT_TMIN_K 3200", "BT_TMIN_K abc"),
        "extra": base.replace("BT_TMIN_K 3200", "BT_TMIN_K 3200 extra"),
        "before_nups": atom_text(2, False).replace(
            "NUPS 1", "BT_TMIN_K 3200\nNUPS 1"),
        "after_density": base.replace("BT_TMIN_K 3200", "DENSITY_RATIO 1\nBT_TMIN_K 3200"),
        "legacy": legacy + "BT_TMIN_K 3200\n",
        "without_bt": table,
    }
    rejected = {name: run_atom("reject_" + name, contents, strength_grid, False)
                for name, contents in malformed.items()}
    after = {str(path): digest(path) for path in sources}
    result = {"passed": before == after, "sources_unchanged": before == after,
              "source_sha256": before, "types": results,
              "old_v2_maximum_relative_strength_error": no_floor_error,
              "legacy_strength": float(legacy_values[0, 1]),
              "malformed_rejected": rejected,
              "temperature_K": STRENGTH_T.tolist(),
              "profile_temperature_K": PROFILE_T.tolist(),
              "densities_cm3": DENSITY}
    (BASE / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    if before != after:
        raise RuntimeError("Source changed during native test")
    print(json.dumps({"passed": result["passed"], "types": len(results),
                      "malformed_rejected": len(rejected)}))


if __name__ == "__main__":
    main()
