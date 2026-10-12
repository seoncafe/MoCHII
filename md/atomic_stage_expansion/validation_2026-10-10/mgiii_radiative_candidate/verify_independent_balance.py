"""Audit saved Mg III native populations against source collisional and radiative flows."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C  # noqa: E402


ATOM = HERE / "nlevel_mg_3_line_strength_candidate.txt"
NATIVE = HERE / "profile_native_full.txt"
SOURCE = Path.home() / "RT_Codes/CHIANTI/dbase/mg/mg_3/mg_3.scups"
RESULT = HERE / "independent_balance_result.json"
DENSITY = np.asarray([1e-6, 1e2, 1e6, 1e10])
TEMPERATURE = np.geomspace(1000.0, 100000.0, 801)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atom_records() -> tuple[dict[int, tuple[float, float]], list[tuple[int, int, float]], int]:
    lines = ATOM.read_text().splitlines()
    begin = lines.index("NLEVEL_VERSION 2") + 1
    fields = lines[begin].split()
    assert fields == ["NLEV", "167"]
    levels = {}
    for line in lines[begin + 1:begin + 168]:
        index, weight, energy = line.split()
        levels[int(index)] = (float(weight), float(energy))
    assert sorted(levels) == list(range(1, 168))
    begin += 168
    assert lines[begin] == "NRAD 3364"
    radiation = []
    for line in lines[begin + 1:begin + 3365]:
        lower, upper, rate = line.split()
        lower, upper, rate = int(lower), int(upper), float(rate)
        assert 0 < lower <= 167 and 0 < upper <= 167
        assert levels[upper][1] > levels[lower][1] and rate > 0
        radiation.append((lower, upper, rate))
    assert lines.count("BT_TMIN_K 1800") == 1
    assert "NODES 207915" in lines and "NUPS 13861" in lines
    return levels, radiation, 1800


def maximum_record(values: np.ndarray, *, level_axis: bool = False) -> dict:
    index = np.unravel_index(int(np.argmax(values)), values.shape)
    record = {
        "value": float(values[index]),
        "temperature_K": float(TEMPERATURE[index[0]]),
        "density_cm3": float(DENSITY[index[1]]),
    }
    if level_axis:
        record["level"] = int(index[2] + 1)
    return record


def main() -> None:
    if RESULT.exists():
        raise FileExistsError(f"Preserve existing ledger result: {RESULT}")
    previous = json.loads((HERE / "recovery_path_fix/audit_completed_results.json").read_text())
    assert previous["passed"] and previous["temperature_points"] == 801
    assert digest(ATOM) == previous["verified_product_hashes"]["full_atom"]
    assert digest(NATIVE) == previous["verified_product_hashes"]["native_full"]

    levels, radiation, floor = atom_records()
    raw_transitions = C.read_scups(SOURCE)
    if len(raw_transitions) != 21736:
        raise ValueError("The complete raw Mg III collision-row count changed")
    transitions = [row for row in raw_transitions
                   if row["ll"] in levels and row["ul"] in levels]
    if len(transitions) != 13861:
        raise ValueError("The raw Mg III collision-row count changed")
    for transition in transitions:
        transition["bt_tmin_K"] = floor
        if transition["ll"] not in levels or transition["ul"] not in levels:
            raise ValueError("Raw collision row is outside the saved atom")
    native = np.loadtxt(NATIVE)
    if native.shape != (3204, 170):
        raise ValueError("Incomplete saved native population profile")
    if not np.allclose(native[:, 0], np.repeat(TEMPERATURE, 4), rtol=1e-13, atol=0):
        raise ValueError("Saved temperature grid differs")
    if not np.array_equal(native[:, 1], np.tile(DENSITY, 801)):
        raise ValueError("Saved density grid differs")
    if not np.all(np.isfinite(native)):
        raise ValueError("Saved native profile is nonfinite")
    populations = native[:, 3:].reshape(801, 4, 167)
    if (np.min(populations) < -1e-13
            or np.max(abs(populations.sum(axis=2) - 1)) > 1e-12):
        raise ValueError("Saved populations are invalid")

    residual = np.zeros_like(populations)
    gross_level_flow = np.zeros_like(populations)
    gross_collision_flow = np.zeros((801, 4))
    collision_power = np.zeros((801, 4))
    gross_collision_power = np.zeros((801, 4))
    root_T = np.sqrt(TEMPERATURE)[:, None]
    for transition in transitions:
        first, second = transition["ll"], transition["ul"]
        lower, upper = sorted((first, second),
                              key=lambda index: (levels[index][1], index))
        g_lower, e_lower = levels[lower]
        g_upper, e_upper = levels[upper]
        gap = (e_upper - e_lower) * C.HC_ERG_CM
        if gap < 0:
            raise ValueError("Collision has negative physical energy gap")
        strength = C.upsilon(transition, TEMPERATURE)[:, None]
        if not np.all(np.isfinite(strength)) or np.any(strength < 0):
            raise ValueError("Raw effective collision strength is invalid")
        common = C.COLL_PREF * strength / root_T
        excitation = common / g_lower * np.exp(-gap / (C.K_B_ERG * TEMPERATURE[:, None]))
        deexcitation = common / g_upper
        forward = DENSITY[None, :] * populations[:, :, lower - 1] * excitation
        reverse = DENSITY[None, :] * populations[:, :, upper - 1] * deexcitation
        net = forward - reverse
        gross = forward + reverse
        residual[:, :, lower - 1] -= net
        residual[:, :, upper - 1] += net
        gross_level_flow[:, :, lower - 1] += gross
        gross_level_flow[:, :, upper - 1] += gross
        gross_collision_flow += gross
        collision_power += net * gap
        gross_collision_power += gross * gap

    radiative_power = np.zeros((801, 4))
    for lower, upper, rate in radiation:
        flow = populations[:, :, upper - 1] * rate
        residual[:, :, lower - 1] += flow
        residual[:, :, upper - 1] -= flow
        gross_level_flow[:, :, lower - 1] += flow
        gross_level_flow[:, :, upper - 1] += flow
        gap = (levels[upper][1] - levels[lower][1]) * C.HC_ERG_CM
        radiative_power += flow * gap

    # Divide only where the corresponding gross physical flow is positive.
    # A zero gross flow must have an exactly zero residual.
    if np.any((gross_level_flow == 0) & (residual != 0)):
        raise ValueError("Stationarity residual with zero gross level flow")
    relative_stationarity = np.zeros_like(residual)
    np.divide(abs(residual), gross_level_flow, out=relative_stationarity,
              where=gross_level_flow > 0)
    if np.any((gross_collision_flow == 0) & (np.max(abs(residual), axis=2) != 0)):
        raise ValueError("Stationarity residual with zero total gross collision flow")
    global_relative_stationarity = np.zeros_like(residual)
    np.divide(abs(residual), gross_collision_flow[:, :, None],
              out=global_relative_stationarity,
              where=gross_collision_flow[:, :, None] > 0)
    significant_level = gross_level_flow > 1e-12 * gross_collision_flow[:, :, None]
    significant_local_stationarity = np.where(significant_level,
                                              relative_stationarity, 0)
    if np.any((gross_collision_power == 0) & (collision_power != radiative_power)):
        raise ValueError("Energy residual with zero gross collision energy flow")
    relative_energy = np.zeros_like(collision_power)
    np.divide(abs(collision_power - radiative_power), gross_collision_power,
              out=relative_energy, where=gross_collision_power > 0)
    cooling = native[:, 2].reshape(801, 4)
    cooling_power = DENSITY[None, :] * cooling
    if np.any((radiative_power == 0) & (cooling_power != 0)):
        raise ValueError("Native cooling is positive with zero radiative power")
    relative_radiation_to_native = np.zeros_like(cooling_power)
    np.divide(abs(radiative_power - cooling_power), radiative_power,
              out=relative_radiation_to_native, where=radiative_power > 0)
    if not all(np.all(np.isfinite(value)) for value in
               (residual, gross_level_flow, collision_power, radiative_power,
                relative_stationarity, global_relative_stationarity,
                significant_local_stationarity, relative_energy,
                relative_radiation_to_native)):
        raise ValueError("Nonfinite ledger value")

    result = {
        "passed": bool(np.max(global_relative_stationarity) < 1e-8
                       and np.max(significant_local_stationarity) < 1e-8
                       and np.max(relative_energy) < 1e-8
                       and np.max(relative_radiation_to_native) < 1e-8),
        "scope": "Saved native full-atom populations with raw source collisions and corrected atom A",
        "temperature_count": 801,
        "temperature_K": [1000.0, 100000.0],
        "density_cm3": DENSITY.tolist(),
        "levels": 167,
        "collision_rows": len(transitions),
        "radiative_rows": len(radiation),
        "adopted_constant_upsilon_below_K": floor,
        "maximum_stationarity_residual_relative_to_gross_level_flow":
            maximum_record(relative_stationarity, level_axis=True),
        "maximum_stationarity_residual_relative_to_total_gross_collision_flow":
            maximum_record(global_relative_stationarity, level_axis=True),
        "maximum_local_stationarity_where_level_flow_exceeds_1e-12_total":
            maximum_record(significant_local_stationarity, level_axis=True),
        "local_stationarity_caveat": (
            "Native zero populations at very weak high-excitation levels can "
            "leave local residual/gross flow near one; global gross-flow "
            "normalization and significant-level local normalization assess "
            "the resolved statistical equilibrium"),
        "maximum_collision_minus_radiation_power_relative_to_gross_collision_power":
            maximum_record(relative_energy),
        "maximum_radiation_vs_native_cooling_relative_error":
            maximum_record(relative_radiation_to_native),
        "minimum_positive_radiative_power_erg_s": float(np.min(
            radiative_power[radiative_power > 0])),
        "input_sha256": {str(path): digest(path) for path in
                         (ATOM, NATIVE, SOURCE, ROOT / "tools/fitting/chianti_cooling.py")},
        "no_solver_or_profile_regeneration": True,
    }
    RESULT.write_text(json.dumps(result, indent=2) + "\n")
    if not result["passed"]:
        raise ValueError("Independent saved-population balance exceeded tolerance")


if __name__ == "__main__":
    main()
