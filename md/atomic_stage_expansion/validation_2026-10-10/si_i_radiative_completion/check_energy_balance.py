"""Check Si I native stationary populations and collisional energy balance."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from run_profiles import (C, DENSITY, EXPECTED_LEVELS, HERE, TEMPERATURE,
                          candidate, check_source, digest, radiative_rows, save_json)


RESULT = HERE / "energy_balance_step.json"


def maximum_record(values: np.ndarray, level_axis: bool = False) -> dict:
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
    report, atom, levels, radiation, collisions = candidate()
    check_source(levels, radiation, collisions)
    native_report = json.loads((HERE / "native_step.json").read_text())
    profile_file = HERE / "native_profiles.txt"
    if (native_report.get("passed") is not True or
            native_report["candidate_sha256"] != report["candidate_sha256"] or
            digest(profile_file) != native_report["profile_sha256"]):
        raise ValueError("Si I native populations changed")
    native = np.loadtxt(profile_file)
    if (native.shape != (1604, EXPECTED_LEVELS + 3) or
            not np.allclose(native[:, 0], np.repeat(TEMPERATURE, 4), rtol=1e-13, atol=0) or
            not np.array_equal(native[:, 1], np.tile(DENSITY, 401)) or
            not np.all(np.isfinite(native))):
        raise ValueError("Incomplete Si I native population grid")
    populations = native[:, 3:].reshape(401, 4, EXPECTED_LEVELS)
    if (np.min(populations) < -1e-13 or
            np.max(abs(populations.sum(axis=2) - 1)) > 1e-12):
        raise ValueError("Invalid Si I native populations")
    residual = np.zeros_like(populations)
    gross_level_flow = np.zeros_like(populations)
    gross_collision_flow = np.zeros((401, 4))
    collision_power = np.zeros((401, 4))
    gross_collision_power = np.zeros((401, 4))
    root_temperature = np.sqrt(TEMPERATURE)[:, None]
    for row in collisions:
        first, second = int(row["ll"]), int(row["ul"])
        lower, upper = sorted((first, second),
                              key=lambda index: (levels[index][1], index))
        g_lower, e_lower = levels[lower]
        g_upper, e_upper = levels[upper]
        gap = (e_upper - e_lower) * C.HC_ERG_CM
        if gap < 0:
            raise ValueError("Negative Si I collision photon gap")
        strength = C.upsilon(row, TEMPERATURE)[:, None]
        if not np.all(np.isfinite(strength)) or np.any(strength < 0):
            raise ValueError("Invalid Si I effective collision strength")
        common = C.COLL_PREF * strength / root_temperature
        excitation = common / g_lower * np.exp(-gap / (C.K_B_ERG * TEMPERATURE[:, None]))
        deexcitation = common / g_upper
        forward = np.asarray(DENSITY)[None, :] * populations[:, :, lower - 1] * excitation
        reverse = np.asarray(DENSITY)[None, :] * populations[:, :, upper - 1] * deexcitation
        net = forward - reverse
        gross = forward + reverse
        residual[:, :, lower - 1] -= net
        residual[:, :, upper - 1] += net
        gross_level_flow[:, :, lower - 1] += gross
        gross_level_flow[:, :, upper - 1] += gross
        gross_collision_flow += gross
        collision_power += net * gap
        gross_collision_power += gross * gap
    radiative_power = np.zeros((401, 4))
    for lower, upper, rate in radiative_rows(radiation):
        flow = populations[:, :, upper - 1] * rate
        residual[:, :, lower - 1] += flow
        residual[:, :, upper - 1] -= flow
        gross_level_flow[:, :, lower - 1] += flow
        gross_level_flow[:, :, upper - 1] += flow
        radiative_power += flow * (levels[upper][1] - levels[lower][1]) * C.HC_ERG_CM
    if np.any((gross_level_flow == 0) & (residual != 0)):
        raise ValueError("Nonzero Si I level residual without physical flow")
    local_stationarity = np.zeros_like(residual)
    np.divide(abs(residual), gross_level_flow, out=local_stationarity,
              where=gross_level_flow > 0)
    if np.any((gross_collision_flow == 0) &
              (np.max(abs(residual), axis=2) != 0)):
        raise ValueError("Nonzero Si I stationarity residual without collision flow")
    global_stationarity = np.zeros_like(residual)
    np.divide(abs(residual), gross_collision_flow[:, :, None],
              out=global_stationarity,
              where=gross_collision_flow[:, :, None] > 0)
    significant_level = gross_level_flow > 1e-12 * gross_collision_flow[:, :, None]
    significant_local = np.where(significant_level, local_stationarity, 0)
    if np.any((gross_collision_power == 0) & (collision_power != radiative_power)):
        raise ValueError("Nonzero Si I energy residual without collision energy flow")
    energy_error = np.zeros_like(collision_power)
    np.divide(abs(collision_power - radiative_power), gross_collision_power,
              out=energy_error, where=gross_collision_power > 0)
    native_cooling_power = np.asarray(DENSITY)[None, :] * native[:, 2].reshape(401, 4)
    if np.any((radiative_power == 0) & (native_cooling_power != 0)):
        raise ValueError("Native Si I cooling without radiative power")
    cooling_power_error = np.zeros_like(radiative_power)
    np.divide(abs(radiative_power - native_cooling_power), radiative_power,
              out=cooling_power_error, where=radiative_power > 0)
    quantities = (residual, gross_level_flow, collision_power, radiative_power,
                  local_stationarity, global_stationarity, significant_local,
                  energy_error, cooling_power_error)
    if not all(np.all(np.isfinite(item)) for item in quantities):
        raise ValueError("Nonfinite Si I energy ledger")
    passed = bool(np.max(global_stationarity) < 1e-8 and
                  np.max(significant_local) < 1e-8 and
                  np.max(energy_error) < 1e-8 and
                  np.max(cooling_power_error) < 1e-8)
    save_json(RESULT, {
        "passed": passed,
        "candidate_sha256": report["candidate_sha256"],
        "native_profile_sha256": digest(profile_file),
        "temperature_points": 401, "densities_cm3": DENSITY,
        "maximum_global_stationarity": maximum_record(global_stationarity, True),
        "maximum_significant_local_stationarity": maximum_record(significant_local, True),
        "maximum_energy_error": maximum_record(energy_error),
        "maximum_cooling_power_error": maximum_record(cooling_power_error),
        "source_accuracy_verified": False,
    })
    if not passed:
        raise ValueError("Si I energy balance did not pass")


if __name__ == "__main__":
    main()
