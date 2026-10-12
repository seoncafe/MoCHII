"""Check stationary populations and electron-to-photon power for full O IV."""

import hashlib
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = HERE / "candidate"
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C  # noqa: E402
from validate_source_excitation_atom import source_subset  # noqa: E402


TEMPERATURE = np.geomspace(1000.0, 100000.0, 401)
DENSITY = np.asarray([1e-6, 100.0, 1e6, 1e10])
LEVEL_COUNT = 147


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def maximum_record(values, level_axis=False):
    index = np.unravel_index(int(np.argmax(values)), values.shape)
    result = {"value": float(values[index]),
              "temperature_K": float(TEMPERATURE[index[0]]),
              "density_cm3": float(DENSITY[index[1]])}
    if level_axis:
        result["level"] = int(index[2] + 1)
    return result


def main():
    output = BASE / "ledger_result.json"
    if output.exists():
        raise FileExistsError(output)
    native_record = json.loads((BASE / "native_result.json").read_text())
    validation = native_record["source_result"]
    profile = BASE / "native/native_profiles.txt"
    atom = BASE / "native/nlevel_o_4.txt"
    if (native_record.get("passed") is not True
        or native_record["products_sha256"].get(str(profile)) != digest(profile)
        or validation["product_sha256"] != digest(atom)):
        raise ValueError("Full O IV saved native profile identity changed")
    levels, radiation, collisions = source_subset("o", 4, validation["selected_ids"])
    if (len(levels), len(radiation), len(collisions)) != (147, 2319, 10731):
        raise ValueError("Full O IV source balance rows changed")
    native = np.loadtxt(profile)
    if (native.shape != (1604, LEVEL_COUNT + 3)
        or not np.allclose(native[:, 0], np.repeat(TEMPERATURE, 4), rtol=1e-13, atol=0)
        or not np.array_equal(native[:, 1], np.tile(DENSITY, 401))
        or not np.all(np.isfinite(native))):
        raise ValueError("Incomplete O IV native populations")
    populations = native[:, 3:].reshape(401, 4, LEVEL_COUNT)
    if (np.min(populations) < -1e-13
        or np.max(abs(populations.sum(axis=2) - 1)) > 1e-12):
        raise ValueError("O IV populations are not normalized")
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
        lower_weight, lower_energy = levels[lower]
        upper_weight, upper_energy = levels[upper]
        gap = (upper_energy - lower_energy) * C.HC_ERG_CM
        if gap < 0:
            raise ValueError("Negative O IV collision energy")
        strength = C.upsilon(row, TEMPERATURE)[:, None]
        if not np.all(np.isfinite(strength)) or np.any(strength < 0):
            raise ValueError("Invalid O IV source strength")
        common = C.COLL_PREF * strength / root_temperature
        excitation = common / lower_weight * np.exp(
            -gap / (C.K_B_ERG * TEMPERATURE[:, None]))
        deexcitation = common / upper_weight
        forward = DENSITY[None, :] * populations[:, :, lower - 1] * excitation
        reverse = DENSITY[None, :] * populations[:, :, upper - 1] * deexcitation
        net, gross = forward - reverse, forward + reverse
        residual[:, :, lower - 1] -= net
        residual[:, :, upper - 1] += net
        gross_level_flow[:, :, lower - 1] += gross
        gross_level_flow[:, :, upper - 1] += gross
        gross_collision_flow += gross
        collision_power += net * gap
        gross_collision_power += gross * gap
    radiative_power = np.zeros((401, 4))
    for lower, upper, _, _, rate in radiation:
        if rate <= 0 or levels[upper][1] <= levels[lower][1]:
            raise ValueError("Invalid O IV radiative direction")
        flow = populations[:, :, upper - 1] * rate
        residual[:, :, lower - 1] += flow
        residual[:, :, upper - 1] -= flow
        gross_level_flow[:, :, lower - 1] += flow
        gross_level_flow[:, :, upper - 1] += flow
        radiative_power += flow * (levels[upper][1] - levels[lower][1]) * C.HC_ERG_CM
    if np.any((gross_level_flow == 0) & (residual != 0)):
        raise ValueError("O IV nonzero level residual without physical flow")
    local_stationarity = np.zeros_like(residual)
    np.divide(abs(residual), gross_level_flow, out=local_stationarity,
              where=gross_level_flow > 0)
    if np.any((gross_collision_flow == 0) &
              (np.max(abs(residual), axis=2) != 0)):
        raise ValueError("O IV residual without collisional flow")
    global_stationarity = np.zeros_like(residual)
    np.divide(abs(residual), gross_collision_flow[:, :, None],
              out=global_stationarity,
              where=gross_collision_flow[:, :, None] > 0)
    significant = gross_level_flow > 1e-12 * gross_collision_flow[:, :, None]
    significant_local = np.where(significant, local_stationarity, 0)
    if np.any((gross_collision_power == 0) & (collision_power != radiative_power)):
        raise ValueError("O IV energy residual without collision energy")
    energy_error = np.zeros_like(collision_power)
    np.divide(abs(collision_power - radiative_power), gross_collision_power,
              out=energy_error, where=gross_collision_power > 0)
    native_power = DENSITY[None, :] * native[:, 2].reshape(401, 4)
    if np.any((radiative_power == 0) & (native_power != 0)):
        raise ValueError("O IV native cooling without radiation")
    cooling_error = np.zeros_like(radiative_power)
    np.divide(abs(radiative_power - native_power), radiative_power,
              out=cooling_error, where=radiative_power > 0)
    if not all(np.all(np.isfinite(array)) for array in (
        residual, gross_level_flow, collision_power, radiative_power,
        global_stationarity, significant_local, energy_error, cooling_error)):
        raise ValueError("Nonfinite O IV energy ledger")
    passed = (np.max(global_stationarity) < 1e-8
              and np.max(significant_local) < 1e-8
              and np.max(energy_error) < 1e-8
              and np.max(cooling_error) < 1e-8)
    report = {"passed": bool(passed), "ion": "OIV", "atom_sha256": digest(atom),
              "native_profile_sha256": digest(profile), "temperature_points": 401,
              "densities_cm3": DENSITY.tolist(),
              "maximum_global_stationarity": maximum_record(global_stationarity, True),
              "maximum_significant_local_stationarity": maximum_record(
                  significant_local, True),
              "maximum_energy_error": maximum_record(energy_error),
              "maximum_cooling_power_error": maximum_record(cooling_error),
              "source_accuracy_verified": False,
              "physical_model_accepted": False}
    output.write_text(json.dumps(report, indent=2) + "\n")
    if not passed:
        raise ValueError("Full O IV stationary or energy balance failed")


if __name__ == "__main__":
    main()
