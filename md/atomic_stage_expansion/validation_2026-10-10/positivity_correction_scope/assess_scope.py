"""Read-only reach and source-cooling impact of selective positive interpolation."""
from pathlib import Path
import json
import sys

import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as source
from fit_cooling import IONS as INSTALLED_COOLING_IONS

INVENTORY = ROOT / "md/atomic_stage_expansion/validation_2026-10-10/nev_collision_fit_diagnostic/positivity_inventory.json"
AUDIT = json.loads(INVENTORY.read_text())
PCHIP = {}
NEGATIVE_NODE = {}
for row in AUDIT["affected_transitions"]:
    key = (row["ion"], tuple(sorted(row["pair"])))
    if row["kind"] == "cubic_overshoot":
        PCHIP[key] = row
    else:
        NEGATIVE_NODE[key] = row


def native_pairs(path):
    lines = path.read_text().splitlines()
    position = next(i for i, line in enumerate(lines) if line.startswith("NUPS "))
    count = int(lines[position].split()[1])
    return [tuple(sorted(map(int, line.split()[:2])))
            for line in lines[position+1:position+1+count]]


def legacy_strength(transition, temperature):
    """Reproduce the prior CubicSpline plus zero-clipping implementation."""
    T = np.atleast_1d(np.asarray(temperature, float))
    if "temperature_grid" in transition:
        return CORRECTED_UPSILON(transition, T)
    de, constant, kind = transition["de"], transition["cups"], transition["ttype"]
    et = source.KB_OVER_RY*T/de
    if kind in (1, 4):
        scaled_temperature = 1-np.log(constant)/np.log(et+constant)
    else:
        scaled_temperature = et/(et+constant)
    scaled_temperature = np.clip(scaled_temperature,
                                 transition["xs"][0], transition["xs"][-1])
    scaled_strength = CubicSpline(transition["xs"], transition["ys"])(scaled_temperature)
    if kind == 1:
        strength = scaled_strength*np.log(et+np.e)
    elif kind == 2:
        strength = scaled_strength
    elif kind == 3:
        strength = scaled_strength/(et+1)
    elif kind == 4:
        strength = scaled_strength*np.log(et+constant)
    elif kind == 5:
        strength = scaled_strength/et
    elif kind == 6:
        strength = 10**scaled_strength
    else:
        raise ValueError(f"Unknown transition type {kind}")
    return np.maximum(strength, 0)


CORRECTED_UPSILON = source.upsilon


def cooling_pair(element, stage, temperature):
    source.upsilon = legacy_strength
    try:
        baseline = source.cooling_low_density_limit(element, stage, temperature)
    finally:
        source.upsilon = CORRECTED_UPSILON
    try:
        corrected = source.cooling_low_density_limit(element, stage, temperature)
    finally:
        source.upsilon = CORRECTED_UPSILON
    return baseline, corrected


def main():
    native = []
    for path in sorted((ROOT / "data/atomic").glob("nlevel_*.txt")):
        stem = path.stem.removeprefix("nlevel_")
        fields = stem.split("_")
        if len(fields) < 2:
            continue
        ion = "_".join(fields[:2])
        pairs = native_pairs(path)
        pchip = sorted(set(pair for pair in pairs if (ion, pair) in PCHIP))
        negative = sorted(set(pair for pair in pairs if (ion, pair) in NEGATIVE_NODE))
        native.append({"file": str(path.relative_to(ROOT)), "ion": ion,
                       "native_collision_rows": len(pairs),
                       "pchip_pair_count": len(pchip), "pchip_pairs": pchip,
                       "negative_source_node_pair_count": len(negative),
                       "negative_source_node_pairs": negative})

    temperature = np.logspace(3, 5, 801)
    cooling = []
    for name, (element, stage) in INSTALLED_COOLING_IONS.items():
        ion = f"{element}_{stage}"
        if not any(key[0] == ion for key in PCHIP) and not any(key[0] == ion for key in NEGATIVE_NODE):
            cooling.append({"registered_name": name, "ion": ion,
                            "affected_source_pairs": 0, "status": "unchanged by selection"})
            continue
        try:
            original, updated = cooling_pair(element, stage, temperature)
        except (ValueError, FileNotFoundError) as error:
            cooling.append({"registered_name": name, "ion": ion,
                            "status": "calculation_error", "error": str(error)})
            continue
        delta = updated-original
        valid = original > 0
        relative = np.full_like(original, np.nan)
        relative[valid] = delta[valid]/original[valid]
        index = int(np.nanargmax(np.abs(relative))) if valid.any() else None
        cooling.append({"registered_name": name, "ion": ion,
                        "affected_source_pairs": sum(key[0] == ion for key in PCHIP),
                        "status": "measured", "max_absolute_change_erg_cm3_s": float(np.max(np.abs(delta))),
                        "max_absolute_change_temperature_K": float(temperature[np.argmax(np.abs(delta))]),
                        "max_relative_change": float(relative[index]) if index is not None else None,
                        "max_relative_change_temperature_K": float(temperature[index]) if index is not None else None,
                        "baseline_at_relative_max": float(original[index]) if index is not None else None,
                        "updated_at_relative_max": float(updated[index]) if index is not None else None,
                        "nonzero_difference_samples": int(np.count_nonzero(delta)),
                        "baseline_zero_samples": int(np.count_nonzero(~valid))})
        print(name, cooling[-1]["max_relative_change"], flush=True)
    registered = {f"{element}_{stage}" for element, stage in INSTALLED_COOLING_IONS.values()}
    expanded = []
    affected_ions = sorted({key[0] for key in PCHIP} | {key[0] for key in NEGATIVE_NODE})
    for ion in affected_ions:
        if ion in registered:
            continue
        element, stage_text = ion.split("_")
        stage = int(stage_text)
        try:
            original, updated = cooling_pair(element, stage, temperature)
        except (ValueError, FileNotFoundError) as error:
            expanded.append({"ion": ion, "status": "calculation_error", "error": str(error)})
            continue
        delta = updated-original
        valid = original > 0
        relative = np.full_like(original, np.nan)
        relative[valid] = delta[valid]/original[valid]
        index = int(np.nanargmax(np.abs(relative))) if valid.any() else None
        expanded.append({"ion": ion, "status": "measured",
                         "affected_source_pairs": sum(key[0] == ion for key in PCHIP),
                         "max_absolute_change_erg_cm3_s": float(np.max(np.abs(delta))),
                         "max_absolute_change_temperature_K": float(temperature[np.argmax(np.abs(delta))]),
                         "max_relative_change": float(relative[index]) if index is not None else None,
                         "max_relative_change_temperature_K": float(temperature[index]) if index is not None else None,
                         "baseline_at_relative_max": float(original[index]) if index is not None else None,
                         "updated_at_relative_max": float(updated[index]) if index is not None else None,
                         "nonzero_difference_samples": int(np.count_nonzero(delta)),
                         "baseline_zero_samples": int(np.count_nonzero(~valid))})
        print("expanded", ion, expanded[-1]["max_relative_change"], flush=True)
    result = {"temperature_K": [float(temperature[0]), float(temperature[-1])],
              "temperature_points": len(temperature),
              "policy": "Previous CubicSpline plus zero clipping versus current source evaluator: PCHIP only for nonnegative-node CHIANTI pairs with negative exact full-source cubic extrema; preserve cubic type 6; reject negative strengths at requested T for negative-node source pairs",
              "installed_native_atoms": native, "registered_cooling": cooling,
              "additional_approved_stage_cooling": expanded}
    Path(__file__).with_name("scope.json").write_text(json.dumps(result, indent=2)+"\n")
    print("native_affected", sum(row["pchip_pair_count"] > 0 or row["negative_source_node_pair_count"] > 0 for row in native),
          "registered_measured", sum(row["status"] == "measured" for row in cooling))


if __name__ == "__main__":
    main()
