#!/usr/bin/env python3
"""Validate continuous-energy HII20/HII40 production calculations.

The hard gates cover run integrity, energy closure, recovery of the C III
tail, and broad non-regression relative to the frozen grouped calculation.
Cloudy c25 profiles are an external diagnostic, not an exact-data identity
target: MoCHII and Cloudy do not use identical recombination and
charge-exchange data.
"""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CLOUDY = REPO / "tests" / "cloudy_c25"
PC_CM = 3.0856776e18
RIN_CM = 3.0e18
ROUT = {"hii20": 3.07874, "hii40": 4.73154}
TAIL_RADII = {
    "hii20": (1.5, 2.0, 2.5),
    "hii40": (4.2, 4.35, 4.5),
}
# A physically sampled tail should be many orders above the grouped
# threshold-bin collapse.  These floors remain deliberately far below Cloudy.
TAIL_FLOOR = {"hii20": 1.0e-5, "hii40": 1.0e-4}
PROFILE_FIELDS = {
    "x_HI": None,
    "x_HeI": None,
    "C II": ("x_c", 1),
    "C III": ("x_c", 2),
    "N II": ("x_n", 1),
    "N III": ("x_n", 2),
    "O II": ("x_o", 1),
    "O III": ("x_o", 2),
    "Ne II": ("x_ne", 1),
    "Ne III": ("x_ne", 2),
}


def dataset(handle, name):
    return handle[f"{name}/data"][:]


def attr(group, name):
    value = np.asarray(group.attrs[name]).reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode().strip()
    if isinstance(value, str):
        return value.strip()
    return value.item() if hasattr(value, "item") else value


def orient_stages(value, nleaf):
    if value.shape[0] == nleaf:
        return value
    if value.shape[1] == nleaf:
        return value.T
    raise ValueError(f"stage array {value.shape} is incompatible with {nleaf} leaves")


def read_mochii(path):
    with h5py.File(path, "r") as handle:
        header = handle["Gamma_HI"]
        xyz = dataset(handle, "LeafXYZ")
        nleaf = xyz.shape[1]
        data = {
            "r": np.sqrt(np.sum(xyz * xyz, axis=0)),
            "T_e": dataset(handle, "T_e"),
            "n_e": dataset(handle, "n_e"),
            "x_HI": dataset(handle, "x_HI"),
            "x_HeI": dataset(handle, "x_HeI"),
            "_attrs": {
                key: attr(header, key)
                for key in (
                    "CONVERGD", "FLDCONS", "NITERDON", "FINALDX", "FINALDTE",
                    "IONEMODE", "NUBINUSE", "ENRGSAMP", "L_EMIT", "L_ABS",
                    "L_ESC",
                )
                if key in header.attrs
            },
        }
        for element in ("c", "n", "o", "ne", "s"):
            data[f"x_{element}"] = orient_stages(
                dataset(handle, f"x_{element}_stages"), nleaf
            )
    return data


def read_cloudy_save(path):
    with path.open() as stream:
        header = stream.readline().lstrip("#").rstrip("\n").split("\t")
    values = np.loadtxt(path, skiprows=1)
    return {name.strip(): values[:, i] for i, name in enumerate(header)}


def stage_column(symbol, stage):
    return symbol if stage == 0 else symbol + "+" + ("" if stage == 1 else str(stage))


def read_cloudy(tag):
    overview = read_cloudy_save(CLOUDY / f"{tag}_c25.ovr")
    radius = (RIN_CM + overview["depth"]) / PC_CM
    keep = radius <= ROUT[tag]
    data = {"r": radius[keep], "T_e": overview["Te"][keep]}
    for element, symbol, nstage in (
        ("h", "H", 2), ("he", "He", 3), ("c", "C", 4),
        ("n", "N", 3), ("o", "O", 3), ("ne", "Ne", 3),
    ):
        table = read_cloudy_save(CLOUDY / f"{tag}_c25.ele_{element}")
        data[f"x_{element}"] = np.column_stack(
            [table[stage_column(symbol, stage)][keep] for stage in range(nstage)]
        )
    data["x_HI"] = data["x_h"][:, 0]
    data["x_HeI"] = data["x_he"][:, 0]
    return data


def shell_median(data, field, radius, half_width=0.035):
    mask = np.abs(data["r"] - radius) < half_width
    if not np.any(mask):
        return np.nan
    return float(np.median(field[mask]))


def radial_profile(data, field, tag, edges):
    values = np.full(edges.size - 1, np.nan)
    for index in range(values.size):
        mask = (data["r"] >= edges[index]) & (data["r"] < edges[index + 1])
        if np.any(mask):
            values[index] = np.median(field[mask])
    return values


def field(data, label):
    selector = PROFILE_FIELDS[label]
    return data[label] if selector is None else data[selector[0]][:, selector[1]]


def profile_error(model, cloudy, label, tag):
    edges = np.linspace(0.97223, ROUT[tag], 33)
    centers = 0.5 * (edges[:-1] + edges[1:])
    actual = radial_profile(model, field(model, label), tag, edges)
    reference = np.interp(centers, cloudy["r"], field(cloudy, label))
    valid = np.isfinite(actual) & np.isfinite(reference)
    return float(np.mean(np.abs(actual[valid] - reference[valid])))


def bool_attr(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return value != 0
    return str(value).strip().lower() in ("t", "true", "1", "yes")


def check_case(tag, continuous_path, grouped_path):
    print(f"=== {tag.upper()} ===")
    if not continuous_path.exists():
        print(f"missing continuous result: {continuous_path}\n")
        return False
    if not grouped_path.exists():
        print(f"missing grouped reference: {grouped_path}\n")
        return False

    continuous = read_mochii(continuous_path)
    grouped = read_mochii(grouped_path)
    cloudy = read_cloudy(tag)
    meta = continuous["_attrs"]
    ok = True

    checks = (
        ("converged", bool_attr(meta.get("CONVERGD", False))),
        ("final field consistency pass", bool_attr(meta.get("FLDCONS", False))),
        ("continuous energy mode", meta.get("IONEMODE") == "continuous"),
        ("diagnostic-only energy bins", meta.get("NUBINUSE") == "diagnostic_only"),
        ("Sobol launch sequence", meta.get("ENRGSAMP") == "sobol"),
    )
    for label, passed in checks:
        print(f"  {label:31s} {'PASS' if passed else 'FAIL'}")
        ok &= passed

    emitted = float(meta.get("L_EMIT", np.nan))
    absorbed = float(meta.get("L_ABS", np.nan))
    escaped = float(meta.get("L_ESC", np.nan))
    closure = abs(absorbed + escaped - emitted) / emitted
    closure_ok = np.isfinite(closure) and closure <= 5.0e-13
    print(f"  energy closure                 {closure:.3e}  "
          f"{'PASS' if closure_ok else 'FAIL'}")
    ok &= closure_ok
    print(f"  iterations/final dx/final dTe  {meta.get('NITERDON', '?')} / "
          f"{float(meta.get('FINALDX', np.nan)):.3e} / "
          f"{float(meta.get('FINALDTE', np.nan)):.3e}")

    print("  C III tail (shell median; Cloudy is interpolated):")
    tail_ok = True
    for radius in TAIL_RADII[tag]:
        c_value = shell_median(continuous, continuous["x_c"][:, 2], radius)
        g_value = shell_median(grouped, grouped["x_c"][:, 2], radius)
        cloudy_value = float(np.interp(radius, cloudy["r"], cloudy["x_c"][:, 2]))
        passed = np.isfinite(c_value) and c_value >= TAIL_FLOOR[tag]
        tail_ok &= passed
        print(
            f"    r={radius:4.2f} pc  continuous={c_value:10.3e}  "
            f"grouped={g_value:10.3e}  Cloudy={cloudy_value:10.3e}  "
            f"{'PASS' if passed else 'FAIL'}"
        )
    ok &= tail_ok

    # Broad non-regression gate: no H/He/metal radial-fraction MAE may worsen
    # against Cloudy by more than 25% plus a 0.02 absolute allowance.  This is
    # intentionally tolerant of atomic-data differences while catching a
    # transport regression in another species.
    print("  radial mean absolute error versus Cloudy:")
    profile_ok = True
    for label in PROFILE_FIELDS:
        c_error = profile_error(continuous, cloudy, label, tag)
        g_error = profile_error(grouped, cloudy, label, tag)
        passed = c_error <= 1.25 * g_error + 0.02
        profile_ok &= passed
        print(f"    {label:7s} continuous={c_error:8.4f}  grouped={g_error:8.4f}  "
              f"{'PASS' if passed else 'FAIL'}")
    ok &= profile_ok

    ionized_c = continuous["x_HI"] < 0.95
    ionized_g = grouped["x_HI"] < 0.95
    te_c = float(np.mean(continuous["T_e"][ionized_c]))
    te_g = float(np.mean(grouped["T_e"][ionized_g]))
    te_shift = abs(te_c / te_g - 1.0)
    te_ok = np.isfinite(te_shift) and te_shift <= 0.10
    print(f"  mean Te(x_HI<0.95)            {te_c:.1f} K vs grouped {te_g:.1f} K"
          f"  shift={te_shift:.2%}  {'PASS' if te_ok else 'FAIL'}")
    ok &= te_ok
    print(f"{tag.upper()} CONTINUOUS PRODUCTION: {'PASS' if ok else 'FAIL'}\n")
    return ok


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hii20", type=Path, default=REPO / "hii20_continuous_rates.h5")
    parser.add_argument("--hii40", type=Path, default=REPO / "hii40_continuous_rates.h5")
    args = parser.parse_args()
    ok20 = check_case(
        "hii20", args.hii20, REPO / "tests/g2_hii/hii20_bench_rates.h5"
    )
    ok40 = check_case(
        "hii40", args.hii40, REPO / "tests/g2_hii/hii40_bench_rates.h5"
    )
    print("CONTINUOUS HII PRODUCTION GATE:", "PASS" if ok20 and ok40 else "FAIL")
    return 0 if ok20 and ok40 else 1


if __name__ == "__main__":
    sys.exit(main())
