#!/usr/bin/env python3
"""Check independent scrambled-Sobol HII production replicates."""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
SEEDS = (101, 202, 303, 404)
ROUT = {"hii20": 3.07874, "hii40": 4.73154}
RADII = {"hii20": (1.5, 2.0, 2.5), "hii40": (4.2, 4.35, 4.5)}


def array(handle, name):
    return handle[f"{name}/data"][:]


def attr(group, name):
    value = np.asarray(group.attrs[name]).reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode().strip()
    return value.item() if hasattr(value, "item") else value


def stages(value, nleaf):
    return value if value.shape[0] == nleaf else value.T


def read(path):
    with h5py.File(path, "r") as handle:
        header = handle["Gamma_HI"]
        xyz = array(handle, "LeafXYZ")
        nleaf = xyz.shape[1]
        return {
            "r": np.sqrt(np.sum(xyz * xyz, axis=0)),
            "x_HI": array(handle, "x_HI"),
            "T_e": array(handle, "T_e"),
            "x_c": stages(array(handle, "x_c_stages"), nleaf),
            "attrs": {name: attr(header, name) for name in
                      ("CONVERGD", "FLDCONS", "IONEMODE", "NUBINUSE",
                       "ENRGSAMP", "L_EMIT", "L_ABS", "L_ESC")},
        }


def shell_median(model, radius):
    mask = np.abs(model["r"] - radius) < 0.035
    return float(np.median(model["x_c"][mask, 2]))


def bool_value(value):
    return value is True or str(value).lower() in ("true", "t", "1")


def check_case(tag, directory):
    paths = [directory / f"{tag}_rqmc_{seed}_rates.h5" for seed in SEEDS]
    if any(not path.exists() for path in paths):
        print(f"[{tag}] missing one or more RQMC outputs")
        return False
    models = [read(path) for path in paths]
    ok = True
    for seed, model in zip(SEEDS, models):
        a = model["attrs"]
        closure = abs(float(a["L_ABS"]) + float(a["L_ESC"]) -
                      float(a["L_EMIT"])) / float(a["L_EMIT"])
        passed = (bool_value(a["CONVERGD"]) and bool_value(a["FLDCONS"]) and
                  a["IONEMODE"] == "continuous" and
                  a["NUBINUSE"] == "diagnostic_only" and
                  a["ENRGSAMP"] == "sobol" and np.isfinite(closure) and
                  closure <= 5.0e-13)
        print(f"[{tag} seed={seed}] metadata/closure {closure:.3e} "
              f"{'PASS' if passed else 'FAIL'}")
        ok &= passed

    te = np.array([np.mean(m["T_e"][m["x_HI"] < 0.95]) for m in models])
    ciii = np.array([[shell_median(m, radius) for radius in RADII[tag]]
                     for m in models])
    te_rel = float(np.std(te, ddof=1) / np.mean(te))
    ciii_rel = np.std(ciii, axis=0, ddof=1) / np.mean(ciii, axis=0)
    print(f"[{tag}] Te mean={np.mean(te):.2f} K, sample relative scatter={te_rel:.3%}")
    for radius, value in zip(RADII[tag], ciii_rel):
        print(f"[{tag}] C III r={radius:.2f} pc relative scatter={value:.3%}")
    # The ensemble is an uncertainty estimate, not a forced equality test.
    # These broad limits catch a failed scramble or a non-reproducible launch
    # stream while retaining a meaningful physical scatter measurement.
    scatter_ok = te_rel <= 0.05 and np.all(ciii_rel <= 0.10)
    print(f"[{tag}] ensemble scatter bounds {'PASS' if scatter_ok else 'FAIL'}")
    ok &= scatter_ok
    print(f"{tag.upper()} RQMC ENSEMBLE: {'PASS' if ok else 'FAIL'}\n")
    return ok


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", type=Path, default=REPO)
    parser.add_argument("--case", choices=("hii20", "hii40", "both"), default="both")
    args = parser.parse_args()
    cases = ("hii20", "hii40") if args.case == "both" else (args.case,)
    ok = all(check_case(tag, args.directory) for tag in cases)
    print("FULL RQMC ENSEMBLE GATE:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
