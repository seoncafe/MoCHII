#!/usr/bin/env python3
"""Summarize the HII40 residual experiments against the Cloudy reference.

The reported fronts use volume-weighted 0.01-pc spherical shells and linear
interpolation between shell centres.  This deliberately matches the diagnostic
used for the Figure 14 residual investigation rather than individual AMR-cell
values.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np


PC = 3.085677581491367e18
RIN = 0.97223
SHELL_EDGES = np.arange(0.98, 4.73 + 0.01, 0.01)
SHELL_CENTRES = 0.5 * (SHELL_EDGES[:-1] + SHELL_EDGES[1:])


def read_mochii(path):
    with h5py.File(path, "r") as handle:
        xyz = handle["LeafXYZ/data"][:]
        size = handle["LeafSize/data"][:]
        radius = np.sqrt(np.sum(xyz * xyz, axis=0))
        volume = (2.0 * size) ** 3
        fields = {
            "H0": handle["x_HI/data"][:],
            "He0": handle["x_HeI/data"][:],
            "C2": handle["x_c_stages/data"][:, 2],
        }
    profiles = {}
    index = np.digitize(radius, SHELL_EDGES) - 1
    for name, values in fields.items():
        result = np.full(SHELL_CENTRES.size, np.nan)
        for ibin in range(result.size):
            select = index == ibin
            if np.any(select):
                result[ibin] = np.average(values[select], weights=volume[select])
        profiles[name] = result
    return profiles


def read_cloudy(prefix):
    profiles = {}
    # ``genfromtxt(names=True)`` sanitises the Cloudy header ``C+2`` to C2.
    for suffix, column, name in (("ele_h", "H", "H0"),
                                 ("ele_he", "He", "He0"),
                                 ("ele_c", "C2", "C2")):
        table = np.genfromtxt(str(prefix) + "." + suffix, names=True)
        radius = (3.0e18 + table["depth"]) / PC
        profiles[name] = np.interp(SHELL_CENTRES, radius, table[column],
                                   left=np.nan, right=np.nan)
    return profiles


def outward_crossing(profile, level):
    valid = np.isfinite(profile)
    r, x = SHELL_CENTRES[valid], profile[valid]
    # H and He neutral fractions increase outward, whereas C2 declines.
    # A sign change handles both conventions without redefining the front.
    hits = np.where((x[:-1] - level) * (x[1:] - level) <= 0.0)[0]
    hits = hits[x[hits] != x[hits + 1]]
    if hits.size == 0:
        return np.nan
    i = hits[np.argmin(np.abs(x[hits] - level) + np.abs(x[hits + 1] - level))]
    return r[i] + (level - x[i]) * (r[i + 1] - r[i]) / (x[i + 1] - x[i])


def at_radius(profile, radius):
    return np.interp(radius, SHELL_CENTRES, profile)


def summarise(label, profiles):
    return [label, outward_crossing(profiles["He0"], 0.5),
            outward_crossing(profiles["C2"], 0.5),
            outward_crossing(profiles["H0"], 0.1),
            at_radius(profiles["C2"], 4.20), at_radius(profiles["C2"], 4.35),
            at_radius(profiles["C2"], 4.50)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path,
                        default=Path("results/continuous_energy/residual_investigation"))
    parser.add_argument("--production", type=Path,
                        default=Path("results/continuous_energy/production/hii40_continuous_rates.h5"))
    parser.add_argument("--cloudy-prefix", type=Path,
                        default=Path("tests/cloudy_c25/hii40_c25"))
    args = parser.parse_args()

    cases = [("MoCHII production", args.production)]
    for tag in ("diffuse_off", "no_hei", "mocrec", "l7", "threshold", "milne"):
        path = args.results / f"hii40_{tag}_rates.h5"
        if path.exists():
            cases.append((tag, path))

    rows = [summarise(label, read_mochii(path)) for label, path in cases]
    rows.append(summarise("Cloudy c25", read_cloudy(args.cloudy_prefix)))
    headings = ("case", "r(He0=0.5)", "r(C2=0.5)", "r(H0=0.1)",
                "C2@4.20", "C2@4.35", "C2@4.50")
    print(" | ".join(headings))
    print(" | ".join("---" for _ in headings))
    for row in rows:
        print(f"{row[0]} | " + " | ".join(f"{item:.5f}" for item in row[1:]))


if __name__ == "__main__":
    main()
