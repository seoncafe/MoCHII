#!/usr/bin/env python3
"""Check grouped HII20/HII40 outputs against the pre-change snapshot.

This is a migration gate, not a Cloudy-accuracy gate.  It protects the legacy
grouped path while continuous-energy infrastructure is added beside it.
"""

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent


def read_dataset(handle, name):
    return handle[f"{name}/data"][:]


def scalar_attr(group, name):
    value = np.asarray(group.attrs[name]).reshape(-1)[0]
    return float(value)


def close_enough(actual, expected, rtol, atol):
    return bool(np.isclose(actual, expected, rtol=rtol, atol=atol))


def report_value(label, actual, expected, rtol, atol):
    ok = close_enough(actual, expected, rtol, atol)
    print(
        f"  {label:30s} actual={actual: .12e} "
        f"reference={expected: .12e}  {'PASS' if ok else 'FAIL'}"
    )
    return ok


def reconstructed_edges(center, width):
    lower = (-width + np.sqrt(width * width + 4.0 * center * center)) / 2.0
    return lower, lower + width


def read_lines(path):
    lines = []
    with path.open() as stream:
        for raw in stream:
            fields = raw.split()
            if raw.startswith("#") or len(fields) < 5:
                continue
            lines.append(
                {
                    "element": fields[0].lower(),
                    "stage": int(fields[1]),
                    "wavelength_A": float(fields[2]),
                    "luminosity": float(fields[3]),
                    "ratio_Hbeta": float(fields[4]),
                }
            )
    return lines


def find_line(lines, expected):
    matches = [
        line
        for line in lines
        if line["element"] == expected["element"].lower()
        and line["stage"] == int(expected["stage"])
        and np.isclose(
            line["wavelength_A"],
            float(expected["wavelength_A"]),
            rtol=0.0,
            atol=0.01,
        )
    ]
    return matches[0] if len(matches) == 1 else None


def check_case(tag, spec, rtol, atol):
    path = REPO / spec["file"]
    if not path.exists():
        print(f"[{tag}] missing grouped rates file: {path}")
        return False

    ok = True
    print(f"=== {tag}: {path.relative_to(REPO)} ===")
    with h5py.File(path, "r") as handle:
        header = handle["Gamma_HI"]
        for name, expected in spec["attributes"].items():
            actual = scalar_attr(header, name)
            ok &= report_value(f"attribute {name}", actual, float(expected), rtol, atol)

        scalar_names = (
            "T_e",
            "n_e",
            "x_HI",
            "x_HeI",
            "x_HeII",
            "Gamma_HI",
            "Gamma_HeI",
            "Gamma_HeII",
            "Heat_HI",
            "Heat_HeI",
            "Heat_HeII",
        )
        arrays = {name: read_dataset(handle, name) for name in scalar_names}
        carbon = read_dataset(handle, "x_c_stages")
        xyz = read_dataset(handle, "LeafXYZ")
        radius = np.sqrt(np.sum(xyz * xyz, axis=0))

        for name in scalar_names:
            ok &= report_value(
                f"mean {name}",
                float(np.mean(arrays[name])),
                float(spec["means"][name]),
                rtol,
                atol,
            )
        for stage, expected in enumerate(spec["means"]["x_c_stages"], start=1):
            ok &= report_value(
                f"mean carbon stage {stage}",
                float(np.mean(carbon[:, stage - 1])),
                float(expected),
                rtol,
                atol,
            )

        half_width = 0.035
        for shell in spec["shells"]:
            rr = float(shell["r_pc"])
            mask = np.abs(radius - rr) < half_width
            nleaf = int(np.count_nonzero(mask))
            shell_label = f"r={rr:.2f} pc"
            if nleaf != int(shell["nleaf"]):
                print(
                    f"  {shell_label:30s} nleaf={nleaf} "
                    f"reference={shell['nleaf']}  FAIL"
                )
                ok = False
                continue
            values = {
                "x_CIII": float(np.median(carbon[mask, 2])),
                "x_He_ionized": float(np.median(1.0 - arrays["x_HeI"][mask])),
                "x_HI": float(np.median(arrays["x_HI"][mask])),
                "T_e": float(np.median(arrays["T_e"][mask])),
            }
            for name, actual in values.items():
                # Tiny grouped C III/He-ionized values outside the He front
                # are protected by an absolute floor rather than a meaningless
                # relative comparison of roundoff-sized numbers.
                local_atol = max(atol, 1.0e-12) if abs(float(shell[name])) < 1.0e-10 else atol
                ok &= report_value(
                    f"{shell_label} {name}",
                    actual,
                    float(shell[name]),
                    rtol,
                    local_atol,
                )

        energy = read_dataset(handle, "E_bin")
        width = read_dataset(handle, "dE_bin")
        lower, upper = reconstructed_edges(energy, width)
        carbon_bin = np.flatnonzero((lower < 24.383) & (upper > 24.383))
        if carbon_bin.size != 1:
            print(f"  C II threshold group count={carbon_bin.size}, expected 1  FAIL")
            ok = False
        else:
            i = int(carbon_bin[0])
            expected_triplet = (23.33, 23.950254904697783, 24.587)
            for name, actual, expected in zip(
                ("window lower", "window center", "window upper"),
                (lower[i], energy[i], upper[i]),
                expected_triplet,
            ):
                ok &= report_value(name, float(actual), expected, 1.0e-12, 1.0e-12)

    line_path = REPO / spec["line_file"]
    if not line_path.exists():
        print(f"  missing line-luminosity file: {line_path}  FAIL")
        ok = False
    else:
        actual_lines = read_lines(line_path)
        for expected in spec["lines"]:
            actual = find_line(actual_lines, expected)
            label = (
                f"line {expected['element']} {expected['stage']} "
                f"{expected['wavelength_A']:.2f} A"
            )
            if actual is None:
                print(f"  {label:30s} missing or ambiguous  FAIL")
                ok = False
                continue
            # The line file itself is formatted to six significant digits.
            ok &= report_value(
                f"{label} luminosity",
                actual["luminosity"],
                float(expected["luminosity"]),
                rtol,
                atol,
            )
            ok &= report_value(
                f"{label} / Hbeta",
                actual["ratio_Hbeta"],
                float(expected["ratio_Hbeta"]),
                rtol,
                atol,
            )

    print(f"{tag} GROUPED BASELINE: {'PASS' if ok else 'FAIL'}\n")
    return ok


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reference",
        type=Path,
        default=HERE / "grouped_baseline.json",
    )
    parser.add_argument("--rtol", type=float, default=5.0e-6)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    args = parser.parse_args()

    with args.reference.open() as stream:
        reference = json.load(stream)

    ok = True
    for tag, spec in reference["cases"].items():
        ok &= check_case(tag, spec, args.rtol, args.atol)

    print("OVERALL GROUPED BASELINE:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
