#!/usr/bin/env python3
"""Require full continuous HII production states to ignore diagnostic bins."""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
NBINS = (8, 16, 32, 64, 128)
EXCLUDED_DATASETS = {"J_nu", "E_bin", "dE_bin"}
REQUIRED_ATTRIBUTES = ("CONVERGD", "FLDCONS", "IONEMODE", "NUBINUSE", "ENRGSAMP")


def value(attribute):
    value = np.asarray(attribute).reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode().strip()
    return value.item() if hasattr(value, "item") else value


def as_bool(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in ("t", "true", "1", "yes")


def data(handle, name):
    return handle[f"{name}/data"][:]


def same_array(reference, candidate):
    return np.array_equal(reference, candidate, equal_nan=True)


def compare_case(tag, directory):
    paths = {nbin: directory / f"{tag}_continuous_n{nbin}_rates.h5" for nbin in NBINS}
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        print(f"[{tag}] missing production outputs:")
        print("\n".join(f"  {path}" for path in missing))
        return False

    ok = True
    with h5py.File(paths[32], "r") as reference:
        reference_names = set(reference.keys()) - EXCLUDED_DATASETS
        for nbin, path in paths.items():
            with h5py.File(path, "r") as candidate:
                header = candidate["Gamma_HI"]
                metadata_ok = (
                    int(value(header.attrs["NNU_ION"])) == nbin
                    and as_bool(value(header.attrs["CONVERGD"]))
                    and as_bool(value(header.attrs["FLDCONS"]))
                    and value(header.attrs["IONEMODE"]) == "continuous"
                    and value(header.attrs["NUBINUSE"]) == "diagnostic_only"
                    and value(header.attrs["ENRGSAMP"]) == "sobol"
                )
                print(f"[{tag} nnu_ion={nbin:3d}] metadata "
                      f"{'PASS' if metadata_ok else 'FAIL'}")
                ok &= metadata_ok
                if nbin == 32:
                    continue
                names_ok = set(candidate.keys()) - EXCLUDED_DATASETS == reference_names
                if not names_ok:
                    print(f"[{tag} nnu_ion={nbin:3d}] dataset set FAIL")
                    ok = False
                    continue
                for name in sorted(reference_names):
                    passed = same_array(data(reference, name), data(candidate, name))
                    if not passed:
                        print(f"[{tag} nnu_ion={nbin:3d}] {name} bitwise FAIL")
                    ok &= passed
                print(f"[{tag} nnu_ion={nbin:3d}] physical arrays "
                      f"{'PASS' if ok else 'FAIL'}")
    print(f"{tag.upper()} DIAGNOSTIC-BIN INDEPENDENCE: {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", type=Path, default=REPO)
    parser.add_argument("--case", choices=("hii20", "hii40", "both"), default="both")
    args = parser.parse_args()
    cases = ("hii20", "hii40") if args.case == "both" else (args.case,)
    ok = all(compare_case(tag, args.directory) for tag in cases)
    print("FULL HII DIAGNOSTIC-BIN GATE:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
