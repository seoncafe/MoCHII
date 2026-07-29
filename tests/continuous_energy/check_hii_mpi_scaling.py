#!/usr/bin/env python3
"""Check rank agreement and report time/memory for the MPI smoke matrix."""

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
RANKS = (1, 2, 3, 5, 8, 68)
EXCLUDED = {"J_nu", "E_bin", "dE_bin", "LeafXYZ", "LeafSize"}
STATE_DATASETS = {
    "Gamma_HI", "Gamma_HeI", "Gamma_HeII", "Heat_HI", "Heat_HeI",
    "Heat_HeII", "x_HI", "x_HeI", "x_HeII", "n_e", "T_e",
    "x_c_stages", "Gamma_c_stages", "Heat_c_stages", "x_n_stages",
    "Gamma_n_stages", "Heat_n_stages", "x_o_stages", "Gamma_o_stages",
    "Heat_o_stages", "x_ne_stages", "Gamma_ne_stages", "Heat_ne_stages",
    "x_s_stages", "Gamma_s_stages", "Heat_s_stages",
}


def array(handle, name):
    return handle[f"{name}/data"][:]


def attr(group, name):
    value = np.asarray(group.attrs[name]).reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode().strip()
    return value.item() if hasattr(value, "item") else value


def read(path):
    with h5py.File(path, "r") as handle:
        header = handle["Gamma_HI"]
        values = {name: array(handle, name) for name in STATE_DATASETS}
        metadata = {name: attr(header, name) for name in
                    ("IONEMODE", "NUBINUSE", "ENRGSAMP", "L_EMIT", "L_ABS", "L_ESC")}
    return values, metadata


def check_case(tag, directory):
    paths = {rank: directory / f"{tag}_mpi_np{rank}_rates.h5" for rank in RANKS}
    if any(not path.exists() for path in paths.values()):
        print(f"[{tag}] missing MPI output")
        return False
    reference, ref_meta = read(paths[68])
    ok = True
    for rank, path in paths.items():
        values, metadata = read(path)
        closure = abs(float(metadata["L_ABS"]) + float(metadata["L_ESC"]) -
                      float(metadata["L_EMIT"])) / float(metadata["L_EMIT"])
        meta_ok = (metadata["IONEMODE"] == "continuous" and
                   metadata["NUBINUSE"] == "diagnostic_only" and
                   metadata["ENRGSAMP"] == "sobol" and
                   np.isfinite(closure) and closure <= 5.0e-12)
        rank_ok = True
        if rank != 68:
            for name in STATE_DATASETS:
                if not np.allclose(values[name], reference[name], rtol=5e-12, atol=5e-14,
                                   equal_nan=True):
                    rank_ok = False
                    print(f"[{tag} np={rank}] dataset mismatch: {name}")
        ok &= meta_ok and rank_ok
        time_path = directory / f"{tag}_mpi_np{rank}.time"
        timing = time_path.read_text().strip() if time_path.exists() else "missing timing"
        print(f"[{tag} np={rank}] closure={closure:.3e} "
              f"metadata={'PASS' if meta_ok else 'FAIL'} "
              f"rank_agreement={'PASS' if rank_ok else 'FAIL'} {timing}")
    print(f"{tag.upper()} MPI SCALING: {'PASS' if ok else 'FAIL'}\n")
    return ok


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", type=Path, default=REPO)
    parser.add_argument("--case", choices=("hii20", "hii40", "both"), default="both")
    args = parser.parse_args()
    cases = ("hii20", "hii40") if args.case == "both" else (args.case,)
    ok = all(check_case(tag, args.directory) for tag in cases)
    print("FULL MPI SCALING GATE:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
