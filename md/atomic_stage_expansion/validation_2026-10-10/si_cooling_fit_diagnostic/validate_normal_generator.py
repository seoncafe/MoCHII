"""Check the rounded S I cooling file generated through the normal fitter."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C
import fit_cooling as F


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    generated = HERE / "normal_generator"
    product = generated / "cooling_s_1.txt"
    lines = [line.split() for line in product.read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    count = int(lines[0][0])
    coefficient = np.array([float(row[0]) for row in lines[1:]])
    excitation = np.array([float(row[1]) for row in lines[1:]])
    if count != coefficient.size or not np.all(np.isfinite(coefficient)):
        raise ValueError("Incomplete or nonfinite rounded coefficient table")
    metrics = {}
    for points in (801, 1601):
        temperature = np.logspace(3, 5, points)
        actual = F.evaluate(coefficient, excitation, temperature)
        source = C.cooling_low_density_limit("s", 1, temperature)
        error = np.abs(actual / source - 1)
        metrics[f"{points}_max_relative_error"] = float(np.max(error))
        metrics[f"{points}_worst_temperature_K"] = float(temperature[np.argmax(error)])
        if not np.all(np.isfinite(actual)) or np.any(actual <= 0):
            raise ValueError("Nonpositive fitted cooling in source interval")
    temperature = np.logspace(0, 9, 10001)
    extended = F.evaluate(coefficient, excitation, temperature)
    if not np.all(np.isfinite(extended)) or np.any(extended < 0):
        raise ValueError("Negative or nonfinite extended cooling")
    fit_temperature = np.logspace(3, 5, 801)
    terms = coefficient[None, :] * np.exp(
        -(excitation[None, :] - np.min(excitation)) / fit_temperature[:, None])
    condition = np.max(np.sum(np.abs(terms), axis=1) / np.sum(terms, axis=1))
    executable = HERE / "native_check/check_cooling_fit.x"
    native_output = generated / "native_801.txt"
    subprocess.run([str(executable), str(product), str(native_output)], check=True)
    native = np.loadtxt(native_output)
    if native.shape != (801, 2):
        raise ValueError("Incomplete linked Fortran output")
    native_source = C.cooling_low_density_limit("s", 1, native[:, 0])
    native_error = np.abs(native[:, 1] / native_source - 1)
    extended_executable = HERE / "native_check/check_si_extended.x"
    extended_output = generated / "native_extended.txt"
    subprocess.run([str(extended_executable), str(product), str(extended_output)], check=True)
    metrics.update({
        "ion": "s_1", "term_count": count,
        "zero_term_count": int(np.count_nonzero(coefficient == 0)),
        "condition_801_max": float(condition),
        "extended_min": float(extended.min()),
        "extended_max": float(extended.max()),
        "extended_nonnegative": bool(np.all(extended >= 0)),
        "native_801_max_relative_error": float(native_error.max()),
        "native_801_worst_temperature_K": float(native[np.argmax(native_error), 0]),
        "native_extended_output": extended_output.read_text().strip(),
        "product_sha256": sha256(product),
        "fitter_sha256": sha256(ROOT / "tools/fitting/fit_cooling.py"),
        "source_evaluator_sha256": sha256(ROOT / "tools/fitting/chianti_cooling.py"),
    })
    if max(metrics["801_max_relative_error"], metrics["1601_max_relative_error"],
           metrics["native_801_max_relative_error"]) >= 0.01 or condition >= 1e5:
        raise ValueError("S I cooling acceptance criterion failed")
    metrics["accepted_low_density_fit"] = True
    (generated / "validation.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
