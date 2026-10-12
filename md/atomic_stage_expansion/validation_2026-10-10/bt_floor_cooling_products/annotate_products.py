"""Copy accepted coefficients with their physical temperature limits."""

import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
SPECIES = {"CIV": ("c", 4, 3200, 32000000),
           "SiV": ("si", 5, 5000, 50000000)}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def numerical_rows(text):
    return [line for line in text.splitlines() if line.strip()
            and not line.lstrip().startswith("#")]


def main():
    candidate = HERE / "candidate"
    candidate.mkdir(exist_ok=False)
    results = {}
    for ion, (element, stage, minimum, maximum) in SPECIES.items():
        original_dir = HERE / ion
        source = original_dir / f"cooling_{element}_{stage}.txt"
        fit = json.loads((original_dir / f"{ion}.json").read_text())
        original = source.read_text()
        addition = [
            f"# source collision temperature interval: {minimum}--{maximum} K",
            f"# adopted cold approximation: constant source Upsilon below {minimum} K; actual T enters Maxwellian excitation and detailed balance",
            "# low-temperature source accuracy is unverified; this fit is a numerical representation over 1000--100000 K",
            "# no high-temperature source continuation above the stated maximum",
        ]
        addition.extend(f"# raw source SHA-256 {path}: {value}"
                        for path, value in sorted(fit["sources"].items()))
        lines = original.splitlines()
        first_data = next(index for index, line in enumerate(lines)
                          if line.strip() and not line.lstrip().startswith("#"))
        annotated = "\n".join(lines[:first_data] + addition + lines[first_data:]) + "\n"
        ion_dir = candidate / ion
        ion_dir.mkdir()
        product = ion_dir / source.name
        product.write_text(annotated)
        if numerical_rows(annotated) != numerical_rows(original):
            raise ValueError("Coefficient values changed")
        opened = ion_dir / "native_open.txt"
        subprocess.run([str(HERE / "build_mpi/check_cooling_fit.x"),
                        str(product), str(opened)], check=True, timeout=120)
        prior_native = original_dir / "native_reader/native_cooling.txt"
        if opened.read_bytes() != prior_native.read_bytes():
            raise ValueError("Header changed native cooling values")
        results[ion] = {
            "original_sha256": digest(source), "candidate_sha256": digest(product),
            "native_values_unchanged": True, "numeric_rows_unchanged": True,
            "maximum_rounded_fit_error": fit["max_relative_fit_error"],
            "source_underflow_points": fit["source_underflow_points"],
            "source_temperature_interval_K": [minimum, maximum],
            "cold_approximation": "constant Upsilon below source minimum; actual T in Maxwellian rates",
        }
    (candidate / "result.json").write_text(json.dumps({"passed": True, "ions": results}, indent=2) + "\n")


if __name__ == "__main__":
    main()
