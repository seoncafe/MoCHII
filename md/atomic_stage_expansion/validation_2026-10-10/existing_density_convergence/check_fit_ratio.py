#!/usr/bin/env python3
"""Measure the low-density fit and finite-density partial/full cooling rules."""

import hashlib
import json
import os
from pathlib import Path
import resource

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
BASE = Path(__file__).resolve().parent
IONS = (("CI", "c", 1), ("CII", "c", 2), ("NI", "n", 1))
DENSITY = [1e-6, 100., 1e6, 1e10]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fit_value(path, temperature):
    rows = [line for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    count = int(rows[0])
    if len(rows) != count + 1:
        raise ValueError(f"Malformed fit {path}")
    values = np.array([[float(v) for v in row.split()] for row in rows[1:]])
    return np.sum(values[:, 0, None] *
                  np.exp(-values[:, 1, None] / temperature), axis=0) / np.sqrt(temperature)


def maximum_relative(actual, expected, temperature):
    normal = expected > 64 * np.finfo(float).tiny
    if not np.any(normal) or np.any(actual[normal] <= 0):
        raise ValueError("No positive normal source cooling")
    relative = actual[normal] / expected[normal] - 1
    index = int(np.argmax(np.abs(relative)))
    return dict(maximum=float(abs(relative[index])),
                signed=float(relative[index]),
                temperature_K=float(temperature[np.flatnonzero(normal)[index]]))


def main():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    result_path = BASE / "result.json"
    initial = json.loads(result_path.read_text())
    files = [Path(__file__), result_path]
    for label, element, stage in IONS:
        files += [BASE / label / "profiles.npz",
                  ROOT / f"data/atomic/cooling_{element}_{stage}.txt"]
    before = {str(path): digest(path) for path in files}
    results = {}
    for label, element, stage in IONS:
        source = np.load(BASE / label / "profiles.npz")
        temperature = source["temperature_K"]
        fit_path = ROOT / f"data/atomic/cooling_{element}_{stage}.txt"
        fit = fit_value(fit_path, temperature)
        reference = source["full_cooling_density_0"]
        installed_reference = source["installed_cooling_density_0"]
        if np.any(fit <= 0) or np.any(reference <= 0) or np.any(installed_reference <= 0):
            raise ValueError(f"{label}: nonpositive low-density cooling")
        fit_error = maximum_relative(fit, reference, temperature)
        fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
        differences, difference_axes = plt.subplots(2, 2, figsize=(9, 7),
                                                  sharex=True)
        for figure in (fig, differences):
            figure.subplots_adjust(left=.16, right=.97, bottom=.11, top=.93,
                                   wspace=.28, hspace=.28)
        density_results = []
        for j, density in enumerate(DENSITY):
            full = source[f"full_cooling_density_{j}"]
            five = source[f"five_cooling_density_{j}"]
            installed = source[f"installed_cooling_density_{j}"]
            partial = np.where(installed_reference <= fit,
                               installed + (fit - installed_reference),
                               (fit / installed_reference) * installed)
            full_ratio = fit * (full / reference)
            if np.any(partial < 0) or np.any(full_ratio <= 0):
                raise ValueError(f"{label}: invalid combined cooling")
            density_results.append({
                "density_cm3": density,
                "current_partial_vs_full_source": maximum_relative(partial, full, temperature),
                "full_ratio_vs_full_source": maximum_relative(full_ratio, full, temperature),
                "installed_vs_full_source": maximum_relative(installed, full, temperature),
                "five_source_vs_full_source": maximum_relative(five, full, temperature)})
            saved = {"temperature_K": temperature, "full": full,
                     "five": five, "installed": installed, "partial": partial,
                     "full_density_ratio": full_ratio}
            np.savez_compressed(BASE / label / f"combined_density_{j}.npz", **saved)
            title = r"$n_e = " + f"{density:g}" + r"\,{\rm cm}^{-3}$"
            axis = axes.flat[j]
            axis.loglog(temperature, full, label="Full bound source")
            axis.loglog(temperature, five, label="Five-state source")
            axis.loglog(temperature, installed, label="Installed atom")
            axis.loglog(temperature, partial, label="Current partial rule")
            axis.loglog(temperature, full_ratio, label="Full density ratio")
            axis.set_title(title)
            axis.set_xlabel(r"$T$ [K]")
            axis.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
            axis.legend(fontsize=6)
            axis = difference_axes.flat[j]
            axis.semilogx(temperature, 100 * (five / full - 1),
                          label="Five-state source")
            axis.semilogx(temperature, 100 * (installed / full - 1),
                          label="Installed atom")
            axis.semilogx(temperature, 100 * (partial / full - 1),
                          label="Current partial rule")
            axis.semilogx(temperature, 100 * (full_ratio / full - 1),
                          label="Full density ratio")
            axis.set_title(title)
            axis.set_xlabel(r"$T$ [K]")
            axis.set_ylabel("Difference [%]")
            axis.legend(fontsize=6)
        fig.savefig(BASE / label / "cooling_profiles_layout.pdf")
        differences.savefig(BASE / label / "cooling_differences_layout.pdf")
        plt.close(fig)
        plt.close(differences)
        results[label] = {"full_levels": initial["ions"][label]["full_levels"],
                          "fit_vs_full_low_density": fit_error,
                          "density_results": density_results}
    after = {str(path): digest(path) for path in files}
    result = {"passed": before == after, "sources_unchanged": before == after,
              "source_sha256": before, "ions": results,
              "scope": "Source convergence and fit ratio; native full atom not yet validated"}
    (BASE / "fit_ratio_result_layout.json").write_text(json.dumps(result, indent=2) + "\n")
    if before != after:
        raise RuntimeError("Source changed during the calculation")
    print(json.dumps({"passed": result["passed"], "ions": list(results)}))


if __name__ == "__main__":
    main()
