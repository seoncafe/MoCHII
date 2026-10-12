#!/usr/bin/env python3
"""Measure N II and S II low-density fits and cooling profiles."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from check_fit_ratio import BASE, ROOT, DENSITY, digest, fit_value, maximum_relative


def main():
    ions = (("NII", "n", 2), ("SII", "s", 2))
    files = [Path(__file__), BASE / "n_s_result.json"]
    files.extend(BASE / label / "profiles.npz" for label, _, _ in ions)
    files.extend(ROOT / f"data/atomic/cooling_{element}_{stage}.txt"
                 for _, element, stage in ions)
    before = {str(path): digest(path) for path in files}
    results = {}
    for label, element, stage in ions:
        source = np.load(BASE / label / "profiles.npz")
        temperature = source["temperature_K"]
        fit = fit_value(ROOT / f"data/atomic/cooling_{element}_{stage}.txt",
                        temperature)
        full_ref = source["full_cooling_density_0"]
        installed_ref = source["installed_cooling_density_0"]
        if np.any(fit <= 0) or np.any(full_ref <= 0) or np.any(installed_ref <= 0):
            raise ValueError(f"{label}: nonpositive low-density cooling")
        figure, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
        figure.subplots_adjust(left=.16, right=.97, bottom=.11, top=.93,
                               wspace=.28, hspace=.28)
        density_results = []
        for j, density in enumerate(DENSITY):
            full = source[f"full_cooling_density_{j}"]
            five = source[f"five_cooling_density_{j}"]
            installed = source[f"installed_cooling_density_{j}"]
            partial = np.where(installed_ref <= fit,
                               installed + (fit - installed_ref),
                               (fit / installed_ref) * installed)
            full_ratio = fit * (full / full_ref)
            density_results.append({
                "density_cm3": density,
                "current_partial_vs_full_source": maximum_relative(partial, full, temperature),
                "full_ratio_vs_full_source": maximum_relative(full_ratio, full, temperature)})
            axis = axes.flat[j]
            axis.loglog(temperature, full, label="Full bound source")
            axis.loglog(temperature, five, label="Five-state source")
            axis.loglog(temperature, installed, label="Installed atom")
            axis.loglog(temperature, partial, label="Current partial rule")
            axis.loglog(temperature, full_ratio, label="Full density ratio")
            axis.set_xlabel(r"$T$ [K]")
            axis.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
            axis.set_title(r"$n_e = " + f"{density:g}" + r"\,{\rm cm}^{-3}$")
            axis.legend(fontsize=6)
        figure.savefig(BASE / label / "cooling_profiles_layout.pdf")
        plt.close(figure)
        results[label] = {"fit_vs_full_low_density": maximum_relative(fit, full_ref, temperature),
                          "density_results": density_results}
    after = {str(path): digest(path) for path in files}
    result = {"passed": before == after, "sources_unchanged": before == after,
              "source_sha256": before, "ions": results}
    (BASE / "n_s_fit_ratio_result_layout.json").write_text(json.dumps(result, indent=2) + "\n")
    if before != after:
        raise RuntimeError("Inputs changed during calculation")
    print(json.dumps({"passed": result["passed"], "ions": list(results)}))


if __name__ == "__main__":
    main()
