"""Validate a rounded existing-format S I signed-exponential candidate."""
from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as source
from fit_cooling import write_ion_file

HERE = Path(__file__).resolve().parent
T = np.logspace(3, 5, 801)


def series(amplitude, excitation, temperature):
    """Match the native summation order and operations."""
    result = np.zeros_like(temperature)
    for a, ti in zip(amplitude, excitation):
        result += a*np.exp(-ti/temperature)
    return result/np.sqrt(temperature)


def main():
    candidates = json.loads((HERE/"constrained_scan.json").read_text())["rows"]
    candidate = next(row for row in candidates if row["count"] == 64
                     and row["step"] == 5000.0 and row["status"] == "solved")
    amplitude = np.array([float(f"{a:.16e}") for a in candidate["coefficient"]])
    excitation = np.array([float(f"{ti:.16e}")
                           for ti in candidate["excitation_temperature"]])
    comparisons = []
    for count in (801,1601):
        temperature = np.logspace(3, 5, count)
        reference = source.cooling_low_density_limit("s", 1, temperature)
        fit = series(amplitude, excitation, temperature)
        error = fit/reference-1
        comparisons.append({"temperature_count": count,
                            "maximum_absolute_relative_error": float(np.abs(error).max()),
                            "worst_temperature_K": float(temperature[np.abs(error).argmax()]),
                            "minimum_fitted_cooling": float(fit.min()),
                            "maximum_fitted_cooling": float(fit.max())})
        if count == 801:
            fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7,7))
            axes[0].loglog(temperature, reference, label="corrected source")
            axes[0].loglog(temperature, fit, linestyle="--", label="64-term fit")
            axes[0].set_ylabel(r"Cooling $\Lambda$ [erg cm$^3$ s$^{-1}$]")
            axes[0].legend()
            axes[1].semilogx(temperature, 100*error)
            axes[1].set_ylabel("Relative difference [%]")
            axes[1].set_xlabel("Electron temperature [K]")
            fig.tight_layout()
            fig.savefig(HERE/"si_cooling_profile.pdf")
            plt.close(fig)
    extended = np.logspace(0, 9, 10001)
    values = series(amplitude, excitation, extended)
    offsets = excitation-excitation.min()
    scaled_terms = amplitude[None,:]*np.exp(-offsets[None,:]/T[:,None])
    condition = np.max(np.sum(np.abs(scaled_terms),axis=1)
                       /np.abs(np.sum(scaled_terms,axis=1)))
    assert np.all(np.isfinite(values)) and np.all(values >= 0)
    assert condition < 1e5
    assert max(row["maximum_absolute_relative_error"] for row in comparisons) < .01
    path = write_ion_file("SI", "s", 1, amplitude, excitation,
                          comparisons[0]["maximum_absolute_relative_error"],
                          comparisons[0]["maximum_absolute_relative_error"], HERE)
    result = {"file":str(Path(path).relative_to(ROOT)),
              "nterm":len(amplitude),"source": "current bound S I low-density cooling",
              "comparisons":comparisons,"condition_801":float(condition),
              "negative_extended_points":int(np.count_nonzero(values<0)),
              "minimum_extended_cooling":float(values.min()),
              "maximum_extended_cooling":float(values.max())}
    (HERE/"validation.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
