"""Compare positive shape-preserving and cubic Ne V 6--15 interpolation."""
from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
from chianti_cooling import KB_OVER_RY, read_scups

SOURCE = Path.home() / "RT_Codes/CHIANTI/dbase/ne/ne_5/ne_5.scups"
TR = next(tr for tr in read_scups(SOURCE) if (tr["ll"], tr["ul"]) == (6, 15))
T = np.logspace(3, 5, 1601)
et = KB_OVER_RY*T/TR["de"]
source_s = et/(et+TR["cups"])
cubic_raw = CubicSpline(TR["xs"], TR["ys"])(source_s)
cubic_clipped = np.maximum(cubic_raw, 0)
pchip = PchipInterpolator(TR["xs"], TR["ys"])(source_s)
first_secant = (TR["ys"][1]-TR["ys"][0])/(TR["xs"][1]-TR["xs"][0])
secant_cubic = CubicSpline(TR["xs"], TR["ys"],
                          bc_type=((1, first_secant), "not-a-knot"))(source_s)


def minimax(c, log_strength=False, target_strength=None):
    target_strength = pchip if target_strength is None else target_strength
    s = et/(et+c)
    lo, hi = float(s[0]), float(s[-1])
    u = 2*(s-lo)/(hi-lo)-1
    matrix = cheb.chebvander(u, 15)
    if log_strength:
        target = np.log10(target_strength)
        relative = matrix
        upper = target
        lower = -target
        scale = 1
    else:
        scale = float(target_strength.max())
        scaled = target_strength/scale
        relative = matrix/scaled[:, None]
        upper = np.ones(len(T))
        lower = -upper
    constraints = np.vstack([np.column_stack([relative, -np.ones(len(T))]),
                             np.column_stack([-relative, -np.ones(len(T))])])
    result = linprog(np.r_[np.zeros(16), 1], A_ub=constraints,
                     b_ub=np.r_[upper, lower],
                     bounds=[(None, None)]*16+[(0, None)], method="highs")
    if not result.success:
        return {"C": c, "success": False, "message": result.message}
    coefficients = result.x[:16]*scale
    fitted = cheb.chebval(u, coefficients)
    if log_strength:
        fitted = 10**fitted
    error = np.abs(fitted/target_strength-1)
    return {"C": c, "success": True, "log_strength": log_strength,
            "max_relative_error_1601": float(error.max()),
            "worst_temperature_K": float(T[error.argmax()]),
            "coordinate_limits": [lo, hi],
            "coefficients": coefficients.tolist()}


def main():
    cases = [minimax(float(c), log_strength) for c in np.r_[np.logspace(-3, 2, 41), TR["cups"]]
             for log_strength in (False, True)]
    cases.sort(key=lambda row: row.get("max_relative_error_1601", 1e10))
    selected = cases[0]
    baseline = [row for row in cases if row["C"] == TR["cups"]]
    secant_cases = [minimax(float(c), log_strength, secant_cubic)
                    for c in np.r_[np.logspace(-3, 0, 31), TR["cups"]]
                    for log_strength in (False, True)]
    secant_cases.sort(key=lambda row: row.get("max_relative_error_1601", 1e10))
    delta = pchip-cubic_clipped
    relative_delta = np.divide(delta, pchip, out=np.zeros_like(delta), where=pchip>0)
    result = {
        "pair": [6, 15], "source_C": TR["cups"], "source_de_Ry": TR["de"],
        "source_type": TR["ttype"], "source_nodes": len(TR["xs"]),
        "temperatures_K": [float(T[0]), float(T[-1])], "temperature_count": len(T),
        "cubic_raw_minimum": float(cubic_raw.min()),
        "cubic_clipped_zero_count": int(np.count_nonzero(cubic_clipped == 0)),
        "pchip_minimum": float(pchip.min()),
        "pchip_maximum": float(pchip.max()),
        "first_secant_slope": float(first_secant),
        "secant_boundary_cubic_minimum": float(secant_cubic.min()),
        "secant_boundary_cubic_best_fit": secant_cases[0],
        "pchip_minus_cubic_max_absolute": float(np.abs(delta).max()),
        "pchip_minus_cubic_max_absolute_temperature_K": float(T[np.abs(delta).argmax()]),
        "pchip_minus_cubic_max_fraction_of_pchip": float(np.abs(relative_delta).max()),
        "source_C_fit": baseline,
        "coordinate_scan": [{k: v for k,v in row.items() if k != "coefficients"}
                            for row in cases], "best_fit": selected}
    Path(__file__).with_name("nev_pchip_profile.json").write_text(json.dumps(result, indent=2)+"\n")
    fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    axes[0].semilogx(T, cubic_raw, label="cubic, before clipping")
    axes[0].semilogx(T, cubic_clipped, label="cubic, clipped")
    axes[0].semilogx(T, pchip, label="PCHIP")
    axes[0].semilogx(T, secant_cubic, label="first-secant boundary cubic",
                     linestyle="--")
    axes[0].set_ylabel(r"Effective collision strength $\Upsilon$")
    axes[0].legend()
    axes[1].semilogx(T, relative_delta)
    axes[1].set_ylabel(r"$(\Upsilon_{\rm PCHIP}-\Upsilon_{\rm cubic})/\Upsilon_{\rm PCHIP}$")
    axes[1].set_xlabel("Electron temperature [K]")
    fig.tight_layout()
    fig.savefig(Path(__file__).with_name("nev_pchip_vs_cubic.pdf"))
    print(json.dumps({"profile": {k: v for k,v in result.items()
                                  if k.startswith("pchip_") or k.startswith("cubic_")},
                      "source_C_fit": [{k: v for k,v in row.items() if k != "coefficients"}
                                       for row in baseline],
                      "best_fit": {k: v for k,v in selected.items() if k != "coefficients"}}))
    print(json.dumps({"secant_boundary_cubic_minimum": result["secant_boundary_cubic_minimum"],
                      "secant_boundary_cubic_best_fit":
                          {k: v for k,v in secant_cases[0].items() if k != "coefficients"}}))


if __name__ == "__main__":
    main()
