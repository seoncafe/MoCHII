"""Read-only diagnostic for Ne V 6--15 effective-strength storage."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, linprog

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "tools/fitting"))
from chianti_cooling import KB_OVER_RY, read_scups, upsilon
from fit_nlevel import TMIN, TMAX, fit_transition, eval_ups_fit

SOURCE = Path.home() / "RT_Codes/CHIANTI/dbase/ne/ne_5/ne_5.scups"
T = np.logspace(np.log10(TMIN), np.log10(TMAX), 1601)
transition = next(row for row in read_scups(SOURCE)
                  if (row["ll"], row["ul"]) == (6, 15))
reference = upsilon(transition, T)
source_spline = CubicSpline(transition["xs"], transition["ys"])


def source_unclipped(temperature):
    et = KB_OVER_RY * temperature / transition["de"]
    return float(source_spline(et / (et + transition["cups"])))


ROOT_T = brentq(source_unclipped, 3000.0, 4000.0)
ROOT_S = KB_OVER_RY * ROOT_T / transition["de"]
ROOT_S = ROOT_S / (ROOT_S + transition["cups"])
ROOT_SLOPE = float(source_spline(ROOT_S, 1)) * (
    KB_OVER_RY / transition["de"] * transition["cups"] /
    (KB_OVER_RY * ROOT_T / transition["de"] + transition["cups"]) ** 2)


def fitted(c, degree=15, anchor_root=False):
    """Minimize maximum relative error where positive; retain the zero region."""
    et = KB_OVER_RY * T / transition["de"]
    scaled = et / (et + c)
    left, right = scaled[0], scaled[-1]
    coordinate = 2 * (scaled - left) / (right - left) - 1
    matrix = cheb.chebvander(coordinate, degree)
    positive = reference > 0
    scale = reference.max()
    target = reference[positive] / scale
    relative_matrix = matrix[positive] / target[:, None]
    constraints = np.vstack([
        np.column_stack([relative_matrix, -np.ones(len(target))]),
        np.column_stack([-relative_matrix, -np.ones(len(target))]),
        np.column_stack([matrix[~positive], np.zeros((~positive).sum())]),
    ])
    if anchor_root:
        et_root = KB_OVER_RY * ROOT_T / transition["de"]
        st_root = et_root / (et_root + c)
        u_root = 2 * (st_root - left) / (right - left) - 1
        basis_root = cheb.chebvander(np.array([u_root]), degree)[0]
        du_dt = (2 / (right - left) * KB_OVER_RY / transition["de"]
                 * c / (et_root + c) ** 2)
        basis_slope = np.array([
            cheb.chebval(u_root, cheb.chebder(np.eye(degree + 1)[k])) * du_dt
            for k in range(degree + 1)])
        equalities = np.vstack([np.r_[basis_root, 0], np.r_[basis_slope, 0]])
        equality_targets = [0, ROOT_SLOPE / scale]
    else:
        equalities = equality_targets = None
    solution = linprog(np.r_[np.zeros(degree + 1), 1.0], A_ub=constraints,
                       b_ub=np.r_[np.ones(len(target)), -np.ones(len(target)),
                                  np.zeros((~positive).sum())],
                       A_eq=equalities, b_eq=equality_targets,
                       bounds=[(None, None)] * (degree + 1) + [(0, None)],
                       method="highs")
    if not solution.success:
        return {"C": c, "success": False, "message": solution.message}
    coeff = solution.x[:-1] * scale
    curve = np.maximum(cheb.chebval(coordinate, coeff), 0)
    error = np.abs(curve[positive] / reference[positive] - 1)
    positive_T = T[positive]
    return {"C": c, "success": True, "degree": degree, "anchor_root": anchor_root,
            "log_strength": False,
            "max_relative_error": float(error.max()),
            "worst_temperature_K": float(positive_T[error.argmax()]),
            "zero_region_false_positives": int(np.count_nonzero(curve[~positive])),
            "coefficients": coeff.tolist(),
            "coordinate_limits": [float(left), float(right)]}


def main():
    baseline = fit_transition(transition)
    positive = reference > 0
    original = np.abs(eval_ups_fit(transition, *baseline[:4], T)[positive] / reference[positive] - 1)
    candidates = [fitted(float(c), 15)
                  for c in np.logspace(-4, 2, 31)]
    anchored_candidates = [fitted(float(c), 15, True)
                           for c in np.logspace(-4, 2, 31)]
    candidates = sorted(candidates, key=lambda row: row.get("max_relative_error", 1e9))
    anchored_candidates = sorted(anchored_candidates,
                                 key=lambda row: row.get("max_relative_error", 1e9))
    selected = fitted(0.04, 15, True)
    selected["C"] = float(f"{selected['C']:.6e}")
    selected["coordinate_limits"] = [float(f"{v:.8e}")
                                     for v in selected["coordinate_limits"]]
    selected["coefficients"] = [float(f"{v:.8e}")
                                for v in selected["coefficients"]]
    et_root = KB_OVER_RY * ROOT_T / transition["de"]
    st_root = et_root / (et_root + selected["C"])
    lo, hi = selected["coordinate_limits"]
    u_root = 2 * (st_root - lo) / (hi - lo) - 1
    coefficients = np.array(selected["coefficients"])
    coefficients[0] = float(f"{coefficients[0] - cheb.chebval(u_root, coefficients):.8e}")
    selected["coefficients"] = coefficients.tolist()
    dense_T = np.logspace(np.log10(TMIN), np.log10(TMAX), 160001)
    dense_reference = upsilon(transition, dense_T)
    et = KB_OVER_RY * dense_T / transition["de"]
    st = et / (et + selected["C"])
    dense_fit = np.maximum(cheb.chebval(2 * (st - lo) / (hi - lo) - 1,
                                        coefficients), 0)
    positive_dense = dense_reference > 0
    dense_error = np.abs(dense_fit[positive_dense] / dense_reference[positive_dense] - 1)
    selected["max_relative_error_160001"] = float(dense_error.max())
    selected["worst_temperature_160001_K"] = float(dense_T[positive_dense][dense_error.argmax()])
    selected["zero_region_false_positives_160001"] = int(np.count_nonzero(dense_fit[~positive_dense]))
    selected["rounded_root_residual"] = float(cheb.chebval(u_root, coefficients))
    selected["native_row"] = (f"  6  15  2 {selected['C']:.6e} {transition['de']:.6e} "
                              f"{lo:.8e} {hi:.8e} 0 16  "
                              + " ".join(f"{v:.8e}" for v in coefficients)
                              + f"  {selected['max_relative_error_160001']*100:.3f}%")
    result = {"source": str(SOURCE), "pair": [6, 15],
              "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "source_de_Ry": transition["de"], "source_C": transition["cups"],
              "source_type": transition["ttype"], "source_nodes": transition["xs"].tolist(),
              "source_scaled_strengths": transition["ys"].tolist(),
              "reference_minmax": [float(reference.min()), float(reference.max())],
              "zero_temperature_range_K": [float(T[~positive].min()), float(T[~positive].max())],
              "positive_temperature_start_K": float(T[positive].min()),
              "source_zero_crossing_K": ROOT_T,
              "source_zero_crossing_slope_per_K": ROOT_SLOPE,
              "original_fit": {"max_error_201": baseline[4],
                               "max_error_1601": float(original.max()),
                               "worst_T_1601": float(T[positive][original.argmax()]),
                               "coefficient_count": len(baseline[2])},
              "coordinate_candidates": [{key: value for key, value in row.items()
                                         if key != "coefficients"} for row in candidates],
              "best": candidates[0],
              "anchored_candidates": [{key: value for key, value in row.items()
                                       if key != "coefficients"} for row in anchored_candidates],
              "best_anchored": anchored_candidates[0],
              "selected_native_storage": selected}
    destination = Path(__file__).with_name("result.json")
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"original": result["original_fit"],
                      "best": {key: val for key, val in candidates[0].items()
                               if key != "coefficients"},
                      "best_anchored": {key: val for key, val in anchored_candidates[0].items()
                                        if key != "coefficients"},
                      "selected": {key: val for key, val in selected.items()
                                   if key not in ("coefficients", "native_row")}}))


if __name__ == "__main__":
    main()
