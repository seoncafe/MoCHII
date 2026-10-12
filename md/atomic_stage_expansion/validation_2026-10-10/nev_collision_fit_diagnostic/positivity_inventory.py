"""Inventory exact negative extrema of bound CHIANTI scaled-strength splines."""
from pathlib import Path
import hashlib
import json
import os
import sys

import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
from chianti_cooling import bound_excitation_atom, ion_dir
from fit_nlevel import bt_scale_t

ROSTER = {"c": 5, "n": 6, "o": 5, "ne": 5, "s": 7, "ar": 7,
          "mg": 4, "fe": 7, "si": 5, "cl": 7, "ca": 6}


def exact_minimum(xs, ys, domain=None):
    """Evaluate all cubic stationary points strictly within each knot interval."""
    spline = CubicSpline(xs, ys)
    c3, c2, c1, c0 = spline.c
    h = np.diff(xs)
    lo, hi = (float(xs[0]), float(xs[-1])) if domain is None else domain
    if hi < xs[0]:
        return float(ys[0]), float(xs[0])
    if lo > xs[-1]:
        return float(ys[-1]), float(xs[-1])
    lo, hi = max(float(xs[0]), lo), min(float(xs[-1]), hi)
    if lo > hi:
        raise ValueError("Invalid scaled-temperature interval")
    knots = np.r_[lo, hi, xs[(xs > lo) & (xs < hi)]]
    values = spline(knots)
    candidates = [(float(values[np.argmin(values)]),
                   float(knots[np.argmin(values)]))]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        qa, qb, qc = 3*c3, 2*c2, c1
        discriminant = qb*qb - 4*qa*qc
        sq = np.sqrt(np.maximum(discriminant, 0))
        r1 = (-qb - sq) / (2*qa)
        r2 = (-qb + sq) / (2*qa)
        linear = -qc/qb
    for roots in (r1, r2):
        absx = xs[:-1] + roots
        mask = (discriminant >= 0) & np.isfinite(roots) & (roots > 0) & (roots < h) & (absx >= lo) & (absx <= hi)
        dx = roots[mask]
        vals = ((c3[mask]*dx + c2[mask])*dx + c1[mask])*dx + c0[mask]
        if vals.size:
            j = int(np.argmin(vals))
            candidates.append((float(vals[j]), float(xs[:-1][mask][j] + dx[j])))
    mask = ((qa == 0) & np.isfinite(linear) & (linear > 0) & (linear < h)
            & (xs[:-1] + linear >= lo) & (xs[:-1] + linear <= hi))
    if mask.any():
        dx = linear[mask]
        vals = (c2[mask]*dx + c1[mask])*dx + c0[mask]
        j = int(np.argmin(vals))
        candidates.append((float(vals[j]), float(xs[:-1][mask][j] + dx[j])))
    return min(candidates)


def main():
    ions = []
    problems = []
    omitted = []
    total = {"bound_transitions": 0, "chianti_transitions": 0,
             "stout_transitions": 0,
             "type6_log_transitions": 0, "negative_source_nodes": 0,
             "negative_cubic_extrema_with_nonnegative_nodes": 0,
             "negative_cubic_in_fit_temperature_range": 0}
    for element, limit in ROSTER.items():
        for stage in range(1, limit + 1):
            label = f"{element}_{stage}"
            try:
                _, _, collisions = bound_excitation_atom(element, stage)
            except FileNotFoundError as error:
                omitted.append({"ion": label, "reason": "no excitation source",
                                "missing": os.path.basename(error.filename)})
                continue
            row = {"ion": label, "bound_transitions": len(collisions),
                   "chianti_transitions": 0, "stout_transitions": 0,
                   "type6_log_transitions": 0,
                   "negative_source_nodes": 0,
                   "negative_cubic_extrema_with_nonnegative_nodes": 0,
                   "negative_cubic_in_fit_temperature_range": 0}
            source = Path(ion_dir(element, stage)) / f"{label}.scups"
            if source.is_file():
                row["source_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
            for tr in collisions:
                if "temperature_grid" in tr:
                    row["stout_transitions"] += 1
                    continue
                if "xs" not in tr:
                    continue
                row["chianti_transitions"] += 1
                if tr["ttype"] == 6:
                    row["type6_log_transitions"] += 1
                    continue
                xs, ys = tr["xs"], tr["ys"]
                node_min = float(np.min(ys))
                minimum, min_x = exact_minimum(xs, ys)
                st = bt_scale_t(tr, np.array([1e3, 1e5]))
                fit_minimum, fit_x = exact_minimum(xs, ys,
                                                  (float(min(st)), float(max(st))))
                in_fit_range = fit_minimum < -max(1e-14, 1e-12*float(np.max(np.abs(ys))))
                if in_fit_range:
                    row["negative_cubic_in_fit_temperature_range"] += 1
                if node_min < 0:
                    row["negative_source_nodes"] += 1
                    problems.append({"ion": label, "pair": [tr["ll"], tr["ul"]],
                                     "type": tr["ttype"], "kind": "negative_source_node",
                                     "node_minimum": node_min,
                                     "cubic_minimum": minimum,
                                     "minimum_scaled_coordinate": min_x,
                                     "negative_in_fit_temperature_range": in_fit_range,
                                     "fit_range_minimum": fit_minimum,
                                     "fit_range_minimum_scaled_coordinate": fit_x})
                elif minimum < -max(1e-14, 1e-12*float(np.max(np.abs(ys)))):
                    row["negative_cubic_extrema_with_nonnegative_nodes"] += 1
                    problems.append({"ion": label, "pair": [tr["ll"], tr["ul"]],
                                     "type": tr["ttype"], "kind": "cubic_overshoot",
                                     "node_minimum": node_min,
                                     "cubic_minimum": minimum,
                                     "minimum_scaled_coordinate": min_x,
                                     "negative_in_fit_temperature_range": in_fit_range,
                                     "fit_range_minimum": fit_minimum,
                                     "fit_range_minimum_scaled_coordinate": fit_x,
                                     "source_C": tr["cups"], "de_Ry": tr["de"]})
            ions.append(row)
            for key in total:
                total[key] += row[key]
            print(label, row["chianti_transitions"],
                  row["negative_cubic_extrema_with_nonnegative_nodes"], flush=True)
    out = {"roster_stage_count": sum(ROSTER.values()), "summary": total,
           "ions": ions, "missing_sources": omitted, "affected_transitions": problems}
    Path(__file__).with_name("positivity_inventory.json").write_text(
        json.dumps(out, indent=2) + "\n")
    print("SUMMARY", total, "missing", len(omitted), flush=True)


if __name__ == "__main__":
    main()
