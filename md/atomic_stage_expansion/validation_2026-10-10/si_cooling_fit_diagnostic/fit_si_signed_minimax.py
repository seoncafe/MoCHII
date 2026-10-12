"""Explore existing-format signed exponential fits of corrected S I cooling."""
from pathlib import Path
import json
import sys

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as source

T = np.logspace(3, 5, 801)
LAMBDA = source.cooling_low_density_limit("s", 1, T)
BASE = 569.8995211793084
TARGET = LAMBDA*np.sqrt(T)*np.exp(BASE/T)
SCALE = float(np.max(TARGET))
NORMALIZED = TARGET/SCALE
EXTENDED = np.logspace(0, 9, 10001)


def candidate(count, low, high, scheme="geometric"):
    if scheme == "geometric":
        offsets = np.r_[0.0, np.geomspace(low, high, count-1)]
    else:
        offsets = np.arange(count, dtype=float)*low
    matrix = np.exp(-offsets[None,:]/T[:,None])/NORMALIZED[:,None]
    constraints = np.vstack([
        np.column_stack([matrix, -np.ones(len(T))]),
        np.column_stack([-matrix, -np.ones(len(T))])])
    result = linprog(np.r_[np.zeros(count), 1.0], A_ub=constraints,
                     b_ub=np.r_[np.ones(len(T)), -np.ones(len(T))],
                     bounds=[(None,None)]*count+[(0,None)], method="highs")
    if not result.success:
        return {"count": count, "low": low, "high": high, "scheme": scheme,
                "status": result.message}
    coefficient = result.x[:count]
    terms = np.exp(-offsets[None,:]/T[:,None])*coefficient[None,:]
    fitted = terms.sum(axis=1)/NORMALIZED
    relative = np.abs(fitted-1)
    condition = np.max(np.sum(np.abs(terms),axis=1)/np.abs(terms.sum(axis=1)))
    # Remove the positive base excitation before checking the sign across
    # the entire native extrapolation interval.
    extended_sum = np.exp(-offsets[None,:]/EXTENDED[:,None])@coefficient
    return {"count": count, "low": low, "high": high, "scheme": scheme, "status": "solved",
            "max_relative_error_801": float(relative.max()),
            "worst_temperature_K": float(T[relative.argmax()]),
            "condition_801": float(condition),
            "min_extended_scaled_sum": float(extended_sum.min()),
            "negative_extended_points": int(np.count_nonzero(extended_sum<0)),
            "coefficient": (coefficient*SCALE).tolist(),
            "excitation_temperature": (BASE+offsets).tolist()}


def main():
    rows = [candidate(count, low, high)
            for count in (14, 18, 22, 24)
            for low in (30.0, 100.0, 300.0, 1000.0)
            for high in (3e4, 1e5, 3e5)]
    rows.extend(candidate(count, step, step*(count-1), "uniform")
                for count in (14, 18, 22, 24, 32, 40)
                for step in (300., 600., 1000., 2000., 3000., 5000., 10000.))
    rows.sort(key=lambda row: row.get("max_relative_error_801", 1e9))
    summary = [{key:value for key,value in row.items()
                if key not in ("coefficient", "excitation_temperature")}
               for row in rows]
    acceptable = [row for row in rows if row["status"] == "solved"
                  and row["max_relative_error_801"] < .01
                  and row["condition_801"] < 1e5
                  and row["negative_extended_points"] == 0]
    selected = min(acceptable, key=lambda row:(row["count"], row["max_relative_error_801"])) if acceptable else None
    output = {"source": "current bound S I steady low-density source",
              "temperature_count": 801, "temperature_K": [1000,100000],
              "source_minmax": [float(LAMBDA.min()), float(LAMBDA.max())],
              "base_excitation_K": BASE,
              "candidate_summary": summary, "selected": selected}
    Path(__file__).with_name("fit_scan.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"best": summary[0],
                      "selected": None if selected is None else
                          {k:v for k,v in selected.items()
                           if k not in ("coefficient", "excitation_temperature")}}))


if __name__ == "__main__":
    main()
