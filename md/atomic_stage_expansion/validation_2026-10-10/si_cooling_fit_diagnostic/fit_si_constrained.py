"""Minimax S I cooling with global positivity and bounded cancellation."""
from pathlib import Path
import json

import numpy as np
from scipy.optimize import linprog

from fit_si_signed_minimax import T, NORMALIZED, SCALE, BASE, EXTENDED

KAPPA = 1e5


def candidate(count, step):
    offsets = np.arange(count, dtype=float)*step
    basis = np.exp(-offsets[None,:]/T[:,None])
    relative = basis/NORMALIZED[:,None]
    positives = np.exp(-offsets[None,:]/np.logspace(0,9,401)[:,None])
    ones = np.ones(len(T))
    zeros = np.zeros(len(T))
    constraints = np.vstack([
        np.column_stack([relative, -relative, -ones]),
        np.column_stack([-relative, relative, -ones]),
        np.column_stack([basis*(1/KAPPA-1), basis*(1/KAPPA+1), zeros]),
        np.column_stack([positives*(1/KAPPA-1),
                         positives*(1/KAPPA+1),
                         np.zeros(len(positives))])])
    rhs = np.r_[ones, -ones, zeros, np.zeros(len(positives))]
    result = linprog(np.r_[np.zeros(2*count), 1], A_ub=constraints,
                     b_ub=rhs, bounds=[(0,None)]*(2*count+1),
                     method="highs")
    if not result.success:
        return {"count":count,"step":step,"status":result.message}
    coeff = result.x[:count]-result.x[count:2*count]
    terms = basis*coeff[None,:]
    fitted = terms.sum(axis=1)/NORMALIZED
    err = np.max(np.abs(fitted-1))
    condition = np.max(np.sum(np.abs(terms),axis=1)/np.abs(terms.sum(axis=1)))
    extended = np.exp(-offsets[None,:]/EXTENDED[:,None])@coeff
    return {"count":count,"step":step,"status":"solved",
            "relative_error_801":float(err),"condition_801":float(condition),
            "negative_extended_points":int(np.count_nonzero(extended<0)),
            "minimum_extended_sum":float(extended.min()),
            "coefficient":(coeff*SCALE).tolist(),
            "excitation_temperature":(BASE+offsets).tolist()}


def main():
    rows=[candidate(count,step)
          for count in (18,24,32,40,48,56,64,80)
          for step in ((1000.,2000.,3000.,5000.,8000.) if count <=48
                       else (3000.,4000.,5000.,6000.,8000.))]
    rows.sort(key=lambda row:row.get("relative_error_801",1e9))
    output={"condition_limit":KAPPA,"rows":rows}
    Path(__file__).with_name("constrained_scan.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps([{k:v for k,v in row.items()
                       if k not in ("coefficient","excitation_temperature")}
                      for row in rows[:5]]))


if __name__=="__main__":main()
