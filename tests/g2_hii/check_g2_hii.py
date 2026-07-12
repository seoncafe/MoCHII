#!/usr/bin/env python3
"""G2 gate: MoCHII vs the MOCASSIN HII20/HII40 baselines.

Compares, for each case:
  - ion-weighted mean temperatures Te(H I), Te(H II)
    (weight n_ion n_e V; MOCASSIN temperature.out rows (1,1), (1,2));
  - mean ionic fractions <X_i> (weight n_e V; MOCASSIN ionratio.out) for
    H, He and the registry metals (C, N, O, Ne, S).

Caveats recorded with the gate: MOCASSIN transports the diffuse field
explicitly (MoCHII G2 is case-B on-the-spot; G3 adds the diffuse field);
the two codes use different atomic-data compilations; and the stored
baselines are reduced runs (1e5 photons, 3 iterations, convergence flag
0.09).  Gate criterion: |Te(H II) - reference| / reference < 5%.
"""
import numpy as np
import h5py
import re, os

MBASE = os.path.expanduser(
    "~/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/tests/baselines")

CASES = [
    ("HII20", "hii20_rates.h5", f"{MBASE}/01_gas_HII20"),
    ("HII40", "hii40_rates.h5", f"{MBASE}/02_gas_HII40"),
]

ELEM_IDX = {"H": 1, "He": 2, "C": 6, "N": 7, "O": 8, "Ne": 10, "S": 16}


def read_mocassin_table(path):
    """Parse 'Element Ion value [value2]' rows -> {(elem, ion): [values]}."""
    out = {}
    with open(path) as fh:
        for ln in fh:
            toks = ln.split()
            if len(toks) >= 3 and toks[0].isdigit() and toks[1].isdigit():
                out[(int(toks[0]), int(toks[1]))] = [float(v) for v in toks[2:]]
    return out


for label, fname, mdir in CASES:
    try:
        f = h5py.File(fname, "r")
    except FileNotFoundError:
        print(f"[{label}] {fname} not found; skipped\n")
        continue
    te = f["T_e"]["data"][:]
    ne = f["n_e"]["data"][:]
    xhi = f["x_HI"]["data"][:]
    # nH from x/ne closure is awkward; gas cells have nH=100 by construction
    nH = np.where(ne > 0, 100.0, 0.0)
    gas = nH > 0
    nHI, nHII = nH*xhi, nH*(1.0 - xhi)

    mt = read_mocassin_table(f"{mdir}/temperature.out")
    mi = read_mocassin_table(f"{mdir}/ionratio.out")

    w1 = (nHI*ne)[gas]
    w2 = (nHII*ne)[gas]
    te_HI = np.sum(te[gas]*w1)/np.sum(w1)
    te_HII = np.sum(te[gas]*w2)/np.sum(w2)
    ref1 = mt[(1, 1)][0]
    ref2 = mt[(1, 2)][0]
    print(f"=== {label} ===")
    print(f"Te(H I)-weighted : MoCHII {te_HI:8.1f} K   MOCASSIN {ref1:8.1f} K"
          f"   dev {te_HI/ref1-1:+.2%}")
    print(f"Te(H II)-weighted: MoCHII {te_HII:8.1f} K   MOCASSIN {ref2:8.1f} K"
          f"   dev {te_HII/ref2-1:+.2%}")
    dev2 = te_HII/ref2 - 1.0
    gate = abs(dev2) < 0.05
    print(f"GATE (Te(H II) within 5% of MOCASSIN row (1,2)): "
          f"{'PASS' if gate else 'FAIL'} ({dev2:+.2%})")

    #--- mean ionic fractions, weight n_e V (MOCASSIN ionratio convention)
    wv = ne[gas]
    print(f"{'ion':10s} {'MoCHII':>10s} {'MOCASSIN':>10s}")
    xh = np.sum(xhi[gas]*wv)/np.sum(wv)
    print(f"{'H I':10s} {xh:10.4f} {mi[(1,1)][0]:10.4f}")
    print(f"{'H II':10s} {1-xh:10.4f} {mi[(1,2)][0]:10.4f}")
    xhe1 = f["x_HeI"]["data"][:];  xhe2 = f["x_HeII"]["data"][:]
    for nm, arr, key in (("He I", xhe1, (2, 1)), ("He II", xhe2, (2, 2))):
        v = np.sum(arr[gas]*wv)/np.sum(wv)
        print(f"{nm:10s} {v:10.4f} {mi[key][0]:10.4f}")
    for el in ("c", "n", "o", "ne", "s"):
        key = f"x_{el}_stages"
        if key not in f:
            continue
        st = f[key]["data"][:]          # (nleaf, nstage)
        Z = ELEM_IDX[el.capitalize() if el != "ne" else "Ne"]
        for i in range(st.shape[1]):
            v = np.sum(st[gas, i]*wv)/np.sum(wv)
            ref = mi.get((Z, i + 1), [np.nan])[0]
            print(f"{el.upper()+' '+('I'*(i+1) if i<3 else 'IV'):10s} "
                  f"{v:10.4f} {ref:10.4f}")
    print()
