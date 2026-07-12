"""End-to-end verification of the Tier-2 n-level data files against PyNeb.

Reads data/atomic/nlevel_<ion>.txt (the fitted product, NOT the CHIANTI
database), solves the n-level statistical equilibrium exactly as the future
Fortran nlevel_mod will, and compares the classic density/temperature
diagnostic ratios against PyNeb:

    [O III] 5007/4363     (Te diagnostic)
    [N II]  6584/5755     (Te diagnostic)
    [O II]  3726/3729     (n_e diagnostic)
    [S II]  6717/6731     (n_e diagnostic)

PyNeb ships its own atomic data (not CHIANTI v11), so this is a
cross-calibration against an independent reference, not a fit-residual
test: agreement at the few-percent level validates the whole chain
(file format -> Upsilon evaluation -> equilibrium solve -> emissivity).

Run:  python3 verify_nlevel_pyneb.py
"""

import os
import numpy as np

K_B_ERG = 1.380649e-16
RY_ERG = 2.1798723611035e-11
KB_OVER_RY = K_B_ERG / RY_ERG
COLL_PREF = 8.629e-6
HC_ERG_CM = 1.9864458571489287e-16   # h*c [erg cm]

ATOMDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "..", "data", "atomic")


def read_nlevel(path):
    """Parse a nlevel_<ion>.txt file into (levels, arad, ups) dicts."""
    levels, arad, ups = {}, {}, []
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    i = 0
    while i < len(lines):
        toks = lines[i].split()
        if not toks:
            i += 1
            continue
        key, n = toks[0], int(toks[1])
        if key == "NLEV":
            for k in range(1, n + 1):
                t = lines[i + k].split()
                levels[int(t[0])] = (float(t[1]), float(t[2]))
        elif key == "NRAD":
            for k in range(1, n + 1):
                t = lines[i + k].split()
                arad[(int(t[0]), int(t[1]))] = float(t[2])
        elif key == "NUPS":
            for k in range(1, n + 1):
                t = lines[i + k].split()
                nc = int(t[8])
                ups.append(dict(ll=int(t[0]), ul=int(t[1]), ttype=int(t[2]),
                                cups=float(t[3]), de=float(t[4]),
                                st_lo=float(t[5]), st_hi=float(t[6]),
                                logf=int(t[7]),
                                coef=np.array(t[9:9 + nc], dtype=float)))
        i += n + 1
    return levels, arad, ups


def ups_eval(tr, T):
    """Evaluate the fitted Upsilon(T) (the Fortran evaluation recipe)."""
    from numpy.polynomial import chebyshev as C
    T = np.asarray(T, dtype=float)
    et = (KB_OVER_RY * T) / tr["de"]
    Cc = tr["cups"]
    if tr["ttype"] in (1, 4):
        st = 1.0 - np.log(Cc) / np.log(et + Cc)
    else:
        st = et / (et + Cc)
    st = np.clip(st, tr["st_lo"], tr["st_hi"])
    if tr["st_hi"] > tr["st_lo"]:
        u = 2.0 * (st - tr["st_lo"]) / (tr["st_hi"] - tr["st_lo"]) - 1.0
    else:
        u = st * 0.0
    y = C.chebval(u, tr["coef"])
    if tr["logf"]:
        y = 10.0 ** y
    tt = tr["ttype"]
    if tt == 1:
        out = y * np.log(et + np.e)
    elif tt == 2:
        out = y
    elif tt == 3:
        out = y / (et + 1.0)
    elif tt == 4:
        out = y * np.log(et + Cc)
    elif tt == 5:
        out = y / et
    else:
        out = 10.0 ** y
    return np.maximum(out, 0.0)


def solve_pops(levels, arad, ups, T, ne):
    """Level populations (normalized) from statistical equilibrium."""
    idx = sorted(levels)
    nl = len(idx)
    pos = {i: k for k, i in enumerate(idx)}
    R = np.zeros((nl, nl))            # R[j,i]: rate i -> j  [s^-1]
    for tr in ups:
        l, u = tr["ll"], tr["ul"]
        gl, el = levels[l]
        gu, eu = levels[u]
        de_erg = (eu - el) * HC_ERG_CM
        q_ul = COLL_PREF / (gu * np.sqrt(T)) * ups_eval(tr, T)
        q_lu = (gu / gl) * q_ul * np.exp(-de_erg / (K_B_ERG * T))
        R[pos[u], pos[l]] += ne * q_lu
        R[pos[l], pos[u]] += ne * q_ul
    for (l, u), a in arad.items():
        R[pos[l], pos[u]] += a
    M = R - np.diag(R.sum(axis=0))    # dn/dt = M n = 0
    M[0, :] = 1.0                     # closure: sum n = 1
    b = np.zeros(nl)
    b[0] = 1.0
    n = np.linalg.solve(M, b)
    return {i: n[pos[i]] for i in idx}


def emissivity(levels, arad, pops, u, l):
    """Line emissivity n_u A_ul dE [erg s^-1 per ion] (no 4pi, no n_e)."""
    de_erg = (levels[u][1] - levels[l][1]) * HC_ERG_CM
    return pops[u] * arad[(l, u)] * de_erg


# (ion file, PyNeb atom args, ratio label, (u,l) numerator, (u,l) denominator,
#  PyNeb wavelengths num/den)
CASES = [
    ("nlevel_o_3.txt", ("O", 3), "[O III] 5007/4363", (4, 3), (5, 4), 5007, 4363),
    ("nlevel_n_2.txt", ("N", 2), "[N II] 6584/5755", (4, 3), (5, 4), 6584, 5755),
    ("nlevel_o_2.txt", ("O", 2), "[O II] 3726/3729", (2, 1), (3, 1), 3726, 3729),
    ("nlevel_s_2.txt", ("S", 2), "[S II] 6717/6731", (3, 1), (2, 1), 6717, 6731),
    # G5: argon
    ("nlevel_ar_3.txt", ("Ar", 3), "[Ar III] 7136/7751", (4, 1), (4, 2), 7136, 7751),
    ("nlevel_ar_4.txt", ("Ar", 4), "[Ar IV] 4740/4711", (2, 1), (3, 1), 4740, 4711),
    # G5: iron
    ("nlevel_fe_3.txt", ("Fe", 3), "[Fe III] 4658/4702", (10, 1), (12, 2), 4658, 4701),
    # Si, Cl, Ca additions
    ("nlevel_si_3.txt", ("Si", 3), "Si III] 1892/1883", (3, 1), (4, 1), 1892, 1883),
    ("nlevel_cl_3.txt", ("Cl", 3), "[Cl III] 5518/5538", (3, 1), (2, 1), 5518, 5538),
    ("nlevel_cl_2.txt", ("Cl", 2), "[Cl II] 8579/9124", (4, 1), (4, 2), 8579, 9124),
]

TE_GRID = [5000.0, 10000.0, 15000.0, 20000.0]
NE_GRID = [1.0e2, 1.0e3, 1.0e4]

if __name__ == "__main__":
    import pyneb as pn
    print(f"PyNeb {pn.__version__} (its own atomic data, not CHIANTI)\n")
    for fname, (elem, spec), label, (nu_u, nu_l), (de_u, de_l), wnum, wden in CASES:
        levels, arad, ups = read_nlevel(os.path.join(ATOMDIR, fname))
        atom = pn.Atom(elem, spec)
        src = pn.atomicData.getDataFile(f"{elem}{spec}")
        print(f"=== {label}  (MoCHII fit file: {fname}) ===")
        print(f"    PyNeb data: {src}")
        print(f"{'Te[K]':>8} {'ne[cm-3]':>9} {'MoCHII':>10} {'PyNeb':>10}"
              f" {'diff':>7}")
        worst = 0.0
        for Te in TE_GRID:
            for ne in NE_GRID:
                pops = solve_pops(levels, arad, ups, Te, ne)
                r_mo = (emissivity(levels, arad, pops, nu_u, nu_l)
                        / emissivity(levels, arad, pops, de_u, de_l))
                r_pn = (atom.getEmissivity(Te, ne, wave=wnum)
                        / atom.getEmissivity(Te, ne, wave=wden))
                d = r_mo / r_pn - 1.0
                worst = max(worst, abs(d))
                print(f"{Te:8.0f} {ne:9.0e} {r_mo:10.4g} {r_pn:10.4g}"
                      f" {d*100:+6.1f}%")
        print(f"--- worst |diff| = {worst*100:.1f}%\n")
