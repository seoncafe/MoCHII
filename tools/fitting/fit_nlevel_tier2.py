"""Tier-2 n-level atom data for the MoCHII generic n-level solver.

For each registry ion, extracts the lowest NLEV levels and produces one
data/atomic/nlevel_<ion>.txt file containing

  1. level energies and statistical weights (.elvlc, observed energies
     preferred, else theoretical),
  2. radiative A-values among those levels (.wgfa; duplicate (l,u) rows
     summed),
  3. per-transition effective collision strengths Upsilon(T) as low-order
     Chebyshev fits in Burgess-Tully (1992) scaled space (.scups).

The Chebyshev fit replaces the cubic-spline-through-table evaluation:
Fortran stays table-free and evaluates

    et = kT/dE ;  st = 1 - ln(C)/ln(et+C)   (type 1,4)
                  st = et/(et+C)            (type 2,3,5,6)
    clip st to [st_lo, st_hi] ;  u = 2(st-st_lo)/(st_hi-st_lo) - 1
    y  = sum_k c_k T_k(u)                   (Chebyshev series)
    Upsilon = y*ln(et+e)  [1] ;  y [2] ;  y/(et+1) [3] ;
              y*ln(et+C)  [4] ;  y/et [5] ;  10^y [6]

Fit domain: the Burgess-Tully image of T in [1e3, 1e5] K intersected with
the tabulated .scups range (outside it the origin spline is clipped anyway).
Fit error is measured on the descaled Upsilon(T) against the direct
CHIANTI spline evaluation (chianti_cooling.upsilon).

Writes:  ../../data/atomic/nlevel_<elem>_<stage>.txt
Run:     python3 fit_nlevel_tier2.py
"""

import os
import numpy as np
from numpy.polynomial import chebyshev as C
from chianti_cooling import (ion_dir, read_elvlc, read_wgfa, read_scups,
                             upsilon, DBASE, KB_OVER_RY)

DATE = "2026-07-11"
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "..", "data", "atomic")
TMIN, TMAX = 1.0e3, 1.0e5
TGRID = np.logspace(np.log10(TMIN), np.log10(TMAX), 201)
ERR_TARGET = 2.0e-3   # per-transition max relative error target
DEG_MAX = 12

# (element, stage, NLEV) — lowest levels needed for the optical/IR
# diagnostics; extra levels are harmless to the solver.
IONS = {
    "OII":   ("o", 2, 5),
    "OIII":  ("o", 3, 6),
    "NII":   ("n", 2, 5),
    "SII":   ("s", 2, 5),
    "SIII":  ("s", 3, 6),
    "NeII":  ("ne", 2, 2),
    "NeIII": ("ne", 3, 5),
    # MOCASSIN HII20/40 line set additions ([C II] 158um, C III] 1907/1909,
    # [N III] 57um, [O I] 6300/63um)
    "CII":   ("c", 2, 5),
    "CIII":  ("c", 3, 5),
    "NIII":  ("n", 3, 5),
    "OI":    ("o", 1, 5),
    # G5: argon ([Ar III] 7136/8.99um; [Ar IV] 4711/4740 density pair).
    # Ar II has no CHIANTI level data (stage tracked without lines).
    "ArIII": ("ar", 3, 5),
    "ArIV":  ("ar", 4, 5),
    # G5: Mg II resonance doublet; [Fe II] lowest terms; [Fe III] with the
    # 3P/3H optical multiplets (4658 etc.) needs 14 levels.
    "MgII":  ("mg", 2, 3),
    "FeII":  ("fe", 2, 13),
    "FeIII": ("fe", 3, 14),
    # Si, Cl, Ca ([Si II] 34.8um; Si III] 1883/1892; Si IV 1394/1403;
    # [Cl II] 8579/9124; [Cl III] 5518/5538 density pair; [Cl IV] 8046;
    # Ca II H&K + IR triplet; [Ca V] 5309/6087).
    "SiII":  ("si", 2, 5),
    "SiIII": ("si", 3, 5),
    "SiIV":  ("si", 4, 3),
    "ClII":  ("cl", 2, 5),
    "ClIII": ("cl", 3, 5),
    "ClIV":  ("cl", 4, 5),
    "CaII":  ("ca", 2, 5),
    "CaV":   ("ca", 5, 5),
}

SPEC_LABEL = {
    "OII": "[O II]", "OIII": "[O III]", "NII": "[N II]",
    "SII": "[S II]", "SIII": "[S III]", "NeII": "[Ne II]",
    "NeIII": "[Ne III]", "CII": "[C II]", "CIII": "C III]",
    "NIII": "[N III]", "OI": "[O I]",
    "ArIII": "[Ar III]", "ArIV": "[Ar IV]",
    "MgII": "Mg II", "FeII": "[Fe II]", "FeIII": "[Fe III]",
    "SiII": "[Si II]", "SiIII": "Si III]", "SiIV": "Si IV",
    "ClII": "[Cl II]", "ClIII": "[Cl III]", "ClIV": "[Cl IV]",
    "CaII": "Ca II", "CaV": "[Ca V]",
}

with open(os.path.join(DBASE, "VERSION")) as fh:
    CHIANTI_VERSION = fh.read().strip()


def bt_scale_t(tr, T):
    """Burgess-Tully scaled temperature st(T) for a .scups transition."""
    et = (KB_OVER_RY * np.asarray(T, float)) / tr["de"]
    Cc = tr["cups"]
    if tr["ttype"] in (1, 4):
        return 1.0 - np.log(Cc) / np.log(et + Cc)
    return et / (et + Cc)


def fit_transition(tr):
    """Chebyshev-fit the scaled Upsilon of one transition over TGRID.

    Returns (st_lo, st_hi, coeffs, logf, maxerr) with maxerr the relative
    error of the descaled Upsilon(T) against the direct spline evaluation.
    logf = 1 means the series stores log10 of the scaled Upsilon (used for
    weak peaked transitions whose table spans decades — a linear-space fit
    cannot hold the relative error there).
    """
    from scipy.interpolate import CubicSpline
    st_t = bt_scale_t(tr, np.array([TMIN, TMAX]))
    st_lo = max(min(st_t), tr["xs"][0])
    st_hi = min(max(st_t), tr["xs"][-1])
    ups_ref = upsilon(tr, TGRID)
    spl = CubicSpline(tr["xs"], tr["ys"])

    if st_hi <= st_lo:               # T range entirely outside the table
        st_lo, st_hi = tr["xs"][0], tr["xs"][-1]
        c0 = float(spl(np.clip(bt_scale_t(tr, TGRID), st_lo, st_hi)).mean())
        cf = np.array([c0])
        return st_lo, st_hi, cf, 0, _maxerr(tr, st_lo, st_hi, cf, 0, ups_ref)

    sdense = np.linspace(st_lo, st_hi, 400)
    ydense = spl(sdense)
    best = None
    for logf in (0, 1):
        if logf and not (ydense > 0.0).all():
            continue
        ytarget = np.log10(ydense) if logf else ydense
        for deg in range(2, DEG_MAX + 1):
            cf = C.Chebyshev.fit(sdense, ytarget, deg,
                                 domain=[st_lo, st_hi]).coef
            err = _maxerr(tr, st_lo, st_hi, cf, logf, ups_ref)
            if best is None or err < best[0]:
                best = (err, cf, logf)
            if err < ERR_TARGET:
                break
        if best[0] < ERR_TARGET:
            break
    return st_lo, st_hi, best[1], best[2], best[0]


def _descale(tr, T, y):
    """Descale a scaled-Upsilon value y at temperature T (BT 1992)."""
    et = (KB_OVER_RY * np.asarray(T, float)) / tr["de"]
    tt = tr["ttype"]
    if tt == 1:
        return y * np.log(et + np.e)
    if tt == 2:
        return y
    if tt == 3:
        return y / (et + 1.0)
    if tt == 4:
        return y * np.log(et + tr["cups"])
    if tt == 5:
        return y / et
    return 10.0 ** y


def eval_ups_fit(tr, st_lo, st_hi, coeffs, logf, T):
    """Evaluate the fitted Upsilon(T) exactly as the Fortran will."""
    st = np.clip(bt_scale_t(tr, T), st_lo, st_hi)
    u = 2.0 * (st - st_lo) / (st_hi - st_lo) - 1.0 if st_hi > st_lo else st * 0.0
    y = C.chebval(u, coeffs)
    if logf:
        y = 10.0 ** y
    return np.maximum(_descale(tr, T, y), 0.0)


def _maxerr(tr, st_lo, st_hi, coeffs, logf, ups_ref):
    fit = eval_ups_fit(tr, st_lo, st_hi, coeffs, logf, TGRID)
    mask = ups_ref > 0
    if not mask.any():
        return 0.0
    return np.abs(fit[mask] / ups_ref[mask] - 1.0).max()


def build_ion(name, elem, stage, nlev, suffix=""):
    d = ion_dir(elem, stage)
    lev_raw = read_elvlc(os.path.join(d, f"{elem}_{stage}.elvlc"))
    levels = {}
    for i in sorted(lev_raw)[:nlev]:
        g, e_obs, e_th = lev_raw[i]
        e = e_obs if (e_obs > 0.0 or i == 1) else e_th
        levels[i] = (g, e)

    arad = {}
    for ll, ul, wl, gf, a in read_wgfa(os.path.join(d, f"{elem}_{stage}.wgfa")):
        if ll in levels and ul in levels and ll < ul and a > 0.0:
            arad[(ll, ul)] = arad.get((ll, ul), 0.0) + a

    ups_rows, worst = [], 0.0
    for tr in read_scups(os.path.join(d, f"{elem}_{stage}.scups")):
        if tr["ll"] not in levels or tr["ul"] not in levels:
            continue
        st_lo, st_hi, cf, logf, err = fit_transition(tr)
        worst = max(worst, err)
        ups_rows.append((tr, st_lo, st_hi, cf, logf, err))

    path = os.path.join(OUTDIR, f"nlevel_{elem}_{stage}{suffix}.txt")
    with open(path, "w") as fh:
        fh.write(f"# MoCHII Tier-2 n-level atom data: {SPEC_LABEL[name]}"
                 f" ({elem}_{stage}), lowest {nlev} levels\n")
        fh.write(f"# source: CHIANTI v{CHIANTI_VERSION} .elvlc/.wgfa/.scups"
                 " (tools/fitting/fit_nlevel_tier2.py)\n")
        fh.write("# levels: index  g  E[cm^-1] (observed preferred)\n")
        fh.write("# rad: l u A[s^-1] (duplicate rows summed)\n")
        fh.write("# ups: l u type C dE[Ry] st_lo st_hi logf ncheb"
                 " c_0..c_{n-1} maxerr\n")
        fh.write("#   logf=1: the series stores log10 of the scaled Upsilon"
                 " (apply 10^y before descaling)\n")
        fh.write("#   Upsilon(T) = Burgess-Tully descaled Chebyshev series;"
                 " see fit_nlevel_tier2.py header for the evaluation recipe\n")
        fh.write(f"# fit range: {TMIN:.1e} - {TMAX:.1e} K;"
                 f" worst transition fit error: {worst*100:.2f}%\n")
        fh.write(f"# generated on {DATE}\n")
        fh.write(f"NLEV {len(levels)}\n")
        for i, (g, e) in levels.items():
            fh.write(f"  {i:3d}  {g:5.1f}  {e:14.4f}\n")
        fh.write(f"NRAD {len(arad)}\n")
        for (ll, ul), a in sorted(arad.items()):
            fh.write(f"  {ll:3d} {ul:3d}  {a:.6e}\n")
        fh.write(f"NUPS {len(ups_rows)}\n")
        for tr, st_lo, st_hi, cf, logf, err in ups_rows:
            cs = " ".join(f"{c:.8e}" for c in cf)
            fh.write(f"  {tr['ll']:3d} {tr['ul']:3d}  {tr['ttype']:1d}"
                     f" {tr['cups']:.6e} {tr['de']:.6e}"
                     f" {st_lo:.8e} {st_hi:.8e} {logf:1d} {len(cf):2d}  {cs}"
                     f"  {err*100:.3f}%\n")
    return path, len(levels), len(arad), len(ups_rows), worst


if __name__ == "__main__":
    print(f"{'ion':6s} {'nlev':>4s} {'nrad':>4s} {'nups':>4s}"
          f" {'worst_ups_err':>13s}")
    for name, (el, st, nlev) in IONS.items():
        path, nl, nr, nu, worst = build_ion(name, el, st, nlev)
        print(f"{name:6s} {nl:4d} {nr:4d} {nu:4d} {worst*100:12.3f}%"
              f"   -> {os.path.basename(path)}")

    # Optional expanded Fe II/III models (par%fe_levels_full): the compact
    # files above keep the solver small for the default run; these _full
    # files add the higher terms for iron-line studies.  Fe III to 34 levels
    # (the full 3d^6 optical multiplets), Fe II to 16.
    FE_FULL = {"FeII": ("fe", 2, 16), "FeIII": ("fe", 3, 34)}
    for name, (el, st, nlev) in FE_FULL.items():
        path, nl, nr, nu, worst = build_ion(name, el, st, nlev, suffix="_full")
        print(f"{name+'*':6s} {nl:4d} {nr:4d} {nu:4d} {worst*100:12.3f}%"
              f"   -> {os.path.basename(path)}")
