"""Tier-1 cooling fits for the first MoCHII HII-region ion set.

Fits the CHIANTI v11 optically-thin collisional line-cooling curves
Lambda(T) per (n_e * n_ion) with the multi-exponential form used in EXHALE
(fit_cno_formulas.py, same form as its Fe II fit):

    Lambda(T) = T^{-1/2} * sum_i A_i * exp(-T_i / T)   [erg cm^3 s^-1]

for the ions of docs/PLAN.md section 8 (first diagnostic set):
O II, O III, N II, S II, S III, Ne II, Ne III.

Population model: lower levels Boltzmann-distributed over the ground-term
fine structure (levels within 0.4 eV of ground); every collisional
excitation radiates (optically thin, low-density limit).  For the single
4S3/2 ground of O II and S II this reduces to the coronal limit.

CAVEAT (carried into each output header): the fits are low-density-limit
cooling.  Transitions with low critical density (e.g. [O II] 3726/3729,
[S II] 6717/6731, the fine-structure IR lines) saturate above n_crit; the
Tier-1 fits overestimate cooling there.  Density-dependent line ratios and
emissivities come from Tier 2 (Upsilon fits + n-level solve).

Writes:  ../../data/atomic/cooling_tier1_<ion>.txt (one file per ion,
         provenance header: CHIANTI version, population model, fit range,
         fit form, max fit error)
         cooling_tier1_fits.png (fit-vs-CHIANTI comparison)
Run:     python3 fit_cooling_tier1.py
"""

import os
import numpy as np
from scipy.optimize import least_squares
from chianti_cooling import cooling_effective, DBASE

DATE = "2026-07-11"
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "..", "data", "atomic")

with open(os.path.join(DBASE, "VERSION")) as fh:
    CHIANTI_VERSION = fh.read().strip()

T = np.logspace(3.0, 5.0, 201)
TMIN, TMAX = T[0], T[-1]

# (element, ion_stage, dE-ladder guesses [K] from the .elvlc level energies:
#  fine structure / forbidden optical / UV blocks)
IONS = {
    # H I: Ly-alpha (118348 K) dominated collisional-excitation cooling;
    # single 1s ground level, so 'ground_term' reduces to coronal.  Serves
    # the G2 H/He thermal balance through the same Tier-1 reader.
    "HI":    ("h", 1, [118000.0, 150000.0]),
    # G2 gate abundance set (MOCASSIN HII20/40) additions: C, N I/III, O I
    "CI":    ("c", 1, [40.0, 14700.0, 31000.0, 100000.0]),
    "CII":   ("c", 2, [90.0, 61900.0, 108000.0, 250000.0]),
    "CIII":  ("c", 3, [75300.0, 147000.0, 300000.0]),
    "NI":    ("n", 1, [27700.0, 41500.0, 120000.0]),
    "NIII":  ("n", 3, [250.0, 82300.0, 145000.0, 300000.0]),
    "OI":    ("o", 1, [300.0, 22800.0, 48600.0, 106000.0]),
    "OII":   ("o", 2, [38600.0, 58200.0, 172000.0]),
    "OIII":  ("o", 3, [300.0, 29200.0, 62100.0, 150000.0]),
    "NII":   ("n", 2, [130.0, 22000.0, 47000.0, 132000.0]),
    "SII":   ("s", 2, [21400.0, 35300.0, 114000.0]),
    "SIII":  ("s", 3, [800.0, 16300.0, 39100.0, 120000.0]),
    # Ne II is one fine-structure line (12.8 um) decaying ~T^-0.35 up to the
    # ~3.1e5 K UV block: the shallow decay needs a T_i ladder to superpose.
    "NeII":  ("ne", 2, [1120.0, 3000.0, 8000.0, 25000.0, 80000.0, 300000.0]),
    "NeIII": ("ne", 3, [1000.0, 37200.0, 80200.0, 250000.0]),
    # G5: argon.  Ar II has no CHIANTI v11 level/collision data (only
    # RR/DR), so the Ar II stage carries no cooling fit ([Ar II] 6.98um
    # omitted); Ar III/IV are complete.
    "ArIII": ("ar", 3, [1600.0, 5000.0, 20200.0, 47900.0, 150000.0]),
    "ArIV":  ("ar", 4, [30400.0, 50200.0, 171000.0]),
    # G5: magnesium and iron (gas-phase; strongly depleted onto grains).
    "MgII":  ("mg", 2, [51400.0, 101000.0, 200000.0]),
    "FeII":  ("fe", 2, [600.0, 3000.0, 12000.0, 20000.0, 60000.0]),
    "FeIII": ("fe", 3, [900.0, 29000.0, 36000.0, 100000.0]),
    # Si, Cl, Ca (5-stage registry; Si V / Cl V and the Ca III/IV gaps
    # carry no cooling — level data absent or excitation negligible at
    # nebular temperatures).
    "SiII":  ("si", 2, [290.0, 45000.0, 63000.0, 150000.0]),
    "SiIII": ("si", 3, [500.0, 76800.0, 110000.0, 148000.0, 300000.0]),
    "SiIV":  ("si", 4, [103000.0, 200000.0]),
    "ClII":  ("cl", 2, [1000.0, 1430.0, 16560.0, 34000.0, 100000.0]),
    "ClIII": ("cl", 3, [26100.0, 41600.0, 100000.0, 250000.0]),
    "ClIV":  ("cl", 4, [700.0, 1900.0, 20700.0, 43800.0, 130000.0]),
    "CaII":  ("ca", 2, [19700.0, 25400.0, 60000.0]),
    "CaV":   ("ca", 5, [2400.0, 5900.0, 27000.0, 55000.0, 150000.0]),
}

SPEC_LABEL = {  # for plot titles (ASCII / LaTeX-safe)
    "HI": "H I", "CI": "[C I]", "CII": "[C II]", "CIII": "C III]",
    "NI": "[N I]", "NIII": "[N III]", "OI": "[O I]",
    "OII": "[O II]", "OIII": "[O III]", "NII": "[N II]",
    "SII": "[S II]", "SIII": "[S III]", "NeII": "[Ne II]",
    "NeIII": "[Ne III]", "ArIII": "[Ar III]", "ArIV": "[Ar IV]",
    "MgII": "Mg II", "FeII": "[Fe II]", "FeIII": "[Fe III]",
    "SiII": "[Si II]", "SiIII": "Si III]", "SiIV": "Si IV",
    "ClII": "[Cl II]", "ClIII": "[Cl III]", "ClIV": "[Cl IV]",
    "CaII": "Ca II", "CaV": "[Ca V]",
}


def fit_multiexp(lam, nterm, Tguess):
    """Lambda = T^-1/2 sum A_i exp(-T_i/T); LM on log Lambda."""
    y = np.log(lam * np.sqrt(T))

    def resid(th):
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            A = np.exp(np.clip(th[:nterm], -700.0, 700.0))
            Ti = np.exp(np.clip(th[nterm:], -700.0, 700.0))
            r = np.log(np.sum(A[:, None] * np.exp(-Ti[:, None] / T[None, :]),
                              axis=0)) - y
        return np.nan_to_num(r, nan=1e3, posinf=1e3, neginf=-1e3)

    th0 = np.concatenate([np.full(nterm, np.log(lam[-1] * np.sqrt(T[-1]) / nterm)),
                          np.log(np.asarray(Tguess, float))])
    sol = least_squares(resid, th0, method="lm", max_nfev=60000)
    A = np.exp(np.clip(sol.x[:nterm], -700.0, 700.0))
    Ti = np.exp(np.clip(sol.x[nterm:], -700.0, 700.0))
    o = np.argsort(Ti)
    return A[o], Ti[o]


def evaluate(A, Ti, T):
    return np.sum(A[:, None] * np.exp(-Ti[:, None] / T[None, :]), axis=0) / np.sqrt(T)


def write_ion_file(name, elem, stage, A, Ti, err, rel_err):
    path = os.path.join(OUTDIR, f"cooling_tier1_{elem}_{stage}.txt")
    with open(path, "w") as fh:
        fh.write(f"# MoCHII Tier-1 cooling coefficients: {SPEC_LABEL[name]}"
                 f" ({elem}_{stage}) collisional line cooling\n")
        fh.write("# form: Lambda(T) = T^-1/2 * sum_i A_i * exp(-T_i/T)"
                 "   [erg cm^3 s^-1] per (n_e * n_ion)\n")
        fh.write(f"# source: CHIANTI v{CHIANTI_VERSION} .elvlc/.scups,"
                 " Burgess-Tully descaling (tools/fitting/chianti_cooling.py)\n")
        fh.write("# population: Boltzmann over ground-term fine structure"
                 " (levels within 0.4 eV); optically thin, low-density limit\n")
        fh.write("# caveat: no density suppression above n_crit;"
                 " density-dependent emissivities are Tier 2\n")
        fh.write(f"# fit range: {TMIN:.1e} - {TMAX:.1e} K;"
                 f" max fit error: {err*100:.2f}%"
                 f" ({rel_err*100:.2f}% where Lambda > 1e-3 max)\n")
        fh.write(f"# generated by tools/fitting/fit_cooling_tier1.py on {DATE}\n")
        fh.write(f"{len(A)}\n")
        for a, t in zip(A, Ti):
            fh.write(f"  {a:.8e}  {t:.6e}\n")
    return path


results = {}
print(f"{'ion':6s} {'nterm':>5s} {'maxerr_fit':>11s} {'maxerr_relevant':>16s}")
for name, (el, st, guess) in IONS.items():
    lam = cooling_effective(el, st, T, pop="ground_term")
    best = None
    for nterm in (len(guess), len(guess) + 1, len(guess) + 2):
        g = guess + [2.0e5] * (nterm - len(guess))
        try:
            A, Ti = fit_multiexp(lam, nterm, g)
        except Exception:
            continue
        fit = evaluate(A, Ti, T)
        err = np.abs(fit / lam - 1.0).max()
        if best is None or err < best[0]:
            best = (err, A, Ti, fit)
    err, A, Ti, fit = best
    mask = lam > 1e-3 * lam.max()
    rel_err = np.abs(fit / lam - 1.0)[mask].max()
    path = write_ion_file(name, el, st, A, Ti, err, rel_err)
    print(f"{name:6s} {len(A):5d} {err*100:10.2f}% {rel_err*100:15.2f}%"
          f"   -> {os.path.relpath(path, os.path.join(OUTDIR, '..', '..'))}")
    results[name] = dict(elem=el, stage=st, lam=lam, fit=fit, A=A, Ti=Ti)

# --- comparison figure -------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

nrow = (len(results) + 3)//4
fig, axes = plt.subplots(nrow, 4, figsize=(15, 3.5*nrow), sharex=True)
for ax, (name, r) in zip(axes.flat, results.items()):
    ax.loglog(T, r["lam"], "k-", lw=1.8, label="CHIANTI")
    ax.loglog(T, r["fit"], "r--", lw=1.4, label="fit")
    ax.set_title(SPEC_LABEL[name])
    ax.set_xlim(TMIN, TMAX)
    lo = max(r["lam"].max() * 1e-8, r["lam"][r["lam"] > 0].min())
    ax.set_ylim(lo, r["lam"].max() * 3)
    if ax is axes.flat[0]:
        ax.legend(frameon=False, fontsize=9)
for ax in axes[-1]:
    ax.set_xlabel(r"$T$ [K]")
for ax in axes[:, 0]:
    ax.set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
for k in range(len(results), axes.size):
    axes.flat[k].axis("off")
fig.tight_layout()
fig.savefig("cooling_tier1_fits.png", dpi=150)
print("wrote cooling_tier1_fits.png")
