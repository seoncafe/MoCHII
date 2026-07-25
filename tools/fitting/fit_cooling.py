"""Cooling fits for the MoCHII HII-region ion set.

Fits the CHIANTI v11 optically-thin collisional line-cooling curves
Lambda(T) per (n_e * n_ion) with the multi-exponential form used in EXHALE
(fit_cno_formulas.py, same form as its Fe II fit):

    Lambda(T) = T^{-1/2} * sum_i A_i * exp(-T_i / T)   [erg cm^3 s^-1]

for the diagnostic ion set: O II, O III, N II, S II, S III, Ne II, Ne III,
and the rest of the registry.

Population model: the n_e -> 0 limit of the n-level statistical equilibrium
(chianti_cooling.cooling_low_density_limit).  At n_e << n_crit every ion
sits in its lowest level, so collisional excitation happens out of level 1
only and the ensuing radiative cascade radiates the full dE_{1u}; the total
cooling is the excitation power out of the ground level.  Transition
energies are the observed .elvlc level differences, the same energies the
Tier-2 n-level solve uses in both the Boltzmann factor and the emitted
photon (the .scups dE is the Burgess-Tully scaling energy, not a physical
transition energy).

CAVEAT (carried into each output header): the fits are low-density-limit
cooling.  Transitions with low critical density (e.g. [O II] 3726/3729,
[S II] 6717/6731, the fine-structure IR lines, and especially [C II] 158um
whose n_crit ~ 50 cm^-3) saturate above n_crit; the cooling fits
overestimate cooling there.  Density-dependent line ratios and emissivities
come from the n-level solve (Upsilon fits + statistical equilibrium).

Writes:  ../../data/atomic/cooling_<ion>.txt (one file per ion,
         provenance header: CHIANTI version, population model, fit range,
         fit form, max fit error)
         cooling_fits.png (fit-vs-CHIANTI comparison)
Run:     python3 fit_cooling.py
"""

import os
import numpy as np
from scipy.optimize import least_squares
from chianti_cooling import (cooling_low_density_limit, ion_dir, read_elvlc,
                             read_scups, upsilon, _level_energy_cm1,
                             transition_energy_erg, DBASE, K_B_ERG, COLL_PREF)

DATE = "2026-07-25"
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "..", "data", "atomic")

with open(os.path.join(DBASE, "VERSION")) as fh:
    CHIANTI_VERSION = fh.read().strip()

T = np.logspace(3.0, 5.0, 201)
TMIN, TMAX = T[0], T[-1]

NTERM_MAX = 8         # most exponential terms allowed in one fit
ERR_TARGET = 1.0e-2   # stop adding terms once the max fit error is below this
NSTART = 6            # restarts per term count (perturbed ladders, fixed seed)

# (element, ion_stage).  The T_i of the fit are excitation temperatures and
# are derived from the level energies (ground_excitation_ladder), not tuned.
IONS = {
    # H I: Ly-alpha (118348 K) dominated collisional-excitation cooling.
    # Serves the H/He thermal balance through the same cooling reader.
    "HI":    ("h", 1),
    # gate abundance set (MOCASSIN HII20/40) additions: C, N I/III, O I
    "CI":    ("c", 1),
    "CII":   ("c", 2),
    "CIII":  ("c", 3),
    "NI":    ("n", 1),
    "NIII":  ("n", 3),
    "OI":    ("o", 1),
    "OII":   ("o", 2),
    "OIII":  ("o", 3),
    "NII":   ("n", 2),
    "SII":   ("s", 2),
    "SIII":  ("s", 3),
    # S IV: the [S IV] 10.51um ground fine-structure line (951 cm^-1 = 1369 K)
    # plus the 3s3p^2 UV blocks; a ~1.2% coolant of the HII40 benchmark.
    "SIV":   ("s", 4),
    "NeII":  ("ne", 2),
    "NeIII": ("ne", 3),
    # G5: argon.  Ar II has no CHIANTI v11 level/collision data (only
    # RR/DR), so the Ar II stage carries no cooling fit ([Ar II] 6.98um
    # omitted); Ar III/IV are complete.
    "ArIII": ("ar", 3),
    "ArIV":  ("ar", 4),
    # G5: magnesium and iron (gas-phase; strongly depleted onto grains).
    "MgII":  ("mg", 2),
    "FeII":  ("fe", 2),
    "FeIII": ("fe", 3),
    # Si, Cl, Ca (5-stage registry; Si V / Cl V and the Ca III/IV gaps
    # carry no cooling — level data absent or excitation negligible at
    # nebular temperatures).
    "SiII":  ("si", 2),
    "SiIII": ("si", 3),
    "SiIV":  ("si", 4),
    "ClII":  ("cl", 2),
    "ClIII": ("cl", 3),
    "ClIV":  ("cl", 4),
    "CaII":  ("ca", 2),
    "CaV":   ("ca", 5),
}

SPEC_LABEL = {  # for plot titles (ASCII / LaTeX-safe)
    "HI": "H I", "CI": "[C I]", "CII": "[C II]", "CIII": "C III]",
    "NI": "[N I]", "NIII": "[N III]", "OI": "[O I]",
    "OII": "[O II]", "OIII": "[O III]", "NII": "[N II]",
    "SII": "[S II]", "SIII": "[S III]", "SIV": "[S IV]", "NeII": "[Ne II]",
    "NeIII": "[Ne III]", "ArIII": "[Ar III]", "ArIV": "[Ar IV]",
    "MgII": "Mg II", "FeII": "[Fe II]", "FeIII": "[Fe III]",
    "SiII": "[Si II]", "SiIII": "Si III]", "SiIV": "Si IV",
    "ClII": "[Cl II]", "ClIII": "[Cl III]", "ClIV": "[Cl IV]",
    "CaII": "Ca II", "CaV": "[Ca V]",
}


def ground_excitation_ladder(elem, stage, nterm):
    """Starting T_i ladder [K] built from the excitation temperatures.

    In the low-density limit the cooling sum is

        Lambda(T) = T^-1/2 sum_u [COLL_PREF/g_1 * Ups_1u(T) * dE_1u]
                                 * exp(-T_u/T) ,   T_u = dE_1u/k ,

    which is the fitted form exactly, apart from the slow temperature
    dependence of Upsilon.  So the T_i are the excitation temperatures of the
    upper levels: group the transitions out of the ground level into nterm
    clusters in log T_u (weighted k-means, the weight being the share of the
    cooling each transition carries across the fit range) and return the
    weighted mean T_u of each cluster.  This replaces hand-tuned starting
    ladders and follows the atomic structure of each ion automatically.
    """
    d = ion_dir(elem, stage)
    lev = _level_energy_cm1(read_elvlc(os.path.join(d, f"{elem}_{stage}.elvlc")))
    Tu, contrib = [], []
    for tr in read_scups(os.path.join(d, f"{elem}_{stage}.scups")):
        if tr["ll"] != 1:
            continue
        de_erg = transition_energy_erg(lev, tr)
        if de_erg <= 0.0:
            continue
        lam_u = (COLL_PREF / (lev[1][0] * np.sqrt(T)) * upsilon(tr, T)
                 * np.exp(-de_erg / (K_B_ERG * T)) * de_erg)
        Tu.append(de_erg / K_B_ERG)
        contrib.append(lam_u)
    contrib = np.asarray(contrib)
    tot = contrib.sum(axis=0)
    wt = (contrib / np.where(tot > 0.0, tot, 1.0)).sum(axis=1)
    Tu = np.asarray(Tu)
    keep = wt > 1.0e-8 * wt.max()
    Tu, wt = Tu[keep], wt[keep]

    x = np.log10(Tu)
    ncl = min(nterm, len(np.unique(np.round(x, 3))))
    # weighted quantiles of the cooling share as the cluster seeds
    o = np.argsort(x)
    cw = np.cumsum(wt[o]) / wt.sum()
    cen = np.interp((np.arange(ncl) + 0.5) / ncl, cw, x[o])
    for _ in range(50):
        lab = np.argmin(np.abs(x[:, None] - cen[None, :]), axis=1)
        new = np.array([(x[lab == k] * wt[lab == k]).sum() / wt[lab == k].sum()
                        if (lab == k).any() else cen[k] for k in range(ncl)])
        if np.allclose(new, cen, rtol=1e-6):
            break
        cen = new
    # Ions with only a few excitation blocks (e.g. the Cl ions have four
    # transitions out of the ground term) still need extra terms: a single
    # exponential has a fixed A, so the temperature dependence of Upsilon has
    # to be carried by superposing neighbouring terms.  Pad by splitting the
    # widest gaps, then by bracketing the ladder.
    cen = list(np.sort(cen))
    while len(cen) < nterm:
        if len(cen) > 1:
            gaps = [(cen[i + 1] - cen[i], i) for i in range(len(cen) - 1)]
            gmax, i = max(gaps)
        else:
            gmax = 0.0
        if gmax > 0.1:
            cen.insert(i + 1, 0.5 * (cen[i] + cen[i + 1]))
        else:
            cen = [cen[0] - 0.3] + cen + [cen[-1] + 0.3]
            cen = sorted(cen)[:nterm] if len(cen) > nterm else cen
    return sorted(10.0 ** np.asarray(cen[:nterm]))


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
    path = os.path.join(OUTDIR, f"cooling_{elem}_{stage}.txt")
    with open(path, "w") as fh:
        fh.write(f"# MoCHII cooling coefficients: {SPEC_LABEL[name]}"
                 f" ({elem}_{stage}) collisional line cooling\n")
        fh.write("# form: Lambda(T) = T^-1/2 * sum_i A_i * exp(-T_i/T)"
                 "   [erg cm^3 s^-1] per (n_e * n_ion)\n")
        fh.write(f"# source: CHIANTI v{CHIANTI_VERSION} .elvlc/.scups,"
                 " Burgess-Tully descaling (tools/fitting/chianti_cooling.py)\n")
        fh.write("# population: true low-density statistical-equilibrium"
                 " populations (n_e -> 0 limit of the n-level solve: all ions"
                 " in the lowest level,\n")
        fh.write("#   excitation out of it only, each excitation radiating the"
                 " full dE_1u through the cascade); optically thin\n")
        fh.write("# energies: observed .elvlc level differences, in both the"
                 " Boltzmann factor and the emitted photon\n")
        fh.write("#   (the .scups dE is the Burgess-Tully scaling energy, not"
                 " a physical transition energy)\n")
        fh.write("# caveat: no density suppression above n_crit; the FIR"
                 " fine-structure lines ([C II] 158um, [O III] 52/88um) and\n")
        fh.write("#   the [O II]/[S II] optical doublets saturate above their"
                 " n_crit, where this limit overestimates the cooling;\n")
        fh.write("#   density-dependent emissivities come from the n-level"
                 " solve\n")
        fh.write(f"# fit range: {TMIN:.1e} - {TMAX:.1e} K;"
                 f" max fit error: {err*100:.2f}%"
                 f" ({rel_err*100:.2f}% where Lambda > 1e-3 max)\n")
        fh.write(f"# generated by tools/fitting/fit_cooling.py on {DATE}\n")
        fh.write(f"{len(A)}\n")
        for a, t in zip(A, Ti):
            fh.write(f"  {a:.8e}  {t:.6e}\n")
    return path


results = {}
print(f"{'ion':6s} {'nterm':>5s} {'maxerr_fit':>11s} {'maxerr_relevant':>16s}")
for name, (el, st) in IONS.items():
    lam = cooling_low_density_limit(el, st, T)
    best = None
    for nterm in range(2, NTERM_MAX + 1):
        base = np.asarray(ground_excitation_ladder(el, st, nterm))
        rng = np.random.default_rng(12345 + nterm)
        starts = [base] + [base * 10.0 ** rng.uniform(-0.25, 0.25, nterm)
                           for _ in range(NSTART - 1)]
        for guess in starts:
            try:
                A, Ti = fit_multiexp(lam, nterm, np.sort(guess))
            except Exception:
                continue
            fit = evaluate(A, Ti, T)
            err = np.abs(fit / lam - 1.0).max()
            if best is None or err < best[0]:
                best = (err, A, Ti, fit)
        if best[0] < ERR_TARGET:
            break
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
fig.savefig("cooling_fits.png", dpi=150)
print("wrote cooling_fits.png")
