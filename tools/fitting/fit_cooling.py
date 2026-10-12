"""Cooling fits for the MoCHII HII-region ion set.

Fits identified CHIANTI or Stout optically thin collisional line-cooling curves
Lambda(T) per (n_e * n_ion) with the multi-exponential form used in EXHALE
(fit_cno_formulas.py, same form as its Fe II fit):

    Lambda(T) = T^{-1/2} * sum_i A_i * exp(-T_i / T)   [erg cm^3 s^-1]

for the diagnostic ion set: O II, O III, N II, S II, S III, Ne II, Ne III,
and the rest of the registry.

Population model: the stationary n-level balance at n_e = 1e-6 cm^-3, the
reference column of the runtime cooling tables (target_log_cooling,
chianti_cooling.reference_log_cooling), not the strict n_e -> 0 limit: a
placeholder spontaneous rate (CHIANTI S V and Ar VII 3s3p 3P0 -> 3s2 1S0,
A = 1e-30 s^-1) drains a level in the strict limit that is collisionally
de-excited at any physical density. Transition
energies are the selected level differences, the same energies the n-level
solve uses in both the Boltzmann factor and the emitted photon. For CHIANTI,
the .scups dE is a Burgess-Tully scaling energy, not a physical transition
energy. Ar II and Ca IV use the published two-level Stout temperature tables;
their Upsilon is held at its first node below the published temperature range
over the 1,000--100,000 K numerical interval. Ar II is not covered by the
paper's better-than-10-percent accuracy estimate. Neither table is continued
above its last source node.

CAVEAT (carried into each output header): the fits are low-density-limit
cooling.  Transitions with low critical density (e.g. [O II] 3726/3729,
[S II] 6717/6731, the fine-structure IR lines, and especially [C II] 158um
whose n_crit ~ 50 cm^-3) saturate above n_crit; the cooling fits
overestimate cooling there.  Density-dependent line ratios and emissivities
come from the n-level solve (Upsilon fits + statistical equilibrium).

Writes:  <output-dir>/cooling_<el>_<stage>.txt, default ../../data/atomic
         (provenance header: source, population model, fit range,
         fit form, max fit error)
         cooling_fits.png (fit-vs-source comparison) unless --no-plot
Run:     python3 fit_cooling.py [--ion NAME | --ion-spec EL:STAGE ...]
                                [--output-dir DIR] [--no-plot]
"""

import os
import argparse
import io
import tempfile
from pathlib import Path
from hashlib import sha256
import numpy as np
from scipy.optimize import least_squares
from chianti_cooling import (cooling_low_density_limit, log_cooling_ground_state_limit, ion_dir, read_elvlc,
                             excitation_levels, HC_ERG_CM, REFERENCE_DENSITY, reference_log_cooling,
                             read_scups, read_wgfa, read_excitation_atom, bound_excitation_atom, upsilon, selected_energies,
                             transition_energy_erg, DBASE, K_B_ERG, COLL_PREF,
                             STOUT_EXCITATION_IONS)
from atomic_output import atomic_write_text

from datetime import datetime
from zoneinfo import ZoneInfo

DATE = datetime.now(ZoneInfo("Asia/Seoul")).date().isoformat()
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "..", "data", "atomic")

with open(os.path.join(DBASE, "VERSION")) as fh:
    CHIANTI_VERSION = fh.read().strip()

T = np.logspace(3.0, 5.0, 201)
TMIN, TMAX = T[0], T[-1]
# Acceptance grid: 801 temperatures over the same interval (every fit node
# T is one of them), the grid tests/ion_stages/check_excitation_atoms.py
# compares each installed file on. A fit accepted on its 201 nodes alone can
# exceed its stated error between them: O IV (o_4), whose source Upsilon is
# held constant below the 3200 K start of the BT interval, met 1.00% on the
# 201 nodes and reached 1.13% at 3199 K.
CHECK_T = np.logspace(3.0, 5.0, 801)

NTERM_MAX = 8         # positive-coefficient exponential terms
SIGNED_NTERM_MAX = 14 # mathematical expansion of total positive cooling
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
    "SI":    ("s", 1),
    "SII":   ("s", 2),
    "SIII":  ("s", 3),
    # S IV: the [S IV] 10.51um ground fine-structure line (951 cm^-1 = 1369 K)
    # plus the 3s3p^2 UV blocks; a ~1.2% coolant of the HII40 benchmark.
    "SIV":   ("s", 4),
    "NeII":  ("ne", 2),
    "NeIII": ("ne", 3),
    "SiI":   ("si", 1),
    # Ar II uses the published ground fine-structure pair from Stout.
    "ArII":  ("ar", 2),
    "ArIII": ("ar", 3),
    "ArIV":  ("ar", 4),
    # G5: magnesium and iron (gas-phase; strongly depleted onto grains).
    "MgII":  ("mg", 2),
    "FeII":  ("fe", 2),
    "FeIII": ("fe", 3),
    # This default selection does not include every installed ion. C IV,
    # Si V, and Fe IV can be selected with --ion-spec. Cl V and Ca III still
    # need accepted excitation sources; an assumed nebular temperature does
    # not establish that their cooling can be omitted.
    "SiII":  ("si", 2),
    "SiIII": ("si", 3),
    "SiIV":  ("si", 4),
    "ClII":  ("cl", 2),
    "ClIII": ("cl", 3),
    "ClIV":  ("cl", 4),
    "CaII":  ("ca", 2),
    "CaIV":  ("ca", 4),
    "CaV":   ("ca", 5),
}

SPEC_LABEL = {  # for plot titles (ASCII / LaTeX-safe)
    "HI": "H I", "CI": "[C I]", "CII": "[C II]", "CIII": "C III",
    "NI": "[N I]", "NIII": "[N III]", "OI": "[O I]",
    "OII": "[O II]", "OIII": "[O III]", "NII": "[N II]",
    "SII": "[S II]", "SIII": "[S III]", "SIV": "[S IV]", "NeII": "[Ne II]",
    "NeIII": "[Ne III]", "ArIII": "[Ar III]", "ArIV": "[Ar IV]",
    "MgII": "Mg II", "FeII": "[Fe II]", "FeIII": "[Fe III]",
    "SiI": "[Si I]",
    "SiII": "[Si II]", "SiIII": "Si III", "SiIV": "Si IV",
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
    lev, radiation, collisions = bound_excitation_atom(elem, stage)
    Tu, contrib = [], []
    for tr in collisions:
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
    # to be carried by superposing neighboring terms.  Pad by splitting the
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


def ground_excitation_threshold(elem, stage):
    """Lowest collisional excitation temperature out of the ground [K]."""
    levels, _, collisions = bound_excitation_atom(elem, stage)
    ground = min(levels, key=lambda i: levels[i][1])
    gaps = [transition_energy_erg(levels, tr) for tr in collisions
            if ground in (tr['ll'], tr['ul'])]
    return min(gap for gap in gaps if gap > 0) / K_B_ERG


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
    sol = least_squares(resid, th0, method="lm", max_nfev=6000)
    A = np.exp(np.clip(sol.x[:nterm], -700.0, 700.0))
    Ti = np.exp(np.clip(sol.x[nterm:], -700.0, 700.0))
    o = np.argsort(Ti)
    return A[o], Ti[o]


def fit_ground_state_log_cooling(log_source, nterm, Tguess):
    """Fit positive exponential terms without underflowing excitation tails."""
    from scipy.special import logsumexp
    target = log_source+.5*np.log(T)
    start_amplitude = target[-1]+min(Tguess)/T[-1]-np.log(nterm)
    def residual(parameters):
        temperatures=np.exp(np.clip(parameters[nterm:],-700,700))
        return logsumexp(parameters[:nterm,None]-temperatures[:,None]/T,axis=0)-target
    initial=np.r_[np.full(nterm,start_amplitude),np.log(Tguess)]
    result=least_squares(residual,initial,method='lm',max_nfev=6000)
    amplitude=np.exp(np.clip(result.x[:nterm],-700,700))
    temperatures=np.exp(np.clip(result.x[nterm:],-700,700))
    order=np.argsort(temperatures)
    return amplitude[order],temperatures[order]


def log_cooling_fit(amplitude, excitation, temperature):
    """Natural log of a positive total fit, preserving excitation underflow."""
    base=float(np.min(excitation))
    scaled=np.sum(amplitude[:,None]*np.exp(-(excitation[:,None]-base)/temperature),axis=0)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.log(scaled)-base/temperature-.5*np.log(temperature)


def log_relative_cooling_error(amplitude, excitation, log_source, temperature=None):
    """Relative total-cooling error with a removed physical excitation factor.

    Signed terms remain permitted, but their scaled sum must be positive.
    The logarithmic source is used only for the one-end-state model.
    """
    log_fit=log_cooling_fit(amplitude,excitation,T if temperature is None else temperature)
    if not np.all(np.isfinite(log_fit)):return float('inf')
    with np.errstate(over='ignore'):
        return float(np.max(np.abs(np.expm1(log_fit-log_source))))


def fit_signed_exponentials(lam, nterm, excitation_temperature, log_source=None):
    """Fit total positive cooling with signed exponential coefficients.

    Coefficients are a mathematical expansion, not individual emission rates.
    Stable-state populations can give temperature dependence that a positive
    exponential sum cannot approximate adequately. Solve linear coefficients
    at each trial temperature ladder, with a common excitation exponential
    removed to retain precision in the cold tail.
    """
    base = float(excitation_temperature)
    if log_source is None:
        target = lam*np.sqrt(T)*np.exp(base/T)
    else:
        target = np.exp(log_source+.5*np.log(T)+base/T)
    scale = target.max()
    target /= scale

    if log_source is not None:
        # Resolve the temperature variation of Upsilon after removing the
        # physical excitation threshold. Multiplicative ladders of the full
        # threshold leave gaps of millions of kelvin for these ions and miss
        # their low-temperature variation. Offsets are numerical basis terms,
        # not additional physical excitation energies.
        best = None
        extended_temperature = np.logspace(0, 9, 10001)
        for lowest in (100., 300., 1000., 3000., 10000.):
            for highest in (1e5, 3e5, 1e6):
                offsets = np.r_[0., np.geomspace(lowest, highest, nterm-1)]
                matrix = np.exp(-offsets[None,:]/T[:,None])/target[:,None]
                coefficient = np.linalg.lstsq(matrix, np.ones(T.size), rcond=1e-10)[0]
                fitted = matrix @ coefficient
                condition = np.max(np.sum(np.abs(matrix*coefficient),axis=1)/np.abs(fitted))
                extended = np.exp(-offsets[None,:]/extended_temperature[:,None]) @ coefficient
                error = float(np.max(np.abs(fitted-1)))
                if condition < 1e5 and np.all(extended >= 0) and (best is None or error < best[0]):
                    best = (error, coefficient*scale, base+offsets)
        if best is not None:
            return best[1], best[2]

    def solution(log_temperature):
        temperatures = np.exp(log_temperature)
        matrix = np.exp(-(temperatures[None,:]-base)/T[:,None])/target[:,None]
        # Remove nearly dependent columns before they create large cancelling terms.
        coefficients = np.linalg.lstsq(matrix, np.ones(T.size), rcond=1e-8)[0]
        return matrix@coefficients-1, coefficients*scale, temperatures

    initial = np.log(np.geomspace(base, max(1e6, 4*base), nterm))
    result = least_squares(lambda x: solution(x)[0], initial,
                           bounds=(np.log(.5*base), np.log(max(1e7, 40*base))),
                           max_nfev=1000)
    residual, coefficients, temperatures = solution(result.x)
    order = np.argsort(temperatures)
    return coefficients[order], temperatures[order]


def constrained_signed_exponentials(elem, stage, excitation_temperature):
    """Minimax total-cooling fit after shorter exponential sums fail.

    Uniform temperature offsets form a mathematical basis, not additional
    excitation thresholds. The corrected physical source is evaluated on 801
    temperatures. A linear program bounds relative error and cancellation in
    the fit interval. Positivity is checked on sampled temperatures from
    1 to 1e9 K; this does not establish extrapolation accuracy.
    """
    from scipy.optimize import linprog
    temperature = np.logspace(np.log10(TMIN), np.log10(TMAX), 801)
    log_source = target_log_cooling(elem, stage, temperature)
    if not np.all(np.isfinite(log_source)):
        return None
    base = float(excitation_temperature)
    log_target = log_source + 0.5*np.log(temperature) + base/temperature
    maximum = float(np.max(log_target))
    normalized = np.exp(log_target - maximum)
    if np.any(normalized <= 0) or not np.all(np.isfinite(normalized)):
        return None
    with np.errstate(over='ignore'):
        amplitude_scale = np.exp(maximum)
    if not np.isfinite(amplitude_scale):
        return None
    extended_temperature = np.logspace(0, 9, 401)
    full_extended_temperature = np.logspace(0, 9, 10001)
    dense_temperature = np.logspace(np.log10(TMIN), np.log10(TMAX), 1601)
    dense_log_source = None
    best = None
    cancellation_limit = 1e5
    for count in (18, 24, 32, 40, 48, 56, 64):
        for step_fraction in (.01, .02, .03, .04, .05, .06, .08):
            offsets = np.arange(count, dtype=float)*(step_fraction*TMAX)
            basis = np.exp(-offsets[None,:]/temperature[:,None])
            relative = basis/normalized[:,None]
            extended_basis = np.exp(-offsets[None,:]/extended_temperature[:,None])
            ones = np.ones(len(temperature))
            zeros = np.zeros(len(temperature))
            constraints = np.vstack([
                np.column_stack([relative, -relative, -ones]),
                np.column_stack([-relative, relative, -ones]),
                np.column_stack([basis*(1/cancellation_limit-1),
                                 basis*(1/cancellation_limit+1), zeros]),
                np.column_stack([extended_basis*(1/cancellation_limit-1),
                                 extended_basis*(1/cancellation_limit+1),
                                 np.zeros(len(extended_temperature))])])
            rhs = np.r_[ones, -ones, zeros, np.zeros(len(extended_temperature))]
            solution = linprog(np.r_[np.zeros(2*count), 1.0],
                               A_ub=constraints, b_ub=rhs,
                               bounds=[(0,None)]*(2*count+1), method='highs')
            if not solution.success:
                continue
            coefficient = solution.x[:count]-solution.x[count:2*count]
            amplitude = np.array([float(f'{value*amplitude_scale:.16e}')
                                  for value in coefficient])
            excitation = np.array([float(f'{base+offset:.16e}') for offset in offsets])
            rounded_offsets = excitation - np.min(excitation)
            terms = amplitude[None,:]*np.exp(-rounded_offsets[None,:]/temperature[:,None])
            total = np.sum(terms, axis=1)
            if np.any(total <= 0):
                continue
            condition = np.max(np.sum(np.abs(terms),axis=1)/total)
            if condition >= cancellation_limit:
                continue
            extended = evaluate(amplitude, excitation, full_extended_temperature)
            if not np.all(np.isfinite(extended)) or np.any(extended < 0):
                continue
            # The rounded file coefficients must meet the actual source on
            # both the fit and independent intermediate temperatures.
            fitted = evaluate(amplitude, excitation, temperature)
            if np.any(fitted <= 0):
                continue
            error = float(np.max(np.abs(np.expm1(np.log(fitted)-log_source))))
            if error >= ERR_TARGET:
                continue
            if dense_log_source is None:
                dense_log_source = target_log_cooling(elem, stage, dense_temperature)
            dense_fitted = evaluate(amplitude, excitation, dense_temperature)
            if np.any(dense_fitted <= 0):
                continue
            dense_error = float(np.max(np.abs(np.expm1(np.log(dense_fitted)-dense_log_source))))
            error = max(error, dense_error)
            if error < ERR_TARGET and (best is None or error < best[0]):
                active = amplitude != 0.0
                best = (error, amplitude[active], excitation[active])
    return best


def evaluate(A, Ti, T):
    return np.sum(A[:, None] * np.exp(-Ti[:, None] / T[None, :]), axis=0) / np.sqrt(T)


def target_log_cooling(elem, stage, temperature):
    """Natural log of the cooling coefficient the fit is made to: the
    stationary balance at REFERENCE_DENSITY = 1e-6 cm^-3.

    Decision of 2026-10-11 (main session): the target is the n_e = 1e-6
    cm^-3 stationary state, the runtime reference column, not the strict
    n_e -> 0 limit. Placeholder spontaneous rates (CHIANTI 11.0.2 A = 1e-30
    s^-1 for the one-photon-forbidden 3s3p 3P0 -> 3s2 1S0 decays of S V and
    Ar VII) make the strict limit radiate every excitation into a level that
    at any physical density (n_e above ~A/q ~ 1e-22 cm^-3) is collisionally
    de-excited instead, and the runtime divides by the 1e-6 column
    (nlevel_cooling_mod::nlevel_cooling_apply), so a fit to the strict limit
    is off by the ratio of the two: S V 0.16-0.38% above 1.8e4 K and 2.7% at
    1e3 K, Ar VII 11.6% at 1e3 K. Where the finite-density cooling underflows
    (below 1e-280 erg cm^3 s^-1) the logarithmic strict ground-state sum
    takes over (chianti_cooling.reference_log_cooling).
    """
    levels, radiation, collisions = bound_excitation_atom(elem, stage)
    return reference_log_cooling(levels, radiation, collisions, temperature)[0]


class AcceptanceReference:
    """The fit target on CHECK_T and the fit errors against it."""

    def __init__(self, elem, stage):
        self.log_target = target_log_cooling(elem, stage, CHECK_T)
        self.significant = self.log_target > np.max(self.log_target)+np.log(1e-3)

    def errors(self, A, Ti):
        """(maximum relative error on CHECK_T, the same where the target
        exceeds 1e-3 of its maximum)."""
        log_fit = log_cooling_fit(A, Ti, CHECK_T)
        if not np.all(np.isfinite(log_fit)):
            return float('inf'), float('inf')
        with np.errstate(over='ignore'):
            error = np.abs(np.expm1(log_fit-self.log_target))
        return float(error.max()), float(error[self.significant].max())


def write_ion_file(name, elem, stage, A, Ti, err, rel_err, output_dir=OUTDIR):
    path = os.path.join(output_dir, f"cooling_{elem}_{stage}.txt")
    fh = io.StringIO()
    fh.write(f"# MoCHII cooling coefficients: {SPEC_LABEL.get(name, name)}"
             f" ({elem}_{stage}) collisional line cooling\n")
    fh.write("# form: Lambda(T) = T^-1/2 * sum_i A_i * exp(-T_i/T)"
             "   [erg cm^3 s^-1] per (n_e * n_ion)\n")
    fh.write("# coefficients represent the total fit; individual terms may be signed\n")
    if (elem.lower(),stage) in [('ar',2),('ca',4)]:
        from stout_fine_structure_atoms import SOURCES, STOUT_ROOT, stout_fine_structure_atom
        key = (elem.lower(),stage)
        source = SOURCES[key]
        _, _, fine_structure_collisions = stout_fine_structure_atom(elem, stage)
        first = fine_structure_collisions[0]['temperature_grid'][0]
        fh.write('# source: Pelan and Berrington (1995), A&AS 110, 209, Table 3; Cloudy c25.00 Stout ground fine-structure pair\n')
        fh.write('# model scope: two ground-term levels; higher radiative levels lack identified electron collision strengths\n')
        fh.write(f'# source validity: log10(T/K) {source["source_log_t"][0]:.1f}--{source["source_log_t"][-1]:.1f}; no Upsilon continuation above the final source node\n')
        fh.write(f'# adopted cold approximation: constant first-node Upsilon below {first:.0f} K; actual T enters Maxwellian rates\n')
        fh.write('# interpolation: log-PCHIP in T and Upsilon between original 11 Stout nodes\n')
        if key == ('ar',2):
            fh.write('# source accuracy: the paper excludes Ar II from its better-than-10-percent estimate\n')
        else:
            fh.write('# source accuracy: the paper estimates better than 10 percent for Ca IV within Table 3 range\n')
        stem = STOUT_ROOT / elem.lower() / f'{elem.lower()}_{stage}' / f'{elem.lower()}_{stage}'
        for extension in ('nrg', 'tp', 'coll'):
            source_path = Path(f'{stem}.{extension}')
            fh.write(f'# source SHA-256 {source_path}: {sha256(source_path.read_bytes()).hexdigest()}\n')
    elif (elem.lower(),stage) == ("cl",7):
        fh.write(f"# source: Cloudy c25.00 Stout {elem}_{stage}.nrg/.tp/.coll; bound states only\n")
        fh.write("# log-PCHIP Upsilon at 9800--9.8e7 K; constant Upsilon below 9800 K\n")
        fh.write("# adopted cold resonance approximation has unverified physical accuracy\n")
    else:
        if (elem.lower(), stage) == ('si', 1):
            from silicon_source_records import NIST_SOURCE
            from silicon_spontaneous_probabilities import INFERRED_BT_FLOOR_K
            source_code = Path(__file__).with_name('silicon_spontaneous_probabilities.py')
            parser_code = Path(__file__).with_name('silicon_source_records.py')
            fh.write(f'# source: CHIANTI v{CHIANTI_VERSION} Si I levels and BT collision strengths; '
                     'original NIST 2002 MCHF spontaneous strengths\n')
            fh.write('# radiation: 185 matched CHIANTI E1 pairs use NIST length-form strengths at selected photon frequencies; '
                     '21 missing NIST pairs are added\n')
            fh.write('# shared E1/M2 contributions are summed; 44 existing non-E1 CHIANTI rates remain approximations\n')
            si_stem = Path(ion_dir('si', 1)) / 'si_1'
            for extension in ('elvlc', 'wgfa', 'scups'):
                source_path = Path(f'{si_stem}.{extension}')
                fh.write(f'# source SHA-256 {source_path}: {sha256(source_path.read_bytes()).hexdigest()}\n')
            fh.write(f'# NIST source SHA-256 {NIST_SOURCE}: {sha256(NIST_SOURCE.read_bytes()).hexdigest()}\n')
            for source_path in (source_code, parser_code):
                fh.write(f'# selection code SHA-256 {source_path}: {sha256(source_path.read_bytes()).hexdigest()}\n')
        else:
            fh.write(f"# source: CHIANTI v{CHIANTI_VERSION} .elvlc/.scups/.wgfa,"
                     " Burgess-Tully descaling (tools/fitting/chianti_cooling.py)\n")
        _, _, source_collisions = read_excitation_atom(elem, stage)
        bt_minima = {row.get('bt_tmin_K') for row in source_collisions}
        if bt_minima and bt_minima != {None}:
            from write_excitation_atom import BT_SOURCE_TMAX_K
            key = (elem.lower(), stage)
            if len(bt_minima) != 1 or None in bt_minima:
                raise ValueError(f'Incomplete BT source temperature interval for {elem}_{stage}')
            bt_minimum = bt_minima.pop()
            if key == ('ne', 1):
                if bt_minimum != 1000.36:
                    raise ValueError('Ne I inferred first-node floor changed')
                fh.write('# radiation: 142 original CHIANTI A rows plus NIST-dated Stout Ne I A(3,15)=900000 and A(5,15)=68200000 s^-1\n')
                fh.write('# CHIANTI rounded BT first finite nodes infer 999.410754--1000.355160 K; not a published physical interval\n')
                fh.write('# adopted numerical Upsilon floor: 1000.36 K; actual gas T enters Maxwellian rates and detailed balance\n')
                fh.write('# numerical atomic-data approximation: 1000--100000 K; source conversion and atomic accuracy unverified\n')
            elif key == ('si', 1):
                if bt_minimum != INFERRED_BT_FLOOR_K:
                    raise ValueError('Si I inferred first-node floor changed')
                fh.write('# CHIANTI rounded BT first finite nodes infer 997.095--1003.020 K; not a published physical interval\n')
                fh.write('# adopted numerical Upsilon floor: 1003.02 K; actual gas T enters Maxwellian rates and detailed balance\n')
                fh.write('# numerical atomic-data approximation: 1000--100000 K; original NIST strength accuracy and final 2005 article methods unverified\n')
            else:
                if key not in BT_SOURCE_TMAX_K:
                    raise ValueError(f'Unidentified BT source upper temperature for {elem}_{stage}')
                bt_maximum = BT_SOURCE_TMAX_K[key]
                if (not np.isfinite(bt_minimum) or bt_minimum <= 0 or
                        not np.isfinite(bt_maximum) or bt_maximum <= bt_minimum):
                    raise ValueError(f'Invalid BT source temperature interval for {elem}_{stage}')
                fh.write(f'# published BT collision temperature interval: '
                         f'{bt_minimum:.0f}--{bt_maximum:.0f} K\n')
                fh.write(f'# adopted cold approximation: constant source Upsilon below '
                         f'{bt_minimum:.0f} K; actual T enters Maxwellian rates and detailed balance\n')
                fh.write('# cold collision-strength accuracy is unverified below the source interval; '
                         'no continuation above its upper bound\n')
    fh.write("# population: stationary electron-collision/spontaneous-radiation"
             f" balance at n_e = {REFERENCE_DENSITY:.0e} cm^-3 (the reference density of the cooling"
             " tables), not the strict n_e -> 0 limit; optically thin\n")
    fh.write("# scope: collision-connected states below the physical ionization edge; no autoionization\n")
    if (elem.lower(),stage) not in STOUT_EXCITATION_IONS:
        fh.write("# interpolation: source cubic unless nonnegative nodes have negative extrema, then PCHIP\n")
        fh.write("# interpolation preserves source nodes; unresolved resonance accuracy is not certified\n")
    fh.write("#   collisions and radiative cascades determine populations of"
             " radiatively stable states; no photoexcitation\n")
    if (elem.lower(),stage) in [('ar',2),('ca',4)]:
        fh.write('# energies: Stout observed ground fine-structure pair, checked below the ionization edge\n')
    elif (elem.lower(),stage) == ('cl',7):
        fh.write('# energies: Stout bound states; above-threshold core states excluded\n')
    elif elem.lower() == 'n' and stage == 4:
        from niv_level_energies import STOUT_SHA256
        fh.write(f"# energies: Cloudy c25.00 Stout n_4.nrg SHA-256 {STOUT_SHA256}\n")
        fh.write("# observed levels 1--29; BSR above; original CHIANTI A and collision indices\n")
        fh.write("# low-temperature resonance uncertainty and convergence above level 87\n")
        fh.write("# remain unquantified; fit interval is numerical, not a physical accuracy claim\n")
    elif (elem.lower(), stage) in [('o',5),('cl',6)]:
        fh.write('# energies: checked observed/ICFT vector in atomic_excitation_energies.py\n')
        fh.write('# O V: observed 1--40, identity-mapped Stout ICFT above; Cl VI: observed 1--44, ICFT above\n')
        fh.write('# unshifted scattering resonances and higher-level convergence remain uncertain\n')
        fh.write('# fit interval is numerical, not a physical accuracy claim\n')
    elif (elem.lower(),stage) in [('o',4),('ne',5),('s',7),('ar',5),('mg',4)]:
        fh.write('# energies: observed complete valence complexes and original theoretical higher states\n')
        fh.write('# checked identities and term centroids in valence_excitation_energies.py\n')
        fh.write('# unshifted scattering resonances have unverified low-temperature accuracy\n')
    elif (elem.lower(),stage)==('s',6):
        fh.write('# energies: unresolved observed centers with statistical-weight-centered theoretical splittings\n')
        fh.write('# no artificial fine-structure gap floor; unshifted scattering resonances remain uncertain\n')
    elif (elem.lower(),stage)==('mg',3):
        fh.write('# energies: selected observed/theoretical values except source ID 167, whose misassigned J=1 observation is rejected for its J=0 state\n')
        fh.write('# source ID 167 uses its original theoretical 634430 cm^-1 energy; all photon and Maxwellian gaps use this selection\n')
        fh.write('# higher E1 probabilities preserve dipole strength with physical frequency cubed and emitting statistical weight; first 27 A values are unchanged\n')
        fh.write('# higher gf=0 A values are unidentified-multipole effective-rate approximations; atomic accuracy remains unverified\n')
    else:
        fh.write("# energies: observed .elvlc values when available, otherwise theoretical;"
                 " the same differences enter excitation and cascade energy\n")
    if (elem.lower(),stage) not in STOUT_EXCITATION_IONS:
        fh.write("#   (the .scups dE is the Burgess-Tully scaling energy, not"
                 " a physical transition energy)\n")
    fh.write("# caveat: no density suppression above n_crit; the FIR"
             " fine-structure lines ([C II] 158um, [O III] 52/88um) and\n")
    fh.write("#   the [O II]/[S II] optical doublets saturate above their"
             " n_crit, where this limit overestimates the cooling;\n")
    fh.write("#   density-dependent emissivities come from the n-level"
             " solve\n")
    fh.write(f"# fit range: {TMIN:.1e} - {TMAX:.1e} K, fitted on {len(T)} source temperatures;"
             f" checked on {len(CHECK_T)} temperatures (the fit nodes included);"
             f" max fit error on these temperatures: {err*100:.2f}%"
             f" ({rel_err*100:.2f}% where Lambda > 1e-3 max)\n")
    fh.write(f"# generated by tools/fitting/fit_cooling.py on {DATE}\n")
    fh.write(f"{len(A)}\n")
    for a, t in zip(A, Ti):
        fh.write(f"  {a:.16e}  {t:.16e}\n")
    atomic_write_text(path, fh.getvalue())
    return path


def generate(ions, output_dir=OUTDIR, plot_path=None):
    results = {}
    print(f"{'ion':6s} {'nterm':>5s} {'maxerr_fit':>11s} {'maxerr_relevant':>16s}")
    for name, (el, st) in ions.items():
        log_target = target_log_cooling(el, st, T)
        if not np.all(np.isfinite(log_target)):
            raise ValueError(f"{name}: no finite target cooling")
        lam = np.exp(log_target)
        log_source = None if np.all(lam > 0) else log_target
        references = AcceptanceReference(el, st)
        best = None
        for nterm in range(2, NTERM_MAX + 1):
            base = np.asarray(ground_excitation_ladder(el, st, nterm))
            rng = np.random.default_rng(12345 + nterm)
            starts = [base] + [base * 10.0 ** rng.uniform(-0.25, 0.25, nterm)
                               for _ in range(NSTART - 1)]
            for guess in starts:
                try:
                    if log_source is None:
                        A, Ti = fit_multiexp(lam, nterm, np.sort(guess))
                    else:
                        A, Ti = fit_ground_state_log_cooling(log_source,nterm,np.sort(guess))
                except Exception:
                    continue
                fit = evaluate(A, Ti, T)
                err = references.errors(A, Ti)[0]
                if best is None or err < best[0]:
                    best = (err, A, Ti, fit)
            if best is not None and best[0] < ERR_TARGET:
                break
        if best is None:
            raise ValueError(f"{name}: no cooling fit converged")
        if best[0] >= ERR_TARGET:
            base = min(ground_excitation_ladder(el, st, NTERM_MAX))
            if log_source is not None:
                # The removed factor exp(-base/T) must be the physical
                # excitation threshold: a padded ladder value far below it
                # (O VII: 3.26e6 K against 6.51e6 K) leaves exp(log_source
                # + base/T) underflowing to zero across the cold part of
                # the range, and the relative least-squares matrix is then
                # infinite. The threshold is used only where the source
                # itself underflows, i.e. in the logarithmic branch.
                base = max(base, ground_excitation_threshold(el, st))
            for nterm in range(6, SIGNED_NTERM_MAX+1):
                A, Ti = fit_signed_exponentials(lam, nterm, base,log_source)
                fit = evaluate(A, Ti, T)
                err = references.errors(A, Ti)[0]
                # Outside the source fit range only positivity is checked;
                # this is not a claim of extrapolation accuracy.
                extended = evaluate(A, Ti, np.logspace(0, 9, 10001))
                terms = A[:,None]*np.exp(-(Ti[:,None]-min(Ti))/T[None,:])
                condition = np.max(np.sum(np.abs(terms),axis=0)/np.abs(np.sum(terms,axis=0)))
                if np.all(np.isfinite(extended)) and np.all(extended >= 0) and condition < 1e5 and err < best[0]:
                    best = (err, A, Ti, fit)
                if best[0] < ERR_TARGET:
                    break
        if best[0] >= ERR_TARGET:
            base = min(ground_excitation_ladder(el, st, NTERM_MAX))
            if log_source is not None:
                base = max(base, ground_excitation_threshold(el, st))
            constrained = constrained_signed_exponentials(el, st, base)
            if constrained is not None:
                _, A, Ti = constrained
                best = (references.errors(A, Ti)[0], A, Ti, evaluate(A, Ti, T))
        if best[0] >= ERR_TARGET:
            raise ValueError(f"{name}: fit error {best[0]:.6g} exceeds {ERR_TARGET}")
        _, A, Ti, fit = best
        err, rel_err = references.errors(A, Ti)
        path = write_ion_file(name, el, st, A, Ti, err, rel_err, output_dir)
        print(f"{name:6s} {len(A):5d} {err*100:10.2f}% {rel_err*100:15.2f}%"
              f"   -> {path}")
        results[name] = dict(elem=el, stage=st, lam=lam, fit=fit, A=A, Ti=Ti)
    if plot_path is not None and results:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        nrow = (len(results) + 3)//4
        fig, axes = plt.subplots(nrow, 4, figsize=(15, 3.5*nrow), sharex=True, squeeze=False)
        for ax, (name, r) in zip(axes.flat, results.items()):
            ax.loglog(T, r["lam"], "k-", lw=1.8, label="Source")
            ax.loglog(T, r["fit"], "r--", lw=1.4, label="fit")
            ax.set_title(SPEC_LABEL.get(name, name))
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
        parent = os.path.dirname(os.path.abspath(plot_path))
        os.makedirs(parent, exist_ok=True)
        handle, temporary = tempfile.mkstemp(prefix=".cooling_fits.", suffix=".png", dir=parent)
        os.close(handle)
        try:
            fig.savefig(temporary, dpi=150)
            os.replace(temporary, plot_path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
            plt.close(fig)
        print(f"wrote {plot_path}")
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ion", action="append", choices=list(IONS),
                        help="Generate a configured ion; repeat for several")
    parser.add_argument("--ion-spec", action="append", metavar="ELEMENT:STAGE",
                        help="Generate an additional ion with an identified excitation source")
    parser.add_argument("--output-dir", default=None, help="Generated file directory")
    parser.add_argument("--no-plot", action="store_true", help="Skip comparison plot")
    args = parser.parse_args(argv)
    ions = {name: IONS[name] for name in (args.ion or ([] if args.ion_spec else IONS))}
    for spec in args.ion_spec or ():
        el, separator, stage = spec.partition(":")
        if not separator or not stage.isdigit() or not el.isalpha():
            parser.error(f"invalid ion specification: {spec}")
        ions[f"{el.upper()} {int(stage)}"] = (el.lower(), int(stage))
    output_dir = args.output_dir or OUTDIR
    plot_path = None if args.no_plot else os.path.join(args.output_dir or ".", "cooling_fits.png")
    generate(ions, output_dir, plot_path)


if __name__ == "__main__":
    main()
