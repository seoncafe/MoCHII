"""n-level atom data for the MoCHII generic n-level solver.

For each registry ion, extracts the lowest NLEV levels and produces one
data/atomic/nlevel_<ion>.txt file containing

  1. level energies and statistical weights (.elvlc, observed energies
     preferred, else theoretical),
  2. radiative A-values among those levels (.wgfa; duplicate (l,u) rows
     summed),
  3. effective collision strengths Upsilon(T) as low-order Chebyshev fits in
     Burgess-Tully (1992) scaled space, one for each transition (.scups).

Transitions are written ORIENTED BY ENERGY: the first index is the
lower-energy level, the second the upper (emitting) level.  CHIANTI does not
order the two index columns consistently -- o_2.wgfa lists 2D3/2 - 2D5/2
(upper level first, because 2D3/2 is the higher term), n_1.wgfa lists
2D5/2 - 2D3/2 (lower level first), and fe_3.scups has nine index-ordered but
energy-inverted rows -- so the level energies, not the column order, decide
which level emits and which one is de-excited collisionally.

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
Run:     python3 fit_nlevel.py
"""

import os
import argparse
import io
import numpy as np
from numpy.polynomial import chebyshev as C
from chianti_cooling import (ion_dir, read_elvlc, read_wgfa, read_scups, read_excitation_atom, bound_excitation_atom,
                             upsilon, DBASE, KB_OVER_RY, scaled_collision_strength_interpolator)
from atomic_output import atomic_write_text
from atomic_excitation_energies import selected_energies
from niv_level_energies import STOUT_SHA256

from datetime import datetime
from zoneinfo import ZoneInfo

DATE = datetime.now(ZoneInfo("Asia/Seoul")).date().isoformat()
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "..", "data", "atomic")
TMIN, TMAX = 1.0e3, 1.0e5
TGRID = np.logspace(np.log10(TMIN), np.log10(TMAX), 201)
ERR_TARGET = 2.0e-3   # max relative error target of one transition
DEG_MAX = 12

# (element, stage, NLEV) — lowest levels needed for the optical/IR
# diagnostics. Additional levels require physical source and convergence checks.
IONS = {
    # Optional collisional H I Balmer component (par%h_coll_effects):
    # the 25-level H I atom solved with collisional excitation from the
    # ground level gives the non-recombination (collisional) line emission.
    "HI":    ("h", 1, 25),
    "OII":   ("o", 2, 5),
    "OIII":  ("o", 3, 6),
    "NII":   ("n", 2, 5),
    "SII":   ("s", 2, 5),
    "SIII":  ("s", 3, 6),
    # S IV: [S IV] 10.51um (the 3s^2 3p ^2P ground doublet) plus the
    # 3s3p^2 ^4P intercombination lines near 1400 A.
    "SIV":   ("s", 4, 5),
    "NeII":  ("ne", 2, 2),
    "NeIII": ("ne", 3, 5),
    # MOCASSIN HII20/40 line set additions ([C II] 158um, C III] 1907/1909,
    # [N III] 57um, [O I] 6300/63um)
    "CII":   ("c", 2, 5),
    "CIII":  ("c", 3, 5),
    "NIII":  ("n", 3, 5),
    "OI":    ("o", 1, 5),
    # Neutral C and N: 5 levels are exactly the 2p^2 / 2p^3 ground
    # configurations, so the FIR/optical coolants of neutral gas are complete.
    # C I: 3P_0,1,2 + 1D_2 + 1S_0 -> [C I] 609.1/370.4um fine structure,
    # 9850/9827/9811 A and 8727 A.  N I: 4S_3/2 + 2D_3/2,5/2 + 2P_1/2,3/2 ->
    # [N I] 5197.9/5200.3 A, 10398-10408 A and 3466 A.
    "CI":    ("c", 1, 5),
    "NI":    ("n", 1, 5),
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
    "HI": "H I",
    "OII": "[O II]", "OIII": "[O III]", "NII": "[N II]",
    "SII": "[S II]", "SIII": "[S III]", "SIV": "[S IV]", "NeII": "[Ne II]",
    "NeIII": "[Ne III]", "CII": "[C II]", "CIII": "C III]",
    "NIII": "[N III]", "OI": "[O I]",
    "CI": "[C I]", "NI": "[N I]",
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


def native_scale_t(tr, T):
    """Scaled temperature of the stored series, which may differ from source C."""
    if "native_cups" not in tr:
        return bt_scale_t(tr, T)
    storage = dict(tr, cups=tr["native_cups"],
                   ttype=tr.get("native_ttype", tr["ttype"]))
    return bt_scale_t(storage, T)


def native_type2_fit(tr, source_strength):
    """Represent unchanged physical Upsilon in a positive type-2 coordinate.

    The CHIANTI source C and transition type remain in tr for upsilon(tr,T).
    The native C/type metadata is used only by the Chebyshev evaluator and
    serialized row. This is attempted after the original-coordinate 16-term
    fit misses its accuracy target.
    """
    from scipy.optimize import linprog
    if np.any(source_strength <= 0) or not np.all(np.isfinite(source_strength)):
        return None
    temperatures = np.logspace(np.log10(TMIN), np.log10(TMAX), 1601)
    target = upsilon(tr, temperatures)
    if np.any(target <= 0) or not np.all(np.isfinite(target)):
        return None
    dense_temperatures = np.logspace(np.log10(TMIN), np.log10(TMAX), 160001)
    dense_target = None
    source_constant = tr["cups"]
    candidates = np.logspace(-3, 2, 11)
    candidates = sorted(candidates,
                        key=lambda value: abs(np.log(value/source_constant)))
    for constant in candidates:
        # Native decimal precision is part of the actual approximation.
        constant = float(f"{constant:.6e}")
        storage = dict(tr, native_ttype=2, native_cups=constant,
                       de=float(f"{tr['de']:.6e}"))
        ends = native_scale_t(storage, np.array([TMIN, TMAX]))
        lower, upper = (float(f"{ends[0]:.8e}"), float(f"{ends[1]:.8e}"))
        scaled = native_scale_t(storage, temperatures)
        coordinate = 2*(scaled-lower)/(upper-lower)-1
        matrix = C.chebvander(coordinate, 15)
        for logarithmic in (1, 0):
            if logarithmic:
                target_series = np.log10(target)
                linear_matrix = matrix
                rhs = np.r_[target_series, -target_series]
            else:
                magnitude = float(np.max(target))
                relative = target/magnitude
                linear_matrix = matrix/relative[:,None]
                rhs = np.r_[np.ones(len(target)), -np.ones(len(target))]
            constraints = np.vstack([
                np.column_stack([linear_matrix, -np.ones(len(target))]),
                np.column_stack([-linear_matrix, -np.ones(len(target))])])
            solution = linprog(np.r_[np.zeros(16), 1.0], A_ub=constraints,
                               b_ub=rhs, bounds=[(None,None)]*16+[(0,None)],
                               method="highs")
            if not solution.success:
                continue
            coefficients = solution.x[:16]
            if not logarithmic:
                coefficients *= magnitude
            coefficients = np.array([float(f"{value:.8e}") for value in coefficients])
            error_1601 = np.max(np.abs(eval_ups_fit(
                storage, lower, upper, coefficients, logarithmic,
                temperatures)/target-1))
            if error_1601 >= ERR_TARGET:
                continue
            if dense_target is None:
                dense_target = upsilon(tr, dense_temperatures)
            if np.any(dense_target <= 0):
                continue
            error_dense = np.max(np.abs(eval_ups_fit(
                storage, lower, upper, coefficients, logarithmic,
                dense_temperatures)/dense_target-1))
            if error_dense < ERR_TARGET:
                error = max(error_1601, error_dense,
                            _maxerr(storage, lower, upper, coefficients,
                                    logarithmic, source_strength))
                return constant, lower, upper, coefficients, logarithmic, error
    return None


def fit_transition(tr):
    """Chebyshev-fit the scaled Upsilon of one transition over TGRID.

    Returns (st_lo, st_hi, coeffs, logf, maxerr) with maxerr the relative
    error of the descaled Upsilon(T) against the direct spline evaluation.
    logf = 1 means the series stores log10 of the scaled Upsilon (used for
    weak peaked transitions whose table spans decades — a linear-space fit
    cannot hold the relative error there).
    """
    tr.pop("native_cups", None)
    tr.pop("native_ttype", None)
    st_t = bt_scale_t(tr, np.array([TMIN, TMAX]))
    st_lo = max(min(st_t), tr["xs"][0])
    st_hi = min(max(st_t), tr["xs"][-1])
    ups_ref = upsilon(tr, TGRID)
    if "temperature_grid" in tr:
        # Invert type-2 coordinate and evaluate the actual adopted table.
        def spl(scaled):
            scaled=np.asarray(scaled)
            temperature=tr["de"]/KB_OVER_RY*tr["cups"]*scaled/(1-scaled)
            return upsilon(tr, np.maximum(temperature,TMIN))
    else:
        spl = scaled_collision_strength_interpolator(tr)

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
    if best[0] >= ERR_TARGET:
        # The native reader already stores up to 16 Chebyshev coefficients.
        # Minimize the largest scaled-strength error on both temperature and
        # scaled-coordinate grids when uniform least squares misses a peak.
        # This improves numerical representation, not source resonance accuracy.
        from scipy.optimize import linprog
        samples = np.unique(np.r_[sdense, np.clip(bt_scale_t(tr,
                     np.logspace(np.log10(TMIN), np.log10(TMAX), 1601)), st_lo, st_hi)])
        coordinate = 2*(samples-st_lo)/(st_hi-st_lo)-1
        strength = spl(samples)
        for deg in range(DEG_MAX, 16):
            matrix = C.chebvander(coordinate, deg)
            for logf in (0, 1):
                if logf and np.any(strength <= 0):
                    continue
                target = np.log10(strength) if logf else strength
                scale = np.ones_like(target) if logf else np.maximum(np.abs(target), 1e-300)
                scaled_matrix = matrix/scale[:,None]
                scaled_target = target/scale
                constraints = np.vstack([
                    np.column_stack([scaled_matrix, -np.ones(len(samples))]),
                    np.column_stack([-scaled_matrix, -np.ones(len(samples))])])
                solution = linprog(np.r_[np.zeros(deg+1), 1],
                                   A_ub=constraints,
                                   b_ub=np.r_[scaled_target, -scaled_target],
                                   bounds=[(None,None)]*(deg+1)+[(0,None)], method='highs')
                if not solution.success:
                    continue
                cf = solution.x[:-1]
                err = _maxerr(tr, st_lo, st_hi, cf, logf, ups_ref)
                if err < best[0]:
                    best = (err, cf, logf)
            if best[0] < ERR_TARGET:
                break
    if best[0] >= ERR_TARGET and "temperature_grid" not in tr:
        native = native_type2_fit(tr, ups_ref)
        if native is not None:
            tr["native_cups"], st_lo, st_hi, coeff, logarithmic, error = native
            tr["native_ttype"] = 2
            return st_lo, st_hi, coeff, logarithmic, error
    return st_lo, st_hi, best[1], best[2], best[0]


def _descale(tr, T, y):
    """Descale a scaled-Upsilon value y at temperature T (BT 1992)."""
    et = (KB_OVER_RY * np.asarray(T, float)) / tr["de"]
    tt = tr.get("native_ttype", tr["ttype"])
    if tt == 1:
        return y * np.log(et + np.e)
    if tt == 2:
        return y
    if tt == 3:
        return y / (et + 1.0)
    if tt == 4:
        return y * np.log(et + tr.get("native_cups", tr["cups"]))
    if tt == 5:
        return y / et
    return 10.0 ** y


def eval_ups_fit(tr, st_lo, st_hi, coeffs, logf, T):
    """Evaluate the fitted Upsilon(T) exactly as the Fortran will."""
    st = np.clip(native_scale_t(tr, T), st_lo, st_hi)
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


def build_ion(name, elem, stage, nlev, suffix="", output_dir=OUTDIR):
    d = ion_dir(elem, stage)
    selected, radiation, collisions = bound_excitation_atom(elem, stage)
    levels = {i: selected[i] for i in sorted(selected)[:nlev]}

    def by_energy(i, j):
        """(lower-energy level, upper-energy level) of an index pair."""
        return (i, j) if levels[j][1] >= levels[i][1] else (j, i)

    arad = {}
    for ll, ul, wl, gf, a in radiation:
        if ll not in levels or ul not in levels or ll == ul or a <= 0.0:
            continue
        if elem == "n" and stage == 4 and levels[ul][1] <= levels[ll][1]:
            raise ValueError(f"N IV radiative pair {ll}->{ul} has inconsistent energy ordering")
        arad[by_energy(ll, ul)] = arad.get(by_energy(ll, ul), 0.0) + a

    ups_rows, worst = [], 0.0
    for tr in collisions:
        if tr["ll"] not in levels or tr["ul"] not in levels:
            continue
        st_lo, st_hi, cf, logf, err = fit_transition(tr)
        worst = max(worst, err)
        # Upsilon is symmetric in the two levels; the orientation matters
        # because the de-excitation rate carries g of the upper level.
        lo, up = by_energy(tr["ll"], tr["ul"])
        ups_rows.append((tr, lo, up, st_lo, st_hi, cf, logf, err))

    path = os.path.join(output_dir, f"nlevel_{elem}_{stage}{suffix}.txt")
    fh = io.StringIO()
    fh.write(f"# MoCHII n-level atom data: {SPEC_LABEL.get(name, name)}"
             f" ({elem}_{stage}), lowest {nlev} levels\n")
    if (elem,stage) == ("cl",7):
        fh.write(f"# source: Cloudy c25.00 Stout {elem}_{stage}.nrg/.tp/.coll; bound states only\n")
        fh.write("# source log-PCHIP Upsilon; type 2 is only the native storage coordinate\n")
        fh.write("# below 9800 K: constant source Upsilon, unverified resonance approximation\n")
        fh.write("# actual T enters Maxwellian rates; above-threshold core states excluded\n")
    else:
        fh.write(f"# source: CHIANTI v{CHIANTI_VERSION} .elvlc/.wgfa/.scups"
                 " (tools/fitting/fit_nlevel.py)\n")
        fh.write('# source cubic unless nonnegative nodes have negative extrema, then PCHIP\n')
    for tr, lo, up, *_ in ups_rows:
        if "native_cups" in tr:
            fh.write(f"# native {lo}->{up} type 2 C {tr['native_cups']:.6e};"
                     f" source type {tr['ttype']} C {tr['cups']:.6e};"
                     " identical physical Upsilon target\n")
    fh.write('# bound collision-connected states only; no above-threshold autoionization\n')
    if elem == "n" and stage == 4:
        fh.write(f"# energies: Cloudy c25.00 Stout n_4.nrg SHA-256 {STOUT_SHA256}\n")
        fh.write("# observed states 1--29; BSR energies above; low-T resonances uncertain\n")
        fh.write("# convergence above state 87 unproven; fit range is numerical only\n")
    elif (elem,stage) in [('o',5),('cl',6)]:
        fh.write('# energies: checked observed/ICFT vector in atomic_excitation_energies.py\n')
        fh.write('# observed O V 1--40 / Cl VI 1--44; original ICFT above\n')
        fh.write('# unshifted scattering resonances uncertain; fit interval is numerical only\n')
    elif (elem,stage) in [('o',4),('ne',5),('s',7),('ar',5),('mg',4)]:
        fh.write('# energies: observed valence complexes, original theoretical higher states\n')
        fh.write('# checked identities/term centroids; unshifted resonance accuracy uncertain\n')
    elif (elem,stage)==('s',6):
        fh.write('# energies: observed term centers with weighted theoretical fine-structure splittings\n')
    fh.write("# levels: index  g  E[cm^-1]\n")
    fh.write("# transitions: l = lower-energy level, u = upper (emitting)"
             " level; the CHIANTI column order is not used\n")
    fh.write("# rad: l u A[s^-1] (duplicate rows summed)\n")
    fh.write("# ups: l u type C dE[Ry] st_lo st_hi logf ncheb"
             " c_0..c_{n-1} maxerr\n")
    fh.write("#   logf=1: the series stores log10 of the scaled Upsilon"
             " (apply 10^y before descaling)\n")
    fh.write("#   Upsilon(T) = Burgess-Tully descaled Chebyshev series;"
             " see fit_nlevel.py header for the evaluation recipe\n")
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
    for tr, lo, up, st_lo, st_hi, cf, logf, err in ups_rows:
        cs = " ".join(f"{c:.8e}" for c in cf)
        fh.write(f"  {lo:3d} {up:3d}  {tr.get('native_ttype', tr['ttype']):1d}"
                 f" {tr.get('native_cups', tr['cups']):.6e} {tr['de']:.6e}"
                 f" {st_lo:.8e} {st_hi:.8e} {logf:1d} {len(cf):2d}  {cs}"
                 f"  {err*100:.3f}%\n")
    atomic_write_text(path, fh.getvalue())
    return path, len(levels), len(arad), len(ups_rows), worst


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ion", action="append", choices=list(IONS),
                        help="Generate a configured ion; repeat for several")
    parser.add_argument("--ion-spec", action="append", metavar="ELEMENT:STAGE:NLEV",
                        help="Generate an additional CHIANTI ion with NLEV levels")
    parser.add_argument("--output-dir", default=OUTDIR,
                        help="Generated file directory")
    parser.add_argument("--full-fe", action="store_true",
                        help="Include expanded Fe II and Fe III level models")
    args = parser.parse_args(argv)
    ions = {name: IONS[name] for name in (args.ion or ([] if args.ion_spec else IONS))}
    for spec in args.ion_spec or ():
        parts = spec.split(":")
        if len(parts) != 3 or not parts[0].isalpha() or not all(
            item.isdigit() for item in parts[1:]
        ) or min(map(int, parts[1:])) < 1:
            parser.error(f"invalid ion specification: {spec}")
        el, st, nlev = parts[0].lower(), int(parts[1]), int(parts[2])
        ions[f"{el.upper()} {st}"] = (el, st, nlev)
    print(f"{'ion':6s} {'nlev':>4s} {'nrad':>4s} {'nups':>4s}"
          f" {'worst_ups_err':>13s}")
    for name, (el, st, nlev) in ions.items():
        path, nl, nr, nu, worst = build_ion(name, el, st, nlev,
                                            output_dir=args.output_dir)
        print(f"{name:6s} {nl:4d} {nr:4d} {nu:4d} {worst*100:12.3f}%"
              f"   -> {os.path.basename(path)}")

    # Historical Fe II/III candidate files for separate source comparisons.
    # These filenames do not select the production excitation model.
    # Production selection uses coherent excitation/<model> directories.
    # Candidate sizes: Fe III 34 levels, Fe II 16 levels.
    FE_FULL = {"FeII": ("fe", 2, 16), "FeIII": ("fe", 3, 34)}
    for name, (el, st, nlev) in FE_FULL.items() if args.full_fe or (not args.ion and not args.ion_spec) else ():
        path, nl, nr, nu, worst = build_ion(name, el, st, nlev, suffix="_full",
                                            output_dir=args.output_dir)
        print(f"{name+'*':6s} {nl:4d} {nr:4d} {nu:4d} {worst*100:12.3f}%"
              f"   -> {os.path.basename(path)}")


if __name__ == "__main__":
    main()
