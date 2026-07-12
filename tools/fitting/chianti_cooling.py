# MoCHII: copied from EXHALE/cooling_data/chianti_cooling.py (2026-07-11); origin tree is read-only reference
"""Minimal, dependency-light CHIANTI reader + collisional line-cooling builder.

Reads the raw CHIANTI ASCII ion files (.elvlc, .wgfa, .scups) directly --- no
ChiantiPy required --- and computes the optically-thin collisional line-cooling
function

    Lambda_ion(T) = sum_{l->u}  q_{lu}(T) * dE_{lu}        [erg cm^3 s^-1]

i.e. cooling per (n_e * n_lower-level) assuming every collisional excitation is
followed by a radiative decay (coronal / low-density limit). The Maxwellian
excitation-rate coefficient is

    q_{lu}(T) = 8.629e-6 / (g_l sqrt(T)) * Upsilon_{lu}(T) * exp(-dE_{lu}/kT)

with the effective collision strength Upsilon obtained by Burgess & Tully (1992)
descaling of the scaled values stored in the .scups file.

Provenance: CHIANTI v11 database at $XUVTOP (here
/home/kiseon/RT_Codes/CHIANTI/dbase). This script is the auditable source for
every cooling coefficient that gets ported into the Fortran (Phase 2).
"""

import os
import numpy as np
from scipy.interpolate import CubicSpline

# Physical constants (CODATA / NIST)
K_B_ERG = 1.380649e-16        # erg / K
RY_ERG = 2.1798723611035e-11  # 1 Rydberg in erg
EV_ERG = 1.602176634e-12      # erg / eV
KB_OVER_RY = K_B_ERG / RY_ERG  # Rydberg per Kelvin = 6.33363e-6
COLL_PREF = 8.629e-6          # cm^3 s^-1 K^1/2 (Maxwellian rate prefactor)

DBASE = "/home/kiseon/RT_Codes/CHIANTI/dbase"


def ion_dir(elem, ion):
    """CHIANTI directory for e.g. ('fe', 2) -> .../fe/fe_2."""
    return os.path.join(DBASE, elem, f"{elem}_{ion}")


def read_elvlc(path):
    """Return {level_index: (g, E_obs_cm1, E_th_cm1)}.

    Robust to configuration labels containing spaces: every CHIANTI elvlc data
    line ends with  ... J  E_obs  E_th , so J is the 3rd-from-last token and the
    statistical weight is g = 2J + 1.
    """
    levels = {}
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n")
            if s.strip() in ("-1", ""):
                break
            toks = s.split()
            if len(toks) < 4:
                continue
            try:
                idx = int(toks[0])
                e_th = float(toks[-1])
                e_obs = float(toks[-2])
                jval = float(toks[-3])
            except ValueError:
                continue
            levels[idx] = (2.0 * jval + 1.0, e_obs, e_th)
    return levels


def read_wgfa(path):
    """Return list of (lower, upper, wavelength_A, gf, A_s^-1)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n")
            if s.strip() in ("-1", ""):
                break
            toks = s.split()
            if len(toks) < 5:
                continue
            try:
                ll = int(toks[0]); ul = int(toks[1])
                wl = float(toks[2]); gf = float(toks[3]); a = float(toks[4])
            except ValueError:
                continue
            rows.append((ll, ul, wl, gf, a))
    return rows


def read_scups(path):
    """Parse the CHIANTI .scups file (scaled effective collision strengths).

    Returns list of dicts with keys:
      ll, ul       lower / upper level index
      de           transition energy [Rydberg]
      gf           weighted oscillator strength
      lim          high-T limit (-1 if unused)
      ttype        Burgess-Tully transition type (1..6)
      cups         scaling parameter C
      xs, ys       scaled temperature and scaled effective collision strength
    """
    trans = []
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    i = 0
    n = len(lines)
    while i < n:
        head = lines[i].split()
        if len(head) < 8 or head[0] == "-1":
            i += 1
            continue
        try:
            ll = int(head[0]); ul = int(head[1])
            de = float(head[2]); gf = float(head[3]); lim = float(head[4])
            nt = int(head[5]); ttype = int(head[6]); cups = float(head[7])
        except (ValueError, IndexError):
            i += 1
            continue
        xs = np.array(lines[i + 1].split(), dtype=float)
        ys = np.array(lines[i + 2].split(), dtype=float)
        if xs.size == nt and ys.size == nt:
            trans.append(dict(ll=ll, ul=ul, de=de, gf=gf, lim=lim,
                              ttype=ttype, cups=cups, xs=xs, ys=ys))
        i += 3
    return trans


def upsilon(tr, T):
    """Burgess & Tully (1992) descaled effective collision strength Upsilon(T).

    T may be a scalar or array [K]. de is in Rydberg.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    de = tr["de"]
    C = tr["cups"]
    tt = tr["ttype"]
    et = (KB_OVER_RY * T) / de            # kT / dE  (dimensionless)
    if tt in (1, 4):
        st = 1.0 - np.log(C) / np.log(et + C)
    else:
        st = et / (et + C)
    spl = CubicSpline(tr["xs"], tr["ys"])
    sups = spl(np.clip(st, tr["xs"][0], tr["xs"][-1]))
    if tt == 1:
        ups = sups * np.log(et + np.e)
    elif tt == 2:
        ups = sups
    elif tt == 3:
        ups = sups / (et + 1.0)
    elif tt == 4:
        ups = sups * np.log(et + C)
    elif tt == 5:
        ups = sups / et
    elif tt == 6:
        ups = 10.0 ** sups
    else:
        raise ValueError(f"unknown BT type {tt}")
    return np.maximum(ups, 0.0)


def cooling_lambda(elem, ion, T, lower_levels=None, upper_levels=None):
    """Optically-thin collisional line-cooling Lambda(T) [erg cm^3 s^-1].

    Sum over .scups transitions of  q_lu * dE_lu, q_lu the Maxwellian
    excitation-rate coefficient out of the lower level. Optionally restrict to a
    set of lower and/or upper levels (e.g. {1} ground, {3,4} the 2p Ly-alpha
    doublet for H I).
    Returns (Lambda_total, per_transition) where per_transition is a list of
    (ll, ul, Lambda_array).
    """
    d = ion_dir(elem, ion)
    lev = read_elvlc(os.path.join(d, f"{elem}_{ion}.elvlc"))
    trans = read_scups(os.path.join(d, f"{elem}_{ion}.scups"))
    T = np.atleast_1d(np.asarray(T, dtype=float))
    total = np.zeros_like(T)
    per = []
    for tr in trans:
        if lower_levels is not None and tr["ll"] not in lower_levels:
            continue
        if upper_levels is not None and tr["ul"] not in upper_levels:
            continue
        g_l = lev[tr["ll"]][0]
        de_erg = tr["de"] * RY_ERG
        ups = upsilon(tr, T)
        q_lu = COLL_PREF / (g_l * np.sqrt(T)) * ups * np.exp(-de_erg / (K_B_ERG * T))
        lam = q_lu * de_erg
        total += lam
        per.append((tr["ll"], tr["ul"], lam))
    return total, per


def black1981_HI_lya(T):
    """Black (1981) H I Ly-alpha collisional-excitation cooling currently in
    Cool_coeff.f90 (coex_rate_HI): 7.5e-19/(1+sqrt(T/1e5)) exp(-118348/T)."""
    T = np.asarray(T, dtype=float)
    return 7.5e-19 / (1.0 + np.sqrt(T / 1.0e5)) * np.exp(-118348.0 / T)


HC_OVER_K = 1.438776877  # h c / k_B  [cm K]; E[cm^-1]*HC_OVER_K/T = E/kT


def _level_energy_cm1(lev):
    """Map {idx:(g,E_obs,E_th)} -> {idx:(g, E_cm1)} using observed energy when
    available (>0 or ground), else the theoretical value."""
    out = {}
    for idx, (g, e_obs, e_th) in lev.items():
        e = e_obs if (e_obs > 0.0 or idx == 1) else e_th
        out[idx] = (g, e)
    return out


def cooling_effective(elem, ion, T, pop="coronal", e_cut_cm1=None,
                      restrict_to_wgfa=False):
    """Effective collisional line-cooling per (n_e * n_ion_total) [erg cm^3 s^-1].

    Sum over .scups transitions of  f_l(T) * q_lu(T) * dE_lu, where f_l is the
    fractional population of the lower level l and q_lu the Maxwellian
    excitation-rate coefficient. Upper levels are assumed to decay radiatively
    (optically thin), so every excitation out of a populated lower level cools.

    pop selects the lower-level population model:
      'coronal'     only the ground level (idx 1) is populated  (f_1 = 1).
      'ground_term' Boltzmann among the lowest levels sharing the ground config
                    label (e.g. the Fe II a6D fine-structure multiplet).
      'boltzmann'   Boltzmann over every level with E < e_cut_cm1 (the
                    metastable manifold); levels >= e_cut are depleted. This is
                    Huang (2023) Section 2.5 for Fe I/Fe II (cut at the lowest
                    odd-parity level that has an E1 decay to the ground term).

    restrict_to_wgfa keeps only transitions that also appear in the .wgfa file
    (i.e. carry a radiative A / oscillator strength), matching Huang's "1105 of
    4339" Fe II selection.
    """
    d = ion_dir(elem, ion)
    lev = _level_energy_cm1(read_elvlc(os.path.join(d, f"{elem}_{ion}.elvlc")))
    trans = read_scups(os.path.join(d, f"{elem}_{ion}.scups"))
    T = np.atleast_1d(np.asarray(T, dtype=float))

    wgfa_pairs = None
    if restrict_to_wgfa:
        rows = read_wgfa(os.path.join(d, f"{elem}_{ion}.wgfa"))
        wgfa_pairs = {(ll, ul) for (ll, ul, wl, gf, a) in rows
                      if (gf != 0.0 or a != 0.0)}

    # build lower-level population fractions f_l(T) (each an array over T)
    if pop == "coronal":
        lows = [1]
    elif pop == "ground_term":
        # Approximate the ground multiplet by the lowest levels within ~0.4 eV
        # of ground (e.g. the Fe II a6D fine-structure block at 0-977 cm^-1).
        thr_cm1 = 0.4 * EV_ERG / (HC_OVER_K * K_B_ERG)  # 0.4 eV -> ~3226 cm^-1
        lows = [i for i, (g, e) in lev.items() if e <= thr_cm1]
    elif pop == "boltzmann":
        if e_cut_cm1 is None:
            raise ValueError("pop='boltzmann' requires e_cut_cm1")
        lows = [i for i, (g, e) in lev.items() if e < e_cut_cm1]
    else:
        raise ValueError(f"unknown pop model {pop!r}")

    if pop == "coronal":
        f = {1: np.ones_like(T)}
    else:
        Z = np.zeros_like(T)
        for i in lows:
            g, e = lev[i]
            Z += g * np.exp(-e * HC_OVER_K / T)
        f = {i: lev[i][0] * np.exp(-lev[i][1] * HC_OVER_K / T) / Z for i in lows}

    lowset = set(lows)
    total = np.zeros_like(T)
    for tr in trans:
        ll = tr["ll"]
        if ll not in lowset:
            continue
        if wgfa_pairs is not None and (ll, tr["ul"]) not in wgfa_pairs:
            continue
        g_l = lev[ll][0]
        de_erg = tr["de"] * RY_ERG
        ups = upsilon(tr, T)
        q_lu = COLL_PREF / (g_l * np.sqrt(T)) * ups * np.exp(-de_erg / (K_B_ERG * T))
        total += f[ll] * q_lu * de_erg
    return total


def free_free_cooling(T, Z=1.0):
    """Huang (2023) Eq. 9 free-free cooling per (n_e * n_ion) [erg cm^3 s^-1]:
    Lambda_ff = 1.9095e-25 * Z^2 * T4^0.55,  T4 = T/1e4 K."""
    T = np.asarray(T, dtype=float)
    return 1.9095e-25 * Z**2 * (T / 1.0e4) ** 0.55


EIGHTPI_RT3 = 8.0 * np.pi / np.sqrt(3.0)  # 14.5104, Van Regemorter prefactor
RY_EV = 13.605693                          # Rydberg [eV] = hydrogen ionization I_H


def van_regemorter_clu_compact(T, dE_eV, f):
    """Compact Van Regemorter form quoted by Huang (2023) Eq. 10:

        C_lu(e) = 2.16 * alpha^-1.68 * exp(-alpha) * T^-3/2 * f ,  alpha = dE/kT.

    WARNING: cross-checking this exact expression against CHIANTI for the Mg II
    and Ca II resonance doublets shows it underestimates the electron-impact
    excitation rate by 3-8x (worse at low T; see cooling_data docs). Huang
    applies it only to the ~1025 permitted Fe I lines that have no CHIANTI
    collision data. Kept for that Fe I use and for documentation; for strong
    resonance lines prefer van_regemorter_rate().
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    alpha = (dE_eV * EV_ERG) / (K_B_ERG * T)
    return 2.16 * alpha ** (-1.68) * np.exp(-alpha) * T ** (-1.5) * f


def van_regemorter_rate(T, dE_eV, gf, g_l, gbar=0.2):
    """Rigorous Van Regemorter (1962) excitation rate coefficient [cm^3 s^-1]:

        Upsilon = (8 pi/sqrt 3) * gf * (I_H/dE) * gbar
        C_lu    = 8.629e-6 / (g_l sqrt(T)) * Upsilon * exp(-dE/kT)

    gf = g_l * f_absorption is the weighted oscillator strength; gbar is the
    effective (Maxwellian-averaged) Gaunt factor. gbar ~ 0.2 for ions: this
    reproduces the CHIANTI effective collision strengths of the Mg II and Ca II
    resonance doublets to within ~20% over 5000-30000 K (back-solved gbar =
    0.19-0.24 for Mg II, 0.14-0.22 for Ca II). Used for Na I (isoelectronic with
    Mg II, no CHIANTI data) and, optionally, Fe I.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    alpha = (dE_eV * EV_ERG) / (K_B_ERG * T)
    ups = EIGHTPI_RT3 * gf * (RY_EV / dE_eV) * gbar
    return COLL_PREF / (g_l * np.sqrt(T)) * ups * np.exp(-alpha)


# Na I resonance D doublet (3s 2S1/2 -> 3p 2P_J); NIST A/f values (Kramida 2020).
#   (dE_eV, f_absorption, label).  Lower level 3s 2S1/2 has g_l = 2.
NAI_GL = 2.0
NaI_D_LINES = [
    (2.10230, 0.320, "D1 5895.92A 3s2S1/2-3p2P1/2"),
    (2.10440, 0.641, "D2 5889.95A 3s2S1/2-3p2P3/2"),
]


def cooling_NaI(T, gbar=0.2):
    """Na I D-doublet collisional line cooling per (n_e * n_NaI) [erg cm^3 s^-1].

    Coronal/optically-thin limit: every electron-impact excitation out of the
    3s ground level is followed by a radiative decay, so Lambda = sum_D C_lu*dE.
    Na I is absent from CHIANTI v11.0.2, but it is isoelectronic with Mg II
    (Na-like 3s 2S -> 3p 2P resonance doublet), so the collision rate uses the
    rigorous Van Regemorter form with the Gaunt factor gbar calibrated against
    CHIANTI Mg II/Ca II (~0.2). This yields Upsilon(3s-3p) ~ 36 near 1e4 K,
    consistent with published Na I R-matrix collision strengths; the anomalous
    two-electron Mg I 1S-1P line (gf=3.4 but Upsilon~2) is NOT a valid analog.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    total = np.zeros_like(T)
    for dE_eV, f, _ in NaI_D_LINES:
        gf = NAI_GL * f
        total += van_regemorter_rate(T, dE_eV, gf, NAI_GL, gbar) * (dE_eV * EV_ERG)
    return total


if __name__ == "__main__":
    T = np.logspace(3.0, 5.0, 60)
    lya, per = cooling_lambda("h", 1, T, lower_levels={1}, upper_levels={3, 4})
    tot, _ = cooling_lambda("h", 1, T, lower_levels={1})
    blk = black1981_HI_lya(T)
    print(f"{'T[K]':>10} {'CHIANTI_Lya':>13} {'CHIANTI_tot':>13} "
          f"{'Black':>12} {'Black/CHIANTI_Lya':>18}")
    for k in range(0, len(T), 4):
        print(f"{T[k]:10.3e} {lya[k]:13.4e} {tot[k]:13.4e} "
              f"{blk[k]:12.4e} {blk[k]/lya[k]:18.2f}")
    # Cross-check the descaled Upsilon(1s-2p) against the pre-extracted file.
    d = ion_dir("h", 1)
    trans = read_scups(os.path.join(d, "h_1.scups"))
    t13 = [t for t in trans if t["ll"] == 1 and t["ul"] == 3][0]
    t14 = [t for t in trans if t["ll"] == 1 and t["ul"] == 4][0]
    for Tc in (1.0e4, 2.0e4):
        u = upsilon(t13, Tc)[0] + upsilon(t14, Tc)[0]
        print(f"Upsilon(1s-2p) at T={Tc:.0e} K : {u:.4f}")
