"""Rate coefficients as MoCHII's Fortran evaluates them.

The analytic gates and the paper's figure generators each build a 1D reference
against which the Monte Carlo solution is judged. Those references have to use
the rates the run used, or the deviation they report is partly a difference
between two atomic datasets rather than a property of the transport. Before this
module they carried their own Hui & Gnedin (1997) case-B coefficients while the
runs used the Badnell total with the Milne ground-level rate, which at 10^4 K
differed by 0.86% for H I (0.33% from the cascade alpha_B used since
2026-09-26) and by 1.2% at the 19-24 kK of the G2a thermal test,
and their free-free cooling used one Gaunt factor for every ion where the code
resolves the net charge, which at 10^4 K is 6.5% on the He III term.

Recombination mirrors `src/recomb_mod.f90` for `par%recomb_model =
'detailed_balance'`, the default (named 'badnell_milne' until 2026-09-26):

  alpha_1  the ground-level rate from the Milne relation applied to the
           transport cross sections, refitted onto Mao & Kaastra's functional
           form (`tools/fitting/fit_alpha1_milne.py`);
  He I     alpha_A the Badnell/Strathclyde TAMOC radiative recombination
           (`rr_badnell`) plus the He II dielectronic term, alpha_B =
           alpha_A - alpha_1;
  H I, He II (hydrogenic)  alpha_B from the nl-resolved cascade tabulated in
           `src/hydrogenic_recombination_data.inc` (read here, interpolated
           by the same Catmull-Rom cubic in log-log as
           `hydrogenic_recombination_mod`), alpha_A = alpha_B + alpha_1.
           Until 2026-09-26 these two took alpha_A from the Badnell fit as
           well, which put alpha_B(H II) 1.13% high at 1e4 K.

The `hui_gnedin_*` functions are kept because `par%recomb_model = 'hui_gnedin'`
selects them and one gate needs to reproduce a run made that way; they are not
the default and must not be used to build a reference for a default-mode run.

The free-free Gaunt factor and the photoionization cross sections follow, each
with its own note on what the gates had before.

Keep this file in step with `src/recomb_mod.f90`, `src/cooling_mod.f90`, and
`src/photo_xsec.f90`. `tests/milne/check_alpha1` guards the Fortran side
against drift; everything below is compared with the Fortran, group by group,
by `tests/rates/check_python_rates.py`.
"""

import os
import re

import numpy as np
from pathlib import Path

#--- cgs, matching src/define.f90
KB_CGS = 1.380649e-16
EV2ERG = 1.602176634e-12

T_TR_HI = 157807.0
T_TR_HEI = 285335.0
T_TR_HEII = 631515.0


def rr_badnell(T, A, B, T0, T1, C=0.0, T2=0.0):
    """Badnell radiative recombination fit; `rr_badnell` in recomb_mod."""
    T = np.asarray(T, dtype=float)
    bp = B + C * np.exp(-T2 / T) if C else np.full_like(T, B)
    tt = np.sqrt(T / T0)
    return A / (tt * (1.0 + tt) ** (1.0 - bp) * (1.0 + np.sqrt(T / T1)) ** (1.0 + bp))


def rr_mao(T, a0, b0, c0, a1, b1, a2, b2):
    """Mao & Kaastra (2016) functional form; `rr_mao` in recomb_mod.

    The coefficients the alpha1_* functions pass in are NOT Mao's -- they come
    from the Milne relation on the transport cross sections.
    """
    te = np.asarray(T, dtype=float) * (KB_CGS / EV2ERG)
    return (a0 * 1.0e-10 * te ** (-b0 - c0 * np.log(te))
            * (1.0 + a2 * te ** (-b2)) / (1.0 + a1 * te ** (-b1)))


def dr_HeII_badnell(T):
    """He II dielectronic recombination; `dr_HeII_badnell` in recomb_mod."""
    T = np.asarray(T, dtype=float)
    return T ** (-1.5) * (1.417e-3 * np.exp(-4.633e5 / T)
                          + 2.235e-4 * np.exp(-5.532e5 / T)
                          - 2.185e-5 * np.exp(-8.887e5 / T))


#--- ground-level rates: Milne relation on the transport cross sections
def alpha1_HII(T):
    return rr_mao(T, 1.87007e-2, 1.29077, -2.17465e-3, 1.29189e1,
                  8.05060e-1, 8.44399e-2, 5.93940e-4)


def alpha1_HeII(T):
    return rr_mao(T, 3.26862e-3, 6.97361e-1, 3.14907e-4, 5.76341e1,
                  1.18534, 2.58417e1, 9.90134e-1)


def alpha1_HeIII(T):
    return rr_mao(T, 2.61525e-1, 1.39527, -3.22079e-4, 6.40708e1,
                  8.97582e-1, 5.01551e-1, 1.21533e-3)


#--- the hydrogenic cascade table, read from the Fortran include
_HYD_INC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                        "src", "hydrogenic_recombination_data.inc")
_HYD = None


def _hydrogenic_table():
    """log T grid and log alpha_B (z = 1, 2) from the DATA statements."""
    global _HYD
    if _HYD is None:
        text = open(_HYD_INC).read().replace("&\n", " ")
        vals = {}
        for m in re.finditer(r"data\s+(\w+(?:\([^)]*\))?)\s*/(.*?)/", text,
                             re.IGNORECASE | re.DOTALL):
            nums = [float(x.replace("_wp", "").replace("E", "e"))
                    for x in m.group(2).replace(",", " ").split()]
            vals[m.group(1).replace(" ", "")] = np.array(nums)
        _HYD = (vals["hyd_logT"], vals["hyd_logaB(:,1)"], vals["hyd_logaB(:,2)"])
    return _HYD


def hydrogenic_alpha_B(z, T):
    """alpha_B of nuclear charge z; `hydrogenic_alpha_B` in
    hydrogenic_recombination_mod (Catmull-Rom in log-log, C1 across every
    node, clamped outside the table)."""
    lt_grid, la1, la2 = _hydrogenic_table()
    la = la1 if z == 1 else la2
    n = lt_grid.size
    dlt = (lt_grid[-1] - lt_grid[0]) / (n - 1)
    x = (np.clip(np.log10(np.asarray(T, dtype=float)), lt_grid[0], lt_grid[-1])
         - lt_grid[0]) / dlt
    i = np.minimum(x.astype(int), n - 2)          # 0-based node below
    w = x - i
    p1, p2 = la[i], la[i + 1]
    # at the table ends the missing neighbor is the linear extrapolation
    p0 = np.where(i > 0, la[np.maximum(i - 1, 0)], 2.0 * p1 - p2)
    p3 = np.where(i < n - 2, la[np.minimum(i + 2, n - 1)], 2.0 * p2 - p1)
    return 10.0 ** (p1 + 0.5 * w * (p2 - p0 + w * (2.0 * p0 - 5.0 * p1 + 4.0 * p2
                                                   - p3 + w * (3.0 * (p1 - p2)
                                                               + p3 - p0))))


#--- case B
def alphaB_HII(T):
    return hydrogenic_alpha_B(1, T)


def alphaB_HeII(T):
    return alphaA_HeII(T) - alpha1_HeII(T)


def alphaB_HeIII(T):
    return hydrogenic_alpha_B(2, T)


#--- totals
def alphaA_HII(T):
    return alphaB_HII(T) + alpha1_HII(T)


def alphaA_HeII(T):
    return (rr_badnell(T, 5.235e-11, 0.6988, 7.301, 4.475e6, 0.0829, 1.682e5)
            + dr_HeII_badnell(T))


def alphaA_HeIII(T):
    return alphaB_HeIII(T) + alpha1_HeIII(T)


#--- par%recomb_model = 'hui_gnedin'; not the default
def hui_gnedin_B(T, TTR, pref, lam0, p1, p2):
    lam = 2.0 * TTR / np.asarray(T, dtype=float)
    return pref * lam ** p1 / (1.0 + (lam / lam0) ** 0.407) ** p2


def hg_alphaB_HII(T):
    return hui_gnedin_B(T, T_TR_HI, 2.753e-14, 2.740, 1.500, 2.242)


def hg_alphaB_HeII(T):
    return 1.260e-14 * (2.0 * T_TR_HEI / np.asarray(T, dtype=float)) ** 0.750


def hg_alphaB_HeIII(T):
    return 2.0 * hui_gnedin_B(T, T_TR_HEII, 2.753e-14, 2.740, 1.500, 2.242)


#===========================================================================
# Photoionization cross sections, as `src/photo_xsec.f90` evaluates them.
#
# The gates carried VFKY96 fits for H I and He II where the code uses the exact
# hydrogenic expression.  The fit is 0.775% high at the H I threshold and 0.05
# to 0.5% elsewhere in the band; since the rate integral weights the threshold
# most heavily for a hot star, a reference built on the fit does not measure the
# same rate the run computed.
#===========================================================================
ETH_HI = 13.598
ETH_HEI = 24.587
ETH_HEII = 54.416


def sigma_hydrogenic(E, E_th, Z):
    """Exact hydrogenic cross section; `sigma_hydrogenic` in photo_xsec."""
    E = np.asarray(E, dtype=float)
    s0 = 6.30e-18 / (Z * Z)
    out = np.zeros_like(E)
    above = E > E_th
    eps = np.sqrt(np.where(above, E / E_th - 1.0, 1.0))
    val = (s0 * (E_th / np.where(above, E, E_th)) ** 4
           * np.exp(4.0 - 4.0 * np.arctan(eps) / eps)
           / (1.0 - np.exp(-2.0 * np.pi / eps)))
    out = np.where(above, val, out)
    return np.where(E == E_th, s0, out)


def sigma_vfky96(E, E_th, E_0, s_0, y_a, P, y_w, y_0, y_1):
    """VFKY96 Eq. (1); `sigma_vfky96` in photo_xsec."""
    E = np.asarray(E, dtype=float)
    x = E / E_0 - y_0
    z = np.sqrt(x * x + y_1 * y_1)
    Q = 5.5 - 0.5 * P
    s = s_0 * ((x - 1.0) ** 2 + y_w ** 2) * z ** (-Q) * (1.0 + np.sqrt(z / y_a)) ** (-P) * 1e-18
    return np.where(E >= E_th, s, 0.0)


def sigma_HI(E):
    """`sigma_HI` in photo_xsec: exact hydrogenic, NOT the VFKY96 fit."""
    return sigma_hydrogenic(E, ETH_HI, 1.0)


def sigma_HeI(E):
    """`sigma_HeI` in photo_xsec: VFKY96, verified against Verner's phfit2."""
    return sigma_vfky96(E, ETH_HEI, 1.361e1, 9.492e2, 1.469, 3.188,
                        2.039, 4.434e-1, 2.136)


def sigma_HeII(E):
    """`sigma_HeII` in photo_xsec: exact hydrogenic with Z = 2."""
    return sigma_hydrogenic(E, ETH_HEII, 2.0)


#===========================================================================
# Thermally averaged free-free Gaunt factor, as `src/cooling_mod.f90` has it.
#
# The gates used one Gaunt factor for every ion, where the code resolves the
# net charge: gbar(Z=2) is 6.5% below gbar(Z=1) at 10^4 K, so the He III
# free-free term of a Z-blind reference is wrong by that much.  The T
# dependence was never the problem -- the 1.1 + 0.34 exp[-(5.5-logT)^2/3]
# approximation the gates carried is within 0.06-0.6% of the code at Z = 1.
#===========================================================================
# gbar = Int_0^inf g_ff(u) exp(-u) du with the Gaunt rule of free_free_mod
# (Hummer 1988 inside its stated range, van Hoof et al. 2014 outside), on the
# grid of cooling_mod::gff_z (log10 T = 2.00 (0.01) 6.00, Z = 1..6). The nodes
# are an independent high-accuracy integral (make_gbar_ff_table.py, Gauss-
# Legendre between the kinks of the integrand, converged to 1.3e-7), not a
# copy of the Fortran: cooling_mod's 400-point trapezoid in ln u on
# [1e-9, 300] reproduces them to 8.1e-6 (Z = 1, 1e6 K; measured 2026-10-11),
# and the nodes dumped from the earlier quadrature on 0.1 dex were 3.3e-4 off.
_GFF = np.loadtxt(Path(__file__).resolve().parent/'gbar_ff_table.txt')
GFF_LOGT = _GFF[:, 0]
GFF_Z = _GFF[:, 1:]


def gbar_ff(T, Z=1):
    """Thermally averaged free-free Gaunt factor; `gbar_ff` in cooling_mod.

    Z is the net charge of the recombining ion -- 1 for H II and He II, 2 for
    He III -- and it is not optional, see the note above.
    """
    iz = min(max(int(Z), 1), GFF_Z.shape[1]) - 1
    lt = np.clip(np.log10(np.asarray(T, dtype=float)), GFF_LOGT[0], GFF_LOGT[-1])
    i = np.minimum(((lt - GFF_LOGT[0]) / (GFF_LOGT[1] - GFF_LOGT[0])).astype(int),
                   GFF_LOGT.size - 2)
    w = (lt - GFF_LOGT[i]) / (GFF_LOGT[i + 1] - GFF_LOGT[i])
    return GFF_Z[i, iz] * (1.0 - w) + GFF_Z[i + 1, iz] * w


#===========================================================================
# Recombination cooling and collisional ionization, as `src/cooling_mod.f90`
# and `src/recomb_mod.f90` evaluate them.  Exported from the Fortran for this
# purpose; see the `public` note above `betaA_HII` in cooling_mod.
#===========================================================================
DLNT_BETA = 1.0e-4


def maxwellian_recombination_cooling(alpha, T):
    """beta = k T [(3/2) alpha + d alpha / d ln T]: the kinetic energy the
    recombining electrons remove, exact for a Maxwellian and any capture cross
    section, from the same alpha the balance uses (centered difference in
    ln T, step DLNT_BETA), as `maxwellian_recombination_cooling` in
    cooling_mod."""
    T = np.asarray(T, dtype=float)
    am = alpha(T * np.exp(-DLNT_BETA))
    ap = alpha(T * np.exp(DLNT_BETA))
    return KB_CGS * T * (1.5 * alpha(T) + (ap - am) / (2.0 * DLNT_BETA))


def betaA_HII(T):
    """`betaA_HII` in cooling_mod: case-A recombination cooling of H II."""
    return maxwellian_recombination_cooling(alphaA_HII, T)


def betaB_HII(T):
    """`betaB_HII` in cooling_mod: the case-B counterpart."""
    return maxwellian_recombination_cooling(alphaB_HII, T)


def beta_HeII(T, caseA=False):
    """`beta_HeII` in cooling_mod."""
    return maxwellian_recombination_cooling(alphaA_HeII if caseA else alphaB_HeII, T)


def beta_HeIII(T, caseA=False):
    """`beta_HeIII` in cooling_mod."""
    return maxwellian_recombination_cooling(alphaA_HeIII if caseA else alphaB_HeIII, T)


#--- collisional ionization: Voronov (1997) with the optional Dere hybrid
#--- rescale.  `ci_dere_ratio` returns 1 unless par%ci_model = 'dere_hybrid',
#--- which is not the default, so these are the Voronov rates as the code uses
#--- them.  The ratio is not reproduced here; a reference built for a
#--- 'dere_hybrid' run would have to add it.
def voronov(T, dE_eV, P, A, X, Kexp):
    """`voronov` in recomb_mod."""
    U = dE_eV * EV2ERG / (KB_CGS * np.asarray(T, dtype=float))
    return A * (1.0 + P * np.sqrt(U)) * U ** Kexp * np.exp(-U) / (X + U)


def ci_HI(T):
    return voronov(T, 13.6, 0.0, 2.91e-8, 0.232, 0.39)


def ci_HeI(T):
    return voronov(T, 24.6, 0.0, 1.75e-8, 0.180, 0.35)


def ci_HeII(T):
    return voronov(T, 54.4, 1.0, 2.05e-9, 0.265, 0.25)
