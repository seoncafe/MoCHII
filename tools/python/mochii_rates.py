"""Rate coefficients as MoCHII's Fortran evaluates them.

The analytic gates and the paper's figure generators each build a 1D reference
against which the Monte Carlo solution is judged. Those references have to use
the rates the run used, or the deviation they report is partly a difference
between two atomic datasets rather than a property of the transport. Before this
module they carried their own Hui & Gnedin (1997) case-B coefficients while the
runs used the Badnell total with the Milne ground-level rate, which at 10^4 K
differ by 0.86% for H I and by 1.2% at the 19-24 kK of the G2a thermal test,
and their free-free cooling used one Gaunt factor for every ion where the code
resolves the net charge, which at 10^4 K is 6.5% on the He III term.

Recombination mirrors `src/recomb_mod.f90` for `par%recomb_model =
'badnell_milne'`, the default:

  alpha_A  Badnell/Strathclyde TAMOC radiative recombination (`rr_badnell`),
           plus the He II dielectronic term for He I;
  alpha_1  the ground-level rate from the Milne relation applied to the
           transport cross sections, refitted onto Mao & Kaastra's functional
           form (`tools/fitting/fit_alpha1_milne.py`);
  alpha_B  alpha_A - alpha_1, exactly as `alphaB_*` computes it.

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

import numpy as np

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


#--- totals
def alphaA_HII(T):
    return rr_badnell(T, 8.318e-11, 0.7472, 2.965, 7.001e5)


def alphaA_HeII(T):
    return (rr_badnell(T, 5.235e-11, 0.6988, 7.301, 4.475e6, 0.0829, 1.682e5)
            + dr_HeII_badnell(T))


def alphaA_HeIII(T):
    return rr_badnell(T, 1.818e-10, 0.7492, 1.017e1, 2.786e6)


#--- case B, as alphaB_* computes it
def alphaB_HII(T):
    return alphaA_HII(T) - alpha1_HII(T)


def alphaB_HeII(T):
    return alphaA_HeII(T) - alpha1_HeII(T)


def alphaB_HeIII(T):
    return alphaA_HeIII(T) - alpha1_HeIII(T)


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
    return np.where(np.isclose(E, E_th), s0, out)


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
# gbar = Int g_ff(u) e^-u du over u = h nu / kT, with g_ff from the Hummer
# (1988) getGauntFF ported into src/gaunt.f90.  cooling_setup builds the table
# once on a log10 T grid for net charge Z = 1..4 and cooling_mod::gbar_ff
# interpolates linearly in log10 T; the nodes below were dumped from that
# routine on 2026-07-30 with par%gaunt_vh14 = .false. (the default), so the two
# sides agree to roundoff rather than to a refit.  Reproducing the nodes also
# carries the code's behavior verbatim, including the Z >= 3 entries below
# ~400 K where the getGauntFF extrapolation leaves the physical range; H II
# regions never go there and the free-free term needs only Z = 1 and 2.
# Regenerate whenever getGauntFF, the u quadrature, or the node grid changes;
# tests/rates/check_python_rates.py compares the two sides.
GFF_LOGT = 2.0 + 4.0 * np.arange(41) / 40.0

#--- rows: log10 T = 2.0 (0.1) 6.0;  columns: net charge Z = 1, 2, 3, 4
GFF_Z = np.array([
    [1.1740538167084194e+00, 1.0900370248407822e+00, 8.6136446328866689e-01, 3.7838339965457707e-01],   # log10 T = 2.0
    [1.1558937225974468e+00, 1.0742649273965772e+00, 8.2366600825820102e-01, 2.8062376679766171e-01],   # log10 T = 2.1
    [1.1397702360009270e+00, 1.0599617514898712e+00, 7.8588934470682148e-01, 1.7849447815848268e-01],   # log10 T = 2.2
    [1.1255351098988695e+00, 1.0471355471989454e+00, 7.4920668261938606e-01, 7.5917049545609960e-02],   # log10 T = 2.3
    [1.1130383801713117e+00, 1.0358188468065039e+00, 7.1504991633346449e-01, -2.2322379368516651e-02],   # log10 T = 2.4
    [1.1021469501741059e+00, 1.0260875477490887e+00, 6.8508527012152298e-01, -1.1070466896253932e-01],   # log10 T = 2.5
    [1.1089104309041755e+00, 1.0552090032194472e+00, 8.4028928316135121e-01, 2.9782148072813475e-01],   # log10 T = 2.6
    [1.1165418213353897e+00, 1.0723652476492969e+00, 9.4073972688734087e-01, 5.8030746783828113e-01],   # log10 T = 2.7
    [1.1250073542520087e+00, 1.0827924331572119e+00, 1.0037104323185022e+00, 7.7100782665137280e-01],   # log10 T = 2.8
    [1.1341677493514934e+00, 1.0899366554150764e+00, 1.0419243523509290e+00, 8.9628148347837422e-01],   # log10 T = 2.9
    [1.1438587062232850e+00, 1.0959384687050586e+00, 1.0645352728817949e+00, 9.7609689177367054e-01],   # log10 T = 3.0
    [1.1539403994745911e+00, 1.1020159002188443e+00, 1.0779432397893911e+00, 1.0253125944747563e+00],   # log10 T = 3.1
    [1.1643240361688501e+00, 1.1087623934832627e+00, 1.0864608131818516e+00, 1.0547515528207321e+00],   # log10 T = 3.2
    [1.1749800124533614e+00, 1.1163757794027722e+00, 1.0928565318684242e+00, 1.0720944735415683e+00],   # log10 T = 3.3
    [1.1859340662168583e+00, 1.1248253682374649e+00, 1.0987870583688888e+00, 1.0826200841320861e+00],   # log10 T = 3.4
    [1.1972553352018029e+00, 1.1339731758964036e+00, 1.1051393000390475e+00, 1.0898064647972103e+00],   # log10 T = 3.5
    [1.2090392138107444e+00, 1.1436546581616847e+00, 1.1122931011654771e+00, 1.0958171480463190e+00],   # log10 T = 3.6
    [1.2213880702909932e+00, 1.1537293964445188e+00, 1.1203189069636654e+00, 1.1018850998949570e+00],   # log10 T = 3.7
    [1.2343926335698818e+00, 1.1641074018019426e+00, 1.1291226852098715e+00, 1.1086150720892216e+00],   # log10 T = 3.8
    [1.2481159983547006e+00, 1.1747576053899718e+00, 1.1385457999356421e+00, 1.1162100795971182e+00],   # log10 T = 3.9
    [1.2625784213469977e+00, 1.1857050810530008e+00, 1.1484299669977154e+00, 1.1246436858649322e+00],   # log10 T = 4.0
    [1.2777484561443668e+00, 1.1970178599199173e+00, 1.1586573851405564e+00, 1.1337787073650130e+00],   # log10 T = 4.1
    [1.2935364165978938e+00, 1.2087910319061235e+00, 1.1691657486212015e+00, 1.1434508603765274e+00],   # log10 T = 4.2
    [1.3097930351808342e+00, 1.2211271396156174e+00, 1.1799518303325836e+00, 1.1535186183647603e+00],   # log10 T = 4.3
    [1.3263106055376741e+00, 1.2341175315002004e+00, 1.1910634622490841e+00, 1.1638908305141742e+00],   # log10 T = 4.4
    [1.3428318940378574e+00, 1.2478258147250847e+00, 1.2025836325072019e+00, 1.1745353813467323e+00],   # log10 T = 4.5
    [1.3590561811636401e+00, 1.2622730475179549e+00, 1.2146135791335710e+00, 1.1854760912030589e+00],   # log10 T = 4.6
    [1.3746533542953725e+00, 1.2774291707373930e+00, 1.2272526121813523e+00, 1.1967806319457219e+00],   # log10 T = 4.7
    [1.3892762113045560e+00, 1.2932060496471798e+00, 1.2405805842734310e+00, 1.2085433102891145e+00],   # log10 T = 4.8
    [1.4025756687149356e+00, 1.3094547212698291e+00, 1.2546425784699788e+00, 1.2208669435980772e+00],   # log10 T = 4.9
    [1.4142161612470012e+00, 1.3259694974607195e+00, 1.2694361152259903e+00, 1.2338430321202978e+00],   # log10 T = 5.0
    [1.4238895696150613e+00, 1.3424933168626414e+00, 1.2849037104412881e+00, 1.2475356597345264e+00],   # log10 T = 5.1
    [1.4313293426334639e+00, 1.3587272087836866e+00, 1.3009278702005507e+00, 1.2619680997055716e+00],   # log10 T = 5.2
    [1.4363220785400348e+00, 1.3743409678822849e+00, 1.3173332505518869e+00, 1.2771105811336618e+00],   # log10 T = 5.3
    [1.4387163836773809e+00, 1.3889871980626045e+00, 1.3338872423297541e+00, 1.2928756080428185e+00],   # log10 T = 5.4
    [1.4384317050720603e+00, 1.4023172401282207e+00, 1.3503120811820493e+00, 1.3091162082838408e+00],   # log10 T = 5.5
    [1.4354600567784255e+00, 1.4139952223912859e+00, 1.3662915358504977e+00, 1.3256283068373904e+00],   # log10 T = 5.6
    [1.4298698986923088e+00, 1.4237120727431556e+00, 1.3814862191997324e+00, 1.3421550807137530e+00],   # log10 T = 5.7
    [1.4218011529690342e+00, 1.4312001971022255e+00, 1.3955466237303407e+00, 1.3583980452437132e+00],   # log10 T = 5.8
    [1.4114599620780590e+00, 1.4362450093278916e+00, 1.4081274705083962e+00, 1.3740277440365305e+00],   # log10 T = 5.9
    [1.3991126614493543e+00, 1.4386939498447564e+00, 1.4189038764431323e+00, 1.3886976359977674e+00],   # log10 T = 6.0
])


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
def betaA_HII(T):
    """`betaA_HII` in cooling_mod: Hui & Gnedin case-A recombination cooling."""
    lam = 2.0 * T_TR_HI / np.asarray(T, dtype=float)
    return 1.778e-29 * T * lam ** 1.965 / (1.0 + (lam / 0.541) ** 0.502) ** 2.697


def betaB_HII(T):
    """`betaB_HII` in cooling_mod: the case-B counterpart."""
    lam = 2.0 * T_TR_HI / np.asarray(T, dtype=float)
    return 3.435e-30 * T * lam ** 1.970 / (1.0 + (lam / 2.250) ** 0.376) ** 3.720


def beta_HeII(T, caseA=False):
    """`beta_HeII` in cooling_mod: k_B T alpha, Hui & Gnedin's recommendation."""
    a = alphaA_HeII(T) if caseA else alphaB_HeII(T)
    return KB_CGS * np.asarray(T, dtype=float) * a


def beta_HeIII(T, caseA=False):
    """`beta_HeIII` in cooling_mod: the hydrogenic Z = 2 form, 2 x the H fit."""
    lam = 2.0 * T_TR_HEII / np.asarray(T, dtype=float)
    if caseA:
        return 2.0 * 1.778e-29 * T * lam ** 1.965 / (1.0 + (lam / 0.541) ** 0.502) ** 2.697
    return 2.0 * 3.435e-30 * T * lam ** 1.970 / (1.0 + (lam / 2.250) ** 0.376) ** 3.720


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
