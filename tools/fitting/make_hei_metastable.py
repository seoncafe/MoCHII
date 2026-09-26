#!/usr/bin/env python3
"""He I 2^3S (triplet metastable) coefficients.

Writes data/atomic/hei_metastable.txt: every coefficient of the He I 2^3S
balance except its recombination source, which the Fortran forms from the
balance's own He II case-B recombination (see below).  MoCHII evaluates what
this file carries; this script is the auditable source of every number in it
(fitted or tabulated offline, evaluated online, as the atomic-data convention
requires).

The 2^3S population n_3 is LINEAR in n_3 (a closed-form quotient in each leaf):

  n_3 = n_e [ alpha_3 n_HeII + (q_13 + q_feed) n_HeI ]
        / [ A_31 + Phi_3 + n_e (q_31g + q_31a + q_31b + q_31s + C_3)
            + n_HI Q_31 ]

Sources:
  alpha_3 : NOT in this file.  alpha_3(T) = f_3S(T) alpha_B(He II, T), the
            triplet share of the case-B cascade (hei_cascade_mod, Benjamin,
            Skillman & Smits 1999 with J. S. Mathis) times the balance's own
            case-B coefficient (recomb_mod).  Every triplet recombination ends
            in 2^3S, and a normalization the balance already fixes must not be
            supplied beside it: the Oklopcic & Hirata (2018) 2.10e-13
            T4^-0.778 the file used to carry is dropped, and so is its
            singlet companion (ALPHA1), which nothing read.
  Electron impact, 1^1S -> n^3L and 2^3S -> n^(2S+1)L, n <= 5: the effective
            collision strengths Upsilon of Bray, Burgess, Fursa & Tully
            (2000, A&AS 146, 481), Table 2 (from 1^1S, triplet columns) and
            Table 3 (from 2^3S), read from the published article (pdftotext
            of the journal PDF, checked against EXHALE's independent
            transcription: identical in every column used here; the one
            misprint of Table 2, 1^1S-4^1P at log T 5.25, is in a singlet
            column this file does not use).  Bray et al. tabulate Upsilon
            (their eq. 2), so the rate coefficient is formed here as
                q = beta Upsilon exp(-dE/kT) / (g_l sqrt T),
            beta = h^2/[(2 pi m_e)^(3/2) k^(1/2)] = 8.629132e-6 (CODATA 2018),
            g_l = 1 (1^1S) or 3 (2^3S); their eq. (1) is the same expression
            with 2 pi^(1/2) a0 hbar/m_e = 3.610e-24 m^3 s^-1 and
            Ry = 13.6058 eV.  dE: the g-weighted OBSERVED term energies of
            CHIANTI he_1.elvlc (Fuhr et al. 1999 / NIST).
            Temperature range.  Bray et al. tabulate 10^3.75-10^5.75 K and
            state why it stops there (p. 483): "Because of the uncertainty
            that pseudo resonances introduce into our cross sections at
            energies close to threshold, we limit the low temperature end of
            our tabulation to about 6 000 degrees."  They also state (Sect. 6)
            that for a neutral target Omega ~ E_j^(1/2) near threshold, so
            Upsilon must fall as T -> 0: holding the 10^3.75 K value would
            contradict their own threshold law.  Below 10^3.75 K this file
            therefore takes the TEMPERATURE DEPENDENCE of the CHIANTI v11.0.2
            he_1 strength of the same term (its J levels summed, Burgess &
            Tully descaled by chianti_cooling.upsilon), anchored to the Bray
            value at 10^3.75 K, at log T = 3.00, 3.25, 3.50; CHIANTI's
            he_1.scups header states those strengths (Sawey & Berrington
            1993, R-matrix, and Bray et al. 2000) are used over
            3.0 < log T < 5.75, so 10^3 K is where that support ends.  Where
            CHIANTI carries no transition to a term from the given lower level,
            the 10^3.75 K value is held.  The Fortran interpolates
            log10 Upsilon linearly in log10 T between nodes and holds the end
            values outside [10^3, 10^5.75] K; the Boltzmann factor is always
            evaluated at the actual T, which is what carries the temperature
            dependence where the held values stand (below 10^3 K every
            excitation here is suppressed by exp(-dE/kT) < exp(-9.2)).
            Roles (the ROLE column of each row):
              Q13     1^1S -> 2^3S, the direct excitation;
              FEED    1^1S -> n^3L, every other triplet term n <= 5: each such
                      atom decays within the triplets to 2^3S, a source of n_3
                      whose excitation energy E(n^3L) - E(2^3S) leaves as the
                      triplet cascade lines;
              Q31A    2^3S -> 2^1S;  Q31B  2^3S -> 2^1P;
              SINGLET 2^3S -> n^1L, n >= 3: a sink into the singlets;
              TRIPLET 2^3S -> n^3L (2^3P and above): the atom cascades back to
                      2^3S (10830 A for 2^3P), no population change, an energy
                      loss dE per excitation.
            q_31g (2^3S -> 1^1S) is the detailed-balance inverse of the DIRECT
            Q13 row, q_31g = q_13 (g_1/g_3) exp(E_13/kT), evaluated with the
            Boltzmann factor removed.
  C_3 : electron-impact ionization of 2^3S, Black (1981, MNRAS 197, 553)
        eq. (12), k = 8.38e-10 T^(1/2) exp(-55338/T) cm^3 s^-1, from the
        Taylor, Kingston & Bell (1979) cross sections, 1e3-1e5 K (read from
        the published scan).
  A_31 = 1.272e-4 s^-1 (NOT the CHIANTI he_1.wgfa 1.73e-4, which is ~36% high).
  Q_31 : Penning + associative ionization He(2^3S)+H0, Garcia Munoz (2025,
         A&A 698 A199, Fig. 5): the Maxwellian average of the Movre & Meyer
         (1997) total ionization cross sections, k = 1e-9 exp(c/T + d1 lnT +
         d2 lnT^2 + d3 lnT^3), <2% over 200-1e4 K; it reproduces the paper's
         Table A.5 (1.02/1.32/1.35/1.27e-9 at 500/2000/5000/1e4 K).  The same
         paper (Sect. 2) partitions it 0.9:0.1 between Penning (He + H+ + e)
         and associative (HeH+ + e) ionization, "an average"; MoCHII carries
         no HeH+, so both branches are counted as H0 -> H+ (the HeH+ ion
         carries the hydrogen's charge) and both deposit E(2^3S) - IP(H) in
         the gas, which for the associative 10% omits the HeH+ binding.
  sigma_3(E): two VFKY96 wings + a log-linear bridge, a ~4% fit to
              Norcross (1971) + TOPbase (EXHALE cross_sec.f90::sigma_HeI23S).

Stage 2 (10833 A absorption) adds the three fine-structure line constants
(2^3S -> 2^3P_{0,1,2}) needed for the line-center opacity kappa_0 =
sqrt(pi) e^2/(m_e c) f (lambda/b) n_3.  The wavelengths, oscillator strengths,
and Einstein A are the NIST ASD values (Drake 2007), taken through the p-winds
public OH18 implementation (p_winds/lines.py::he_3_properties) as the
cross-check.  See the file header for the NIST-vs-CHIANTI resolution.

The Black (1981) 2^3S -> 2^3P fit (the former Q10833 row, recorded and never
read) is superseded by the Bray et al. 2^3S -> 2^3P row and is dropped.

The module-level functions read_coefficients / triplet_rates evaluate the
FILE (not this script's constants), with the same formulas and operations as
hei_metastable_mod, for the Python gates.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ATOM = os.path.join(HERE, "..", "..", "data", "atomic")
OUT = os.path.join(ATOM, "hei_metastable.txt")

# Boltzmann constant [eV/K] and hc [eV cm], both formed with the SAME MoCHII
# constants (kboltz_cgs, ev2erg, h_planck_cgs, clight_cgs) and the SAME
# operations the Fortran module uses, so the agreement gate is bit-exact.
KB_EV = 1.380649e-16 / 1.602176634e-12                       # 8.617333e-5
HC_EVCM = 6.62607015e-27 * 2.99792458e10 / 1.602176634e-12   # 1.239841984e-4

# Clamp range [K] of the two closed forms with a stated validity range (C_3,
# Q_31); the collision strengths carry their own node grid.
# validity ranges of the two closed forms, where rates are evaluated at T
# clamped to them: Black (1981) eq. 12 is stated for 1e3-1e5 K, the
# Garcia Munoz (2025) fit to <2% over 200-1e4 K
C3_TRANGE = (1.0e3, 1.0e5)
Q31_TRANGE = (200.0, 1.0e4)

H_CGS, ME_CGS, KB_CGS = 6.62607015e-27, 9.1093837015e-28, 1.380649e-16
BETA_COLL = H_CGS**2 / ((2.0*np.pi*ME_CGS)**1.5 * np.sqrt(KB_CGS))  # 8.629132e-6

# Ionization potentials [eV] (NIST ASD), as in ion_balance_mod.
IP_HE, IP_H = 24.587387, 13.598434

# ------------------------------------------------------------ Bray et al. 2000
# Upsilon at log10 T = 3.75, 4.00, ..., 5.75, keyed (lower, n, 2S+1, L) with
# lower 1 = 1^1S (Table 2, triplet columns) and 2 = 2^3S (Table 3).
BRAY_LOGT = 3.75 + 0.25*np.arange(9)
BRAY = {
    (1, 2, 3, "S"): (0.06198, 0.06458, 0.06387, 0.06157, 0.05832, 0.0532, 0.04787, 0.04018, 0.03167),
    (1, 3, 3, "S"): (0.01642, 0.01584, 0.0153, 0.01486, 0.01434, 0.01346, 0.01198, 0.009968, 0.007741),
    (1, 4, 3, "S"): (0.005317, 0.005233, 0.005183, 0.005179, 0.005142, 0.004936, 0.004461, 0.003738, 0.002904),
    (1, 5, 3, "S"): (0.002593, 0.002542, 0.002477, 0.002448, 0.002421, 0.002325, 0.002108, 0.001774, 0.001384),
    (1, 2, 3, "P"): (0.01716, 0.02233, 0.02826, 0.03477, 0.04128, 0.04635, 0.048, 0.04503, 0.03808),
    (1, 3, 3, "P"): (0.006041, 0.006906, 0.007931, 0.009178, 0.01059, 0.01185, 0.0124, 0.01184, 0.0102),
    (1, 4, 3, "P"): (0.002074, 0.002468, 0.002983, 0.003585, 0.004236, 0.004796, 0.005048, 0.004832, 0.004169),
    (1, 5, 3, "P"): (0.001193, 0.001466, 0.001762, 0.00207, 0.002378, 0.002625, 0.002704, 0.002545, 0.002166),
    (1, 3, 3, "D"): (0.001863, 0.00206, 0.002202, 0.002322, 0.002396, 0.002372, 0.002195, 0.00186, 0.001434),
    (1, 4, 3, "D"): (0.00103, 0.001125, 0.001247, 0.001354, 0.001407, 0.001389, 0.001283, 0.001087, 0.0008368),
    (1, 5, 3, "D"): (0.0008023, 0.000863, 0.0009306, 0.0009867, 0.0009957, 0.000948, 0.0008485, 0.0007044, 0.0005371),
    (1, 4, 3, "F"): (0.0005082, 0.0004547, 0.0003748, 0.0002945, 0.0002268, 0.0001724, 0.0001279, 9.081e-05, 6.114e-05),
    (1, 5, 3, "F"): (0.0005039, 0.0004501, 0.0003664, 0.0002801, 0.0002078, 0.0001525, 0.0001111, 7.878e-05, 5.335e-05),
    (1, 5, 3, "G"): (0.0003037, 0.0002389, 0.0001717, 0.0001157, 7.485e-05, 4.734e-05, 2.968e-05, 1.86e-05, 1.197e-05),
    (2, 2, 1, "S"): (2.389, 2.456, 2.275, 1.916, 1.496, 1.111, 0.8003, 0.566, 0.3944),
    (2, 3, 1, "S"): (0.3544, 0.3295, 0.2832, 0.228, 0.173, 0.1247, 0.08624, 0.05765, 0.03747),
    (2, 4, 1, "S"): (0.1199, 0.1062, 0.08877, 0.07059, 0.05351, 0.0387, 0.02682, 0.01793, 0.01165),
    (2, 5, 1, "S"): (0.05624, 0.0493, 0.04129, 0.03305, 0.02517, 0.01823, 0.01261, 0.008402, 0.005419),
    (2, 3, 3, "S"): (2.41, 2.286, 2.235, 2.37, 2.761, 3.397, 4.187, 5.013, 5.755),
    (2, 4, 3, "S"): (0.6651, 0.5911, 0.5383, 0.5328, 0.5919, 0.7092, 0.8591, 1.016, 1.163),
    (2, 5, 3, "S"): (0.2321, 0.2154, 0.2001, 0.2009, 0.2265, 0.2732, 0.3302, 0.3889, 0.4443),
    (2, 2, 1, "P"): (0.7965, 0.9579, 1.042, 1.015, 0.895, 0.7265, 0.5516, 0.3948, 0.2677),
    (2, 3, 1, "P"): (0.1306, 0.1499, 0.1543, 0.1439, 0.1234, 0.09896, 0.07494, 0.05367, 0.03639),
    (2, 4, 1, "P"): (0.04941, 0.05791, 0.06089, 0.0574, 0.0493, 0.03939, 0.02964, 0.02106, 0.01419),
    (2, 5, 1, "P"): (0.03272, 0.03599, 0.03587, 0.03259, 0.0273, 0.02141, 0.01586, 0.01114, 0.007434),
    (2, 2, 3, "P"): (15.08, 25.8, 41.85, 64.55, 95.26, 135.1, 184.2, 240.1, 297.6),
    (2, 3, 3, "P"): (1.606, 1.611, 1.576, 1.552, 1.615, 1.868, 2.41, 3.274, 4.382),
    (2, 4, 3, "P"): (0.4476, 0.4465, 0.4476, 0.4578, 0.4938, 0.5842, 0.7603, 1.037, 1.408),
    (2, 5, 3, "P"): (0.1833, 0.1936, 0.2015, 0.2105, 0.2291, 0.2707, 0.3493, 0.4719, 0.6315),
    (2, 3, 1, "D"): (0.2479, 0.2591, 0.2558, 0.2375, 0.2071, 0.1694, 0.1292, 0.09188, 0.06137),
    (2, 4, 1, "D"): (0.08385, 0.08375, 0.08145, 0.07635, 0.06833, 0.05755, 0.04502, 0.03264, 0.02212),
    (2, 5, 1, "D"): (0.04285, 0.0424, 0.04053, 0.03781, 0.03398, 0.02874, 0.02253, 0.01635, 0.01108),
    (2, 3, 3, "D"): (1.392, 1.954, 2.76, 3.9, 5.496, 7.599, 9.99, 12.26, 14.08),
    (2, 4, 3, "D"): (0.3528, 0.4605, 0.6158, 0.8407, 1.182, 1.672, 2.262, 2.841, 3.319),
    (2, 5, 3, "D"): (0.1678, 0.2013, 0.2496, 0.3309, 0.4661, 0.6649, 0.9058, 1.143, 1.339),
    (2, 4, 1, "F"): (0.04106, 0.04019, 0.03838, 0.03475, 0.02954, 0.02341, 0.01719, 0.01175, 0.007577),
    (2, 5, 1, "F"): (0.02766, 0.02809, 0.02713, 0.0245, 0.02068, 0.01628, 0.01192, 0.008147, 0.005257),
    (2, 4, 3, "F"): (0.3325, 0.3889, 0.462, 0.5534, 0.6671, 0.7814, 0.8556, 0.8707, 0.8385),
    (2, 5, 3, "F"): (0.184, 0.2165, 0.2564, 0.3066, 0.3692, 0.4329, 0.4769, 0.4903, 0.4774),
    (2, 5, 1, "G"): (0.01409, 0.01175, 0.009224, 0.0069, 0.004932, 0.003359, 0.002181, 0.00136, 0.0008284),
    (2, 5, 3, "G"): (0.08969, 0.08533, 0.07928, 0.07392, 0.06903, 0.06279, 0.0544, 0.04498, 0.03616),
}
# extension nodes below the table: the CHIANTI shape anchored at 10^3.75 K
EXT_LOGT = np.array([3.00, 3.25, 3.50])
UPS_LOGT = np.concatenate([EXT_LOGT, BRAY_LOGT])
LNUM = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}


def role(key):
    lower, n, mult, L = key
    if lower == 1:
        return "Q13" if (n, mult, L) == (2, 3, "S") else "FEED"
    if (n, mult, L) == (2, 1, "S"):
        return "Q31A"
    if (n, mult, L) == (2, 1, "P"):
        return "Q31B"
    return "SINGLET" if mult == 1 else "TRIPLET"


# ---------------------------------------------------------------- CHIANTI he_1
def chianti_terms():
    """{(n, 2S+1, L): [level indices]} and {level: (g, E_obs_cm1)}."""
    sys.path.insert(0, HERE)
    import chianti_cooling as cc
    d = cc.ion_dir("he", 1)
    lev = cc.read_elvlc(os.path.join(d, "he_1.elvlc"))
    terms = {}
    with open(os.path.join(d, "he_1.elvlc")) as fh:
        for line in fh:
            t = line.split()
            if not t or t[0] == "-1":
                break
            cfg = t[1]
            n = int(cfg.split(".")[-1][0]) if "." in cfg else 1
            terms.setdefault((n, int(t[2]), t[3]), []).append(int(t[0]))
    trs = cc.read_scups(os.path.join(d, "he_1.scups"))
    return cc, lev, terms, trs


def term_energy_cm1(lev, terms, term):
    us = terms[term]
    g = sum(lev[u][0] for u in us)
    return sum(lev[u][0]*lev[u][1] for u in us) / g


def extension(cc, terms, trs, key):
    """Upsilon at EXT_LOGT: the CHIANTI T dependence of the term (J levels
    summed) anchored to the Bray value at 10^3.75 K; None -> held."""
    lower, n, mult, L = key
    ll = {1: 1, 2: 2}[lower]                  # CHIANTI levels of 1^1S, 2^3S
    T = 10.0**np.concatenate([EXT_LOGT, [BRAY_LOGT[0]]])
    tot, found = np.zeros_like(T), False
    for u in terms[(n, mult, L)]:
        t = [x for x in trs if x["ll"] == ll and x["ul"] == u]
        if t:
            tot = tot + cc.upsilon(t[0], T)
            found = True
    if not found or tot[-1] <= 0.0:
        return None
    return BRAY[key][0] * tot[:-1] / tot[-1]


def build_rows():
    cc, lev, terms, trs = chianti_terms()
    E23S = term_energy_cm1(lev, terms, (2, 3, "S"))
    rows = []
    for key, ups in BRAY.items():
        lower, n, mult, L = key
        E_u = term_energy_cm1(lev, terms, (n, mult, L))
        dE = (E_u - (0.0 if lower == 1 else E23S))*HC_EVCM
        ext = extension(cc, terms, trs, key)
        held = ext is None
        if held:
            ext = np.full(EXT_LOGT.size, ups[0])
        rows.append(dict(key=key, role=role(key), gl=1.0 if lower == 1 else 3.0,
                         dE=dE, ups=np.concatenate([ext, ups]), held=held))
    return rows, E23S*HC_EVCM


# 2^3S photoionization, Penning, C_3, A_31, line constants (unchanged).
CION3 = (8.38e-10, 0.5, 55338.0)
A31 = 1.272e-4
GM25 = dict(A=1.0e-9, c=-8.64804e1, d1=-2.86766e-1, d2=8.68445e-2,
            d3=-5.73001e-3)
GM25_TABLE_A5 = ((500.0, 1.02e-9), (2000.0, 1.32e-9), (5000.0, 1.35e-9),
                 (1.0e4, 1.27e-9))
LINE10833 = (
    #  lambda_vac [A]   f          A [s^-1]     (J' = 0, 1, 2)
    (10832.057, 5.9902e-2, 1.0216e7),
    (10833.217, 1.7974e-1, 1.0216e7),
    (10833.306, 2.9958e-1, 1.0216e7),
)
SIG3A = (4.78,  2.645,  2.08e1,  1.0e12,   3.42,  2.681, 1.956,    2.603)
SIG3B = (45.59, 49.68,  1.052e3, 4.393e-2, 2.941, 1.717, 5.488e-5, 1.118)
SIG3_E3 = HC_EVCM / (357.340e-8)
SIG3_E4 = HC_EVCM / (271.940e-8)


def sigma_vfky96(E, E_th, E_0, s_0, y_a, P, y_w, y_0, y_1):
    """VFKY96 full form; s_0 in Mb, returns cm^2, 0 below threshold."""
    E = np.asarray(E, dtype=float)
    x = E / E_0 - y_0
    z = np.sqrt(x * x + y_1 * y_1)
    Q = 5.5 - 0.5 * P
    Fy = ((x - 1.0) ** 2 + y_w ** 2) * z ** (-Q) * (1.0 + np.sqrt(z / y_a)) ** (-P)
    return np.where(E < E_th, 0.0, s_0 * Fy * 1.0e-18)


def sigma3(E):
    """He 2^3S photoionization [cm^2]: two VFKY96 wings + log-linear bridge."""
    E = np.atleast_1d(np.asarray(E, dtype=float))
    out = np.zeros_like(E)
    logE = np.where(E > 0.0, np.log10(np.maximum(E, 1e-300)), -300.0)
    lx3, lx4 = np.log10(SIG3_E3), np.log10(SIG3_E4)
    wa = logE <= lx3
    wb = logE >= lx4
    wm = (~wa) & (~wb)
    out[wa] = sigma_vfky96(E[wa], *SIG3A)
    out[wb] = sigma_vfky96(E[wb], *SIG3B)
    if np.any(wm):
        y3 = np.log10(sigma_vfky96(SIG3_E3, *SIG3A))
        y4 = np.log10(sigma_vfky96(SIG3_E4, *SIG3B))
        mb = (y4 - y3) / (lx4 - lx3)
        out[wm] = 10.0 ** (mb * (logE[wm] - lx3) + y3)
    return out


# ------------------------------------------------- evaluation of the FILE
def read_coefficients(path=OUT):
    """Parse the data file into a dict: scalar rows by keyword, the UPSILON
    rows as a list of dicts."""
    c = {"UPSILON": []}
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        t = ln.split()
        if t[0] == "UPSILON":
            c["UPSILON"].append(dict(role=t[1], lower=int(t[2]), n=int(t[3]),
                                     mult=int(t[4]), l=int(t[5]),
                                     gl=float(t[6]), dE=float(t[7]),
                                     ups=np.array([float(x) for x in t[8:]])))
        elif t[0] == "UPSGRID":
            c["UPSGRID"] = np.array([float(x) for x in t[2:]])
        else:
            c.setdefault(t[0], [])
            c[t[0]] = [float(x) for x in t[1:]] if t[0] != "LINE10833" \
                else c[t[0]] + [[float(x) for x in t[1:]]]
    return c


def upsilon_at(T, grid, ups):
    """log10 Upsilon linear in log10 T, end values held outside the grid."""
    x = np.clip(np.log10(T), grid[0], grid[-1])
    return 10.0**np.interp(x, grid, np.log10(ups))


def triplet_rates(T, c):
    """Rate coefficients [cm^3 s^-1] and energy-weighted sums
    [eV cm^3 s^-1] of the file at T, grouped as hei_metastable_mod does."""
    r = dict(q13=0.0, qfeed=0.0, cool13=0.0, q31a=0.0, q31b=0.0, q31s=0.0,
             q33=0.0, cool31=0.0, q31g=0.0)
    beta = c["BETA"][0]
    grid = c["UPSGRID"]
    for row in c["UPSILON"]:
        ups = upsilon_at(T, grid, row["ups"])
        q = beta*ups*np.exp(-row["dE"]/(KB_EV*T))/(row["gl"]*np.sqrt(T))
        ro = row["role"]
        if ro == "Q13":
            r["q13"] += q
            r["cool13"] += q*row["dE"]
            r["q31g"] = beta*ups/(3.0*np.sqrt(T))
            r["E13"] = row["dE"]
        elif ro == "FEED":
            r["qfeed"] += q
            r["cool13"] += q*row["dE"]
        else:
            r["cool31"] += q*row["dE"]
            key = {"Q31A": "q31a", "Q31B": "q31b", "SINGLET": "q31s",
                   "TRIPLET": "q33"}[ro]
            r[key] += q
    Tc = min(max(T, c["C3TRANGE"][0]), c["C3TRANGE"][1])
    A, p, T0 = c["CION3"]
    r["c3"] = A*Tc**p*np.exp(-T0/Tc)
    Tc = min(max(T, c["Q31TRANGE"][0]), c["Q31TRANGE"][1])
    A, cc_, d1, d2, d3 = c["PENNING"]
    lt = np.log(Tc)
    r["pen"] = A*np.exp(cc_/Tc + lt*(d1 + lt*(d2 + lt*d3)))
    r["A31"] = c["A31"][0]
    return r


def penning(T):
    T = np.asarray(T, dtype=float)
    lnT = np.log(T)
    return GM25["A"] * np.exp(GM25["c"] / T + GM25["d1"] * lnT
                              + GM25["d2"] * lnT ** 2 + GM25["d3"] * lnT ** 3)


# ------------------------------------------------------------- file emission
def write_file(rows, E23S):
    with open(OUT, "w") as f:
        w = lambda s: f.write(s + "\n")
        w("# MoCHII He I 2^3S metastable (triplet) coefficients.")
        w("# generated by tools/fitting/make_hei_metastable.py (the auditable")
        w("# source; see its docstring for every choice made here)")
        w("#")
        w("# Balance, linear in n_3 (hei_metastable_mod::hei_triplet_population):")
        w("#   n_3 = n_e [alpha_3 n_HeII + (q_13 + q_feed) n_HeI]")
        w("#       / [A_31 + Phi_3 + n_e (q_31g + q_31a + q_31b + q_31s + C_3)")
        w("#          + n_HI Q_31]")
        w("# alpha_3 is NOT in this file: alpha_3 = f_3S(T) alpha_B(He II), the")
        w("#   triplet share of the case-B cascade (hei_cascade_mod) times the")
        w("#   balance's own case-B coefficient (recomb_mod).  The earlier")
        w("#   ALPHA3/ALPHA1 rows (Oklopcic & Hirata 2018) are gone, and a file")
        w("#   carrying them is refused.")
        w("#")
        w("# Provenance:")
        w("#  UPSILON rows: effective collision strengths of Bray, Burgess, Fursa")
        w("#     & Tully (2000, A&AS 146, 481), Tables 2 (1^1S -> n^3L) and 3")
        w("#     (2^3S -> n^(2S+1)L), n <= 5, at log T = 3.75 ... 5.75 as")
        w("#     printed; at log T = 3.00, 3.25, 3.50 the temperature dependence")
        w("#     of the CHIANTI v11.0.2 he_1 strength of the same term (J levels")
        w("#     summed; Burgess & Tully descaled) anchored to the Bray value at")
        w("#     10^3.75 K, or the 10^3.75 K value held where CHIANTI carries no")
        w("#     level of the term.  Rate q = BETA Upsilon exp(-dE/kT)/(g_l sqrt T),")
        w("#     BETA = h^2/[(2 pi m_e)^1.5 k^0.5] (CODATA 2018); dE from the")
        w("#     g-weighted observed term energies of CHIANTI he_1.elvlc.")
        w("#     Upsilon: log-log linear between nodes, end values held outside.")
        w("#     q_31g (2^3S -> 1^1S) = q_13 (g_1/g_3) exp(E_13/kT) from the Q13 row.")
        w("#  C_3 : e-impact ionization of 2^3S, Black (1981, MNRAS 197, 553) eq. 12,")
        w("#     8.38e-10 T^0.5 exp(-55338/T), Taylor, Kingston & Bell (1979) cross")
        w("#     sections, 1e3-1e5 K.")
        w("#  A_31 = 1.272e-4 s^-1 (2^3S->1^1S, lifetime ~7870 s).  PITFALL: the")
        w("#     CHIANTI he_1.wgfa A(1-2)=1.73e-4 is ~36% high; do NOT use it.")
        w("#  Q_31 : Penning + associative ionization He(2^3S)+H0, Garcia Munoz")
        w("#     (2025, A&A 698 A199) Fig. 5: Maxwellian average of the Movre &")
        w("#     Meyer (1997) cross sections, k = 1e-9 exp(c/T + d1 lnT + d2 lnT^2")
        w("#     + d3 lnT^3), <2% over 200-1e4 K; reproduces the paper's Table A.5.")
        w("#     c is DIVIDED BY T (the paper's network file prints exp(c)).")
        w("#     Penning : associative = 0.9 : 0.1 (the paper's average); both")
        w("#     counted as H0 -> H+ in MoCHII, which carries no HeH+.")
        w("#  sigma_3(E) : He 2^3S photoionization, two VFKY96 wings + a log-linear")
        w("#     bridge; ~4% fit to Norcross (1971) + TOPbase (EXHALE sigma_HeI23S);")
        w("#     threshold 4.78 eV.  Cross-checked vs p-winds (same Norcross data).")
        w("#  LINE10833 : 2^3S->2^3P_{0,1,2} fine-structure line constants for the")
        w("#     10833 A line-center absorption opacity (Stage 2).  Vacuum")
        w("#     wavelength [A], oscillator strength f, Einstein A [s^-1].  NIST")
        w("#     ASD (Drake 2007) values, via p-winds (p_winds/lines.py).  A-value")
        w("#     RESOLUTION: NIST/p-winds give A=1.0216e7 s^-1 for each component; the")
        w("#     CHIANTI he_1.wgfa ~1.15e7 (~13% high) is NOT used.  All three")
        w("#     components share A (f split 1:3:5 by g_u cancels g_l/g_u in A_ul).")
        w("#     J'=1,2 (10833.217, 10833.306) blend (0.089 A << Doppler ~0.4 A at")
        w("#     1e4 K) = the observed 'doublet'; J'=0 (10832.057) is resolved.")
        w("#")
        w("# Keyword formats (T in K, kT = KB_EV*T in eV, rates cm^3 s^-1):")
        w("#  C3TRANGE  Tmin Tmax      - C_3 is evaluated at T clamped to this range")
        w("#  Q31TRANGE Tmin Tmax      - Q_31 is evaluated at T clamped to this range")
        w("#  A31     value            - s^-1")
        w("#  BETA    value            - Maxwellian constant [cm^3 s^-1 K^0.5]")
        w("#  E23S    value            - E(2^3S) above 1^1S [eV] (observed)")
        w("#  UPSGRID n log10T(1..n)   - node temperatures of the UPSILON rows")
        w("#  UPSILON role lower n 2S+1 l g_lower dE_eV Upsilon(1..n)")
        w("#     role Q13 | FEED (lower 1 = 1^1S) or Q31A | Q31B | SINGLET |")
        w("#     TRIPLET (lower 2 = 2^3S); upper term n^(2S+1)L with l = L;")
        w("#     dE = E(upper) - E(lower) [eV]")
        w("#  CION3   A p T0           - C_3 = A T^p exp(-T0/T)")
        w("#  PENNING A c d1 d2 d3")
        w("#     Q = A exp(c/T + d1 lnT + d2 (lnT)^2 + d3 (lnT)^3)")
        w("#  SIGMA3A/SIGMA3B  E_th E_0 s_0 y_a P y_w y_0 y_1   (VFKY96; s_0 in Mb)")
        w("#  SIGBRIDGE  E3 E4          - bridge boundaries [eV] (log-linear join)")
        w("#  LINE10833  lambda_vac_A f A   - 10833 fine-structure component")
        w("#     (one line for each 2^3P_J component; kappa_0 = sqrt(pi) e^2/(m_e c)")
        w("#      f (lambda/b) n_3, b = sqrt(2 k T_e/m_He + v_turb^2))")
        w("#")
        w("C3TRANGE   % .6e  % .6e" % C3_TRANGE)
        w("Q31TRANGE  % .6e  % .6e" % Q31_TRANGE)
        w("A31        % .6e" % A31)
        w("BETA       % .9e" % BETA_COLL)
        w("E23S       % .9e" % E23S)
        w("UPSGRID    %d  " % UPS_LOGT.size
          + "  ".join("%.2f" % x for x in UPS_LOGT))
        for r in rows:
            lower, n, mult, L = r["key"]
            w("#   %s -> %d^%d%s%s"
              % ("1^1S" if lower == 1 else "2^3S", n, mult, L,
                 "   (below 10^3.75 K: held, no CHIANTI transition from 2^3S)" if r["held"]
                 else ""))
            w("UPSILON  %-7s %d %d %d %d  %.1f  %.9e  " %
              (r["role"], lower, n, mult, LNUM[L], r["gl"], r["dE"])
              + "  ".join("%.6e" % u for u in r["ups"]))
        w("CION3      % .6e  % .6e  % .6e" % CION3)
        w("PENNING    % .6e  % .6e  % .6e  % .6e  % .6e"
          % (GM25["A"], GM25["c"], GM25["d1"], GM25["d2"], GM25["d3"]))
        w("SIGMA3A    " + "  ".join("% .6e" % v for v in SIG3A))
        w("SIGMA3B    " + "  ".join("% .6e" % v for v in SIG3B))
        w("SIGBRIDGE  % .6e  % .6e" % (SIG3_E3, SIG3_E4))
        for lam, fval, a in LINE10833:
            w("LINE10833  % .6e  % .6e  % .6e" % (lam, fval, a))
    print("wrote", OUT)


# ------------------------------------------------------------- self checks
def self_check(rows):
    print("\n=== self-checks ===")
    c = read_coefficients()
    # (i) the file reproduces Bray et al. at every printed node.
    worst = 0.0
    for row in c["UPSILON"]:
        key = (row["lower"], row["n"], row["mult"], "SPDFG"[row["l"]])
        for j, lt in enumerate(BRAY_LOGT):
            u = upsilon_at(10.0**lt, c["UPSGRID"], row["ups"])
            worst = max(worst, abs(u/BRAY[key][j] - 1.0))
    print("  Bray nodes reproduced from the file: worst |ratio-1| = %.1e" % worst)
    assert worst < 1e-6, worst
    held = [r["key"] for r in rows if r["held"]]
    print("  held below 10^3.75 K (no CHIANTI transition from 2^3S):", held)

    # (ii) Penning reproduces Garcia Munoz (2025) Table A.5.
    for Tk, kref in GM25_TABLE_A5:
        k = float(penning(Tk))
        assert abs(k / kref - 1.0) < 0.01, (Tk, k, kref)
    print("  Q_31 reproduces GM25 Table A.5 to 1%")

    # (iii) the rate set at 1e4 K against EXHALE's evaluation of the same
    #       Bray tables (cross-check, not a source).
    exhale = dict(q13=5.720e-19, qfeed=6.108e-20, q31a=2.804e-8, q31b=5.438e-9,
                  q31s=6.47e-10, q31g=1.858e-9)
    r = triplet_rates(1.0e4, c)
    print("  1e4 K rates [cm^3 s^-1] (MoCHII file / EXHALE):")
    for k, v in exhale.items():
        print("    %-6s %.4e  %.4e  ratio %.4f" % (k, r[k], v, r[k]/v))
    for Tk in (3.0e3, 5.0e3, 8.0e3, 1.0e4, 1.5e4, 2.0e4):
        r = triplet_rates(Tk, c)
        dsum = r["q31g"] + r["q31a"] + r["q31b"] + r["q31s"] + r["c3"]
        print("  T=%6.0f K: feed/q13 %.4f  q31s/(q31a+q31b) %.4f  "
              "n_crit = A31/(e sinks) %.0f cm^-3  triplet cooling "
              "%.3e eV cm^3/s" % (Tk, r["qfeed"]/max(r["q13"], 1e-300),
                                  r["q31s"]/(r["q31a"] + r["q31b"]),
                                  r["A31"]/dsum, r["cool31"]))
    print("=== self-checks PASSED ===")


if __name__ == "__main__":
    rows, E23S = build_rows()
    write_file(rows, E23S)
    self_check(rows)
