#!/usr/bin/env python3
"""He I 2^3S (triplet metastable) diagnostic coefficients.

Writes data/atomic/hei_metastable.txt: the closed-form rate fits and the two
VFKY96 photoionization-cross-section wings for the He I 2^3S metastable state,
transcribed VERBATIM from the EXHALE triplet implementation (same author,
already verified in production).  MoCHII only evaluates the closed forms; this
script is the auditable source of every coefficient (fitted offline, evaluated
online, as the atomic-data convention requires).

The 2^3S population n_3 follows the OH18 (Oklopcic & Hirata 2018, ApJ 855 L11)
balance, LINEAR in n_3 (a closed-form quotient in each leaf):

  n_3 = [ alpha_3 n_e n_HeII + q_13 n_e n_HeI ]
        / [ A_31 + Phi_3 + n_e (q_31a + q_31b) + n_HI Q_31 ]

Sources (see the data-file header for the full provenance):
  alpha_3, alpha_1  : OH18 recombination into the triplet / singlet system.
  q_13, q_31a, q_31b: Bray et al. (2000) collision strengths via OH18/Lampon
                      (2020); q_31a = 2^3S -> 2^1S, q_31b = 2^3S -> 2^1P.
  A_31 = 1.272e-4 s^-1 (NOT the CHIANTI he_1.wgfa 1.73e-4, which is ~36% high).
  Q_31 : Taylor et al. (2025, ApJ 989 68) Table 2 Penning + associative
         ionization He(2^3S)+H0; Garcia Munoz (2025, A&A 698 A199) is recorded
         as the alternative fit (Taylor is 1.2-2.2x above GM25 over 300-1e4 K).
  sigma_3(E): two VFKY96 wings + a log-linear bridge, a ~4% fit to
              Norcross (1971) + TOPbase (EXHALE cross_sec.f90::sigma_HeI23S).

Stage 2 (10833 A absorption) adds the three fine-structure line constants
(2^3S -> 2^3P_{0,1,2}) needed for the line-center opacity kappa_0 =
sqrt(pi) e^2/(m_e c) f (lambda/b) n_3.  The wavelengths, oscillator strengths,
and Einstein A are the NIST ASD values (Drake 2007), taken through the p-winds
public OH18 implementation (p_winds/lines.py::he_3_properties) as the
cross-check.  See the file header for the NIST-vs-CHIANTI resolution.

EXHALE source routines (do not modify there):
  radiation/Cool_coeff.f90::{rec_HeII_23S, rec_HeII_11S, coex_HeI_1S_23S,
      coex_HeI_23S_21S, coex_HeI_23S_21P, penning_HeI_23S,
      coex_rate_HeI23S_10833}
  util_ion_eq.f90::HeITR_coeffs (A_31)
  functions/cross_sec.f90::sigma_HeI23S
"""
import os
import numpy as np

ATOM = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "data", "atomic")
OUT = os.path.join(ATOM, "hei_metastable.txt")

# Boltzmann constant [eV/K] and hc [eV cm], both formed with the SAME MoCHII
# constants (kboltz_cgs, ev2erg, h_planck_cgs, clight_cgs) and the SAME
# operations the Fortran module uses, so the agreement gate is bit-exact.
KB_EV = 1.380649e-16 / 1.602176634e-12                       # 8.617333e-5
HC_EVCM = 6.62607015e-27 * 2.99792458e10 / 1.602176634e-12   # 1.239841984e-4

# Rate-evaluation validity range [K]; T is clamped to this in BOTH the Fortran
# module and the Python reference so the agreement gate is exact.
T_MIN, T_MAX = 100.0, 3.0e4

# ------------------------------------------------------------------ rate fits
# alpha_i(T) = C (T/1e4)^p  [cm^3 s^-1]
ALPHA3 = (2.10e-13, -0.778)     # rec_HeII_23S
ALPHA1 = (1.54e-13, -0.486)     # rec_HeII_11S (bookkeeping only)

# q(T) = prefac sqrt(sconst/kT) exp(-Eexc/kT) (a exp(b T) + c exp(d T)) / div
# [cm^3 s^-1]; prefac=2.10e-8, sconst=13.60 eV common to the three (Bray 2000).
QEXC13  = dict(prefac=2.10e-8, sconst=13.60, Eexc=19.81, div=1.0,
               a=0.06703, b=-2.673e-6, c=-0.03584, d=-3.880e-4)   # 1^1S -> 2^3S
QEXC31A = dict(prefac=2.10e-8, sconst=13.60, Eexc=0.80,  div=3.0,
               a=2.847,   b=-1.252e-5, c=-1.953,   d=-3.558e-4)   # 2^3S -> 2^1S
QEXC31B = dict(prefac=2.10e-8, sconst=13.60, Eexc=1.40,  div=3.0,
               a=1.185,   b=-4.749e-6, c=-0.9131,  d=-1.669e-4)   # 2^3S -> 2^1P

A31 = 1.272e-4                  # s^-1  (2^3S -> 1^1S radiative; lifetime ~7870 s)

# Penning + associative ionization He(2^3S)+H0 [cm^3 s^-1]; Taylor et al. 2025.
# Q(T) = A_lo (300/T)^p_lo for T <= Tsplit, else A_hi (300/T)^p_hi.
PENNING = dict(Tsplit=4.0e3, A_lo=1.9e-9, p_lo=0.07, A_hi=9.1e-9, p_hi=0.50)
# Garcia Munoz (2025) alternative (recorded only): k = 1e-9 exp(c + d1 lnT +
# d2 lnT^2 + d3 lnT^3), c=-8.64804e1, d1=-2.86766e-1, d2=8.68445e-2,
# d3=-5.73001e-3 (<2% over 200-1e4 K).
GM25 = dict(c=-8.64804e1, d1=-2.86766e-1, d2=8.68445e-2, d3=-5.73001e-3)

# 2^3S -> 2^3P collisional excitation (10833 emission, Stage 2, recorded only):
# q(T) = 1.16e-20 sqrt(T) exp(-13179/T)  [cm^3 s^-1].
Q10833 = (1.16e-20, 13179.0)

# --- Stage 2: 10833 A fine-structure absorption constants (2^3S -> 2^3P_J).
# NIST ASD values (Drake 2007), via p-winds (p_winds/lines.py::he_3_properties).
# Vacuum wavelengths [A]; oscillator strengths f (absorption); Einstein A [s^-1].
# The three components carry EQUAL A: A_ul = (8 pi^2 e^2 / m_e c lambda^2)
# (g_l/g_u) f_lu with g_l = 3 (2^3S) and g_u = 1,3,5 (2^3P_{0,1,2}), and the
# NIST f split 1:3:5 by g_u makes (g_l/g_u) f identical, so A ~ 1.0216e7 for
# each.  The J'=1,2 pair (10833.217, 10833.306) sits 0.089 A apart, far inside
# the ~0.4 A thermal Doppler width at 1e4 K, so they blend into one observed
# feature (the "doublet"); the J'=0 line at 10832.057 is ~1.2 A away, resolved.
# RESOLUTION of the A-value discrepancy: NIST ASD / Drake (2007) give
# A(2^3S-2^3P) = 1.0216e7 s^-1 for each fine-structure component, and the public
# OH18 code p-winds adopts the same value.  The CHIANTI he_1.wgfa entry of
# ~1.15e7 s^-1 is ~13% higher and is NOT used; the NIST/p-winds value is
# adopted here.  (This is the 10833 EMISSION rate 2^3S->2^3P; it is distinct
# from A_31 = 1.272e-4 s^-1, the forbidden 2^3S->1^1S decay above.)
LINE10833 = (
    #  lambda_vac [A]   f          A [s^-1]     (J' = 0, 1, 2)
    (10832.057, 5.9902e-2, 1.0216e7),
    (10833.217, 1.7974e-1, 1.0216e7),
    (10833.306, 2.9958e-1, 1.0216e7),
)

# sigma_3(E) [Mb] VFKY96 wings [E_th, E_0, s_0, y_a, P, y_w, y_0, y_1]:
SIG3A = (4.78,  2.645,  2.08e1,  1.0e12,   3.42,  2.681, 1.956,    2.603)
SIG3B = (45.59, 49.68,  1.052e3, 4.393e-2, 2.941, 1.717, 5.488e-5, 1.118)
# bridge boundaries (EXHALE: x3=hc/357.340 A ~ 34.70 eV, x4=hc/271.940 A):
SIG3_E3 = HC_EVCM / (357.340e-8)
SIG3_E4 = HC_EVCM / (271.940e-8)


# --------------------------------------------------------- closed-form models
def alpha_power(T, C, p):
    return C * (T / 1.0e4) ** p


def qexc(T, d):
    kT = KB_EV * T
    ups = d["a"] * np.exp(d["b"] * T) + d["c"] * np.exp(d["d"] * T)
    return (d["prefac"] * np.sqrt(d["sconst"] / kT)
            * np.exp(-d["Eexc"] / kT) * ups / d["div"])


def penning(T):
    T = np.asarray(T, dtype=float)
    lo = PENNING["A_lo"] * (300.0 / T) ** PENNING["p_lo"]
    hi = PENNING["A_hi"] * (300.0 / T) ** PENNING["p_hi"]
    return np.where(T <= PENNING["Tsplit"], lo, hi)


def penning_gm25(T):
    lnT = np.log(T)
    return 1.0e-9 * np.exp(GM25["c"] + GM25["d1"] * lnT
                           + GM25["d2"] * lnT ** 2 + GM25["d3"] * lnT ** 3)


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


def n3_ratio_over_heii(T, ne, nHeI=0.0, nHI=0.0, Phi3=0.0):
    """n_3/n_HeII for the closed-form balance (source term uses n_HeII)."""
    a3 = alpha_power(T, *ALPHA3)
    q13 = qexc(T, QEXC13)
    q31a = qexc(T, QEXC31A)
    q31b = qexc(T, QEXC31B)
    Q31 = penning(T)
    # per unit n_HeII: alpha_3 ne + q_13 ne (nHeI/nHeII); use full n form outside.
    src = a3 * ne + q13 * ne * (nHeI)      # caller supplies nHeI as ratio if wanted
    snk = A31 + Phi3 + ne * (q31a + q31b) + nHI * Q31
    return src / snk


# ------------------------------------------------------------- file emission
def fmt(*vals):
    return "  ".join(("%-9s" % v[0]) if i == 0 else ("% .6e" % v)
                     for i, v in enumerate(vals))


def write_file():
    with open(OUT, "w") as f:
        w = lambda s: f.write(s + "\n")
        w("# MoCHII He I 2^3S metastable (triplet) diagnostic coefficients.")
        w("# Closed-form rate fits + VFKY96 photoionization cross section for the")
        w("# 2^3S balance (Oklopcic & Hirata 2018, ApJ 855 L11 = OH18):")
        w("#   n_3 = [alpha_3 ne n_HeII + q_13 ne n_HeI]")
        w("#       / [A_31 + Phi_3 + ne (q_31a + q_31b) + n_HI Q_31]")
        w("# Transcribed from the EXHALE triplet implementation")
        w("#   (EXHALE/src/modules/radiation/Cool_coeff.f90, util_ion_eq.f90;")
        w("#    functions/cross_sec.f90::sigma_HeI23S).  Fortran evaluates the")
        w("#   closed forms; this file is the auditable coefficient source.")
        w("# generated by tools/fitting/make_hei_metastable.py")
        w("#")
        w("# Provenance:")
        w("#  alpha_3, alpha_1 : OH18 recombination into 2^3S / 1^1S systems.")
        w("#  q_13, q_31a, q_31b : Bray et al. (2000) collision strengths via")
        w("#     OH18 / Lampon et al. (2020); q_31a = 2^3S->2^1S, q_31b = 2^3S->2^1P")
        w("#     (spin-changing e-impact into the singlets, NOT 2^3S->1^1S).")
        w("#  A_31 = 1.272e-4 s^-1 (2^3S->1^1S, lifetime ~7870 s).  PITFALL: the")
        w("#     CHIANTI he_1.wgfa A(1-2)=1.73e-4 is ~36%% high; do NOT use it.")
        w("#  Q_31 : Penning + associative ionization He(2^3S)+H0, Taylor et al.")
        w("#     (2025, ApJ 989 68) Table 2 (MB-averaged Morgner & Niehaus 1979 +")
        w("#     Cohen & Lane 1971).  Two-branch power law, DISCONTINUOUS at 4000 K")
        w("#     (jump ~1.6e-9 -> 2.5e-9).  ALTERNATIVE (recorded, not used):")
        w("#     Garcia Munoz (2025, A&A 698 A199) k = 1e-9 exp(c + d1 lnT +")
        w("#     d2 lnT^2 + d3 lnT^3), c=-8.64804e1 d1=-2.86766e-1 d2=8.68445e-2")
        w("#     d3=-5.73001e-3; Taylor is 1.2-2.2x above GM25 over 300-1e4 K -")
        w("#     the documented systematic uncertainty of the n_HI Q_31 sink.")
        w("#  sigma_3(E) : He 2^3S photoionization, two VFKY96 wings + a log-linear")
        w("#     bridge; ~4%% fit to Norcross (1971) + TOPbase (EXHALE sigma_HeI23S);")
        w("#     threshold 4.78 eV.  Cross-checked vs p-winds (same Norcross data).")
        w("#  q_10833 : 2^3S->2^3P collisional excitation (Stage 2 emission source,")
        w("#     recorded only; Black 1981 form).")
        w("#  LINE10833 : 2^3S->2^3P_{0,1,2} fine-structure line constants for the")
        w("#     10833 A line-center absorption opacity (Stage 2).  Vacuum")
        w("#     wavelength [A], oscillator strength f, Einstein A [s^-1].  NIST")
        w("#     ASD (Drake 2007) values, via p-winds (p_winds/lines.py).  A-value")
        w("#     RESOLUTION: NIST/p-winds give A=1.0216e7 s^-1 for each component; the")
        w("#     CHIANTI he_1.wgfa ~1.15e7 (~13%% high) is NOT used.  All three")
        w("#     components share A (f split 1:3:5 by g_u cancels g_l/g_u in A_ul).")
        w("#     J'=1,2 (10833.217, 10833.306) blend (0.089 A << Doppler ~0.4 A at")
        w("#     1e4 K) = the observed 'doublet'; J'=0 (10832.057) is resolved.")
        w("#")
        w("# Keyword formats (T in K, kT = KB_EV*T in eV, rates cm^3 s^-1):")
        w("#  TRANGE  Tmin Tmax          - clamp range for rate evaluation")
        w("#  A31     value              - s^-1")
        w("#  ALPHA3/ALPHA1  C p         - alpha = C (T/1e4)^p")
        w("#  QEXC13/QEXC31A/QEXC31B  prefac sconst Eexc div a b c d")
        w("#     q = prefac sqrt(sconst/kT) exp(-Eexc/kT) (a exp(bT)+c exp(dT))/div")
        w("#  PENNING Tsplit A_lo p_lo A_hi p_hi")
        w("#     Q = A_lo (300/T)^p_lo (T<=Tsplit) else A_hi (300/T)^p_hi")
        w("#  SIGMA3A/SIGMA3B  E_th E_0 s_0 y_a P y_w y_0 y_1   (VFKY96; s_0 in Mb)")
        w("#  SIGBRIDGE  E3 E4          - bridge boundaries [eV] (log-linear join)")
        w("#  Q10833  A E0              - q = A sqrt(T) exp(-E0/T)  (Stage 2)")
        w("#  LINE10833  lambda_vac_A f A   - 10833 fine-structure component")
        w("#     (one line for each 2^3P_J component; kappa_0 = sqrt(pi) e^2/(m_e c)")
        w("#      f (lambda/b) n_3, b = sqrt(2 k T_e/m_He + v_turb^2))")
        w("#")
        w("TRANGE     % .6e  % .6e" % (T_MIN, T_MAX))
        w("A31        % .6e" % A31)
        w("ALPHA3     % .6e  % .6e" % ALPHA3)
        w("ALPHA1     % .6e  % .6e" % ALPHA1)
        for key, d in (("QEXC13", QEXC13), ("QEXC31A", QEXC31A),
                       ("QEXC31B", QEXC31B)):
            w("%-10s % .6e  % .6e  % .6e  % .6e  % .6e  % .6e  % .6e  % .6e"
              % (key, d["prefac"], d["sconst"], d["Eexc"], d["div"],
                 d["a"], d["b"], d["c"], d["d"]))
        w("PENNING    % .6e  % .6e  % .6e  % .6e  % .6e"
          % (PENNING["Tsplit"], PENNING["A_lo"], PENNING["p_lo"],
             PENNING["A_hi"], PENNING["p_hi"]))
        w("SIGMA3A    " + "  ".join("% .6e" % v for v in SIG3A))
        w("SIGMA3B    " + "  ".join("% .6e" % v for v in SIG3B))
        w("SIGBRIDGE  % .6e  % .6e" % (SIG3_E3, SIG3_E4))
        w("Q10833     % .6e  % .6e" % Q10833)
        for lam, fval, a in LINE10833:
            w("LINE10833  % .6e  % .6e  % .6e" % (lam, fval, a))
    print("wrote", OUT)


# ------------------------------------------------------------- self checks
def self_check():
    print("\n=== self-checks ===")
    T = np.logspace(2, np.log10(3e4), 60)

    # (i) all rates finite and positive over the grid.
    for name, fn in (("alpha_3", lambda t: alpha_power(t, *ALPHA3)),
                     ("alpha_1", lambda t: alpha_power(t, *ALPHA1)),
                     ("q_13", lambda t: qexc(t, QEXC13)),
                     ("q_31a", lambda t: qexc(t, QEXC31A)),
                     ("q_31b", lambda t: qexc(t, QEXC31B)),
                     ("Q_31", penning)):
        v = np.asarray(fn(T))
        # finite and non-negative everywhere; the collisional rates underflow
        # to exactly 0 at low T (Boltzmann suppression, physically correct).
        ok = np.all(np.isfinite(v)) and np.all(v >= 0.0) and np.any(v > 0.0)
        print("  %-8s finite & non-negative over 1e2-3e4 K : %s" % (name, ok))
        assert ok, name

    # (ii) analytic limits at ne=100, T=1e4 K.
    T0 = 1.0e4
    a3 = alpha_power(T0, *ALPHA3)
    q31a = qexc(T0, QEXC31A)
    q31b = qexc(T0, QEXC31B)
    # low-density radiative limit n_3/n_HeII = alpha_3 ne / A_31
    lowd = a3 * 100.0 / A31
    print("  alpha_3(1e4 K)              = %.4e  (spec 2.10e-13)" % a3)
    print("  n_3/n_HeII low-dens ne=100  = %.4e  (spec 1.65e-7)" % lowd)
    assert abs(lowd / 1.65e-7 - 1.0) < 0.03, lowd
    # high-density plateau alpha_3/(q31a+q31b); crossover n_crit = A31/(q31a+q31b)
    plateau = a3 / (q31a + q31b)
    ncrit = A31 / (q31a + q31b)
    print("  high-dens plateau           = %.4e" % plateau)
    print("  n_crit = A31/(q31a+q31b)    = %.4e cm^-3 (few x 1e3 expected)" % ncrit)
    assert 1e3 < ncrit < 1e4, ncrit

    # (iii) alpha_3 vs 0.75 alpha_B(He II) consistency (diffuse_mod HEI_F3S).
    # alpha_B(He II, badnell_milne) at 1e4 K = 2.807e-13 -> 0.75x = 2.105e-13.
    print("  alpha_3(1e4)/(0.75 alpha_B_HeII=2.11e-13) = %.4f (expect ~1.00)"
          % (a3 / 2.105e-13))

    # (iv) sigma_3 sample values + p-winds cross-check.
    for E in (10.0, 21.0, 34.7, 45.6, 50.0):
        print("  sigma_3(%5.1f eV) = %.4e cm^2" % (E, sigma3(E)[0]))
    try:
        import sys
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "..", "ExoAtmosphere", "p-winds"))
        from p_winds import microphysics as mp
        wl, a3lam = mp.helium_triplet_cross_section()   # A, cm^2 (Norcross)
        Epw = HC_EVCM / (wl * 1.0e-8)
        print("  --- p-winds (Norcross) cross-check ---")
        for Etest in (7.0, 10.0, 15.0, 25.0):
            # nearest p-winds sample
            j = int(np.argmin(np.abs(Epw - Etest)))
            se = sigma3(Epw[j])[0]
            print("    E=%6.2f eV: EXHALE-fit %.3e vs p-winds %.3e (ratio %.3f)"
                  % (Epw[j], se, a3lam[j], se / a3lam[j]))
    except Exception as e:
        print("  p-winds cross-check skipped:", e)

    print("=== self-checks PASSED ===")


if __name__ == "__main__":
    write_file()
    self_check()
