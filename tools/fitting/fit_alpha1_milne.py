#!/usr/bin/env python3
"""Ground-level (n=1) radiative-recombination coefficients from the Milne relation.

WHAT THIS PRODUCES
------------------
The seven `rr_mao` coefficients of `alpha1_HII`, `alpha1_HeII` and
`alpha1_HeIII` in `src/recomb_mod.f90`.  The fit form is unchanged (Mao &
Kaastra 2016, A&A 587, A84); only the coefficients come from here.

WHY THE MILNE RELATION AND NOT A LEVEL-RESOLVED RATE FIT
--------------------------------------------------------
MoCHII forms case B as alpha_B = alpha_A - alpha_1, with alpha_A a Badnell
(2023) total-RR fit.  Taking alpha_1 from a different data set (Mao &
Kaastra's level-resolved fit) makes that difference mix two determinations:
the residual is only as accurate as the *agreement* between them, and for
He I the Mao n=1 value sits 2.4% below the Milne value, so 2.4% of alpha_1
leaks into alpha_B.

Detailed balance (the Milne relation) fixes alpha_1 exactly from the
ground-state photoionization cross section, and MoCHII already carries that
cross section: the transport absorbs with it, and the diffuse ground continua
are sampled from E^2 sigma_pi(E) exp[-(E-E_th)/kT] built from it
(`src/milne_recomb_spectrum_mod.f90`).  Deriving alpha_1 the same way makes
the opacity, the emitted spectrum and the emission rate one consistent set
instead of three independent inputs.

    alpha_1(T) = (g_n / 2 g_+) (2/sqrt(pi)) (kT)^(-3/2) sqrt(2/m_e)
                 / (m_e c^2) Int E^2 sigma_pi(E) exp[-(E-E_th)/kT] dE

CROSS SECTIONS (the transport ones, copied from src/photo_xsec.f90)
------------------------------------------------------------------
  H I  (1s)    exact hydrogenic, Z = 1  (Osterbrock & Ferland 2006, Eq. 2.4)
  He I (1^1S)  VFKY96 analytic fit (Verner, Ferland, Korista & Yakovlev 1996,
               ApJ 465, 487), the PH2 parameters of the authors' phfit2.f
  He II(1s)    exact hydrogenic, Z = 2

Statistical weights g_n/(2 g_+): H I(1s) 1.0, He I(1^1S) 0.25, He II(1s) 1.0
(the same `gfac` as `milne_setup`).

VALIDATION
----------
Cloudy c25.00 derives its H-like/He-like recombination the same way — level
by level from its own photoionization cross sections through the Milne
relation (`source/iso_radiative_recomb.cpp`, tables `data/h_iso_recomb.dat`
and `data/he_iso_recomb.dat`, whose row 0 is the ground level).  At 1e4 K:

    channel   this script    Cloudy c25.00   previous Mao coefficients
    H I       1.5790e-13     1.5840e-13      1.5873e-13
    He I      1.6086e-13     1.6032e-13      1.5650e-13
    He II     6.5116e-13     6.5145e-13      6.5366e-13

so the Milne path agrees with Cloudy's level-resolved values to <= 0.4% on
all three channels, while the Mao n=1 fit is 2.4% low on He I.

Run:  python3 fit_alpha1_milne.py
      (prints the Fortran coefficient lines, the fit residuals, and the
       comparison table; writes nothing -- the coefficients are pasted into
       src/recomb_mod.f90, whose comments carry the same provenance)
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

# --- physical constants, byte-for-byte those of src/define.f90 -------------
PI = 3.141592653589793238462643383279502884197
CLIGHT_CGS = 2.99792458e10        # cm/s
EV2ERG = 1.602176634e-12          # erg/eV
KBOLTZ_CGS = 1.380649e-16         # erg/K
ME_CGS = 9.1093837015e-28         # g
ETH_HI = 13.598                   # eV
ETH_HEI = 24.587                  # eV
ETH_HEII = 54.416                 # eV

MB2CM2 = 1.0e-18

# Upper limit of the x = (E - E_th)/kT integration, the XMAX of
# milne_recomb_spectrum_mod: exp(-60) of the near-threshold emissivity, far
# below every other error in the channel.
XMAX = 60.0

# Fitted temperature range.  Production runs bisect Te inside 3e3-3e4 K, but
# par%te_min / par%te_max allow 1e3-5e4 K, so the fit covers 1e3-1e5 K.
T_FIT_LO, T_FIT_HI = 1.0e3, 1.0e5
NT_FIT = 241

NSTART = 60        # perturbed restarts per channel (fixed seed, reproducible)
ERR_STOP = 5.0e-4  # stop restarting once the max relative residual is below this

# Cloudy c25.00 ground-level (row 0) Milne values at 1e4 K [cm^3 s^-1], from
# data/h_iso_recomb.dat and data/he_iso_recomb.dat.
CLOUDY_1E4 = {"H I": 1.5840e-13, "He I": 1.6032e-13, "He II": 6.5145e-13}

# Coefficients replaced by this script, kept for the printed comparison
# (Mao & Kaastra 2016 level fit, J/A+A/587/A84).
MAO_COEFF = {
    "H I": (2.3390e-2, 1.2310, 1.0380e-2, 1.6430e1, 7.5080e-1,
            8.9970e-2, 3.5560e-1),
    "He I": (3.7790e-2, 1.1400, 3.9760e-3, 1.2030e2, 9.9590e-1,
             3.6800, 4.1030e-1),
    "He II": (1.6140e-1, 1.1180, 1.4260e-2, 3.7370e1, 7.3750e-1,
              4.4050e-1, 3.6220e-1),
}


# --------------------------------------------------------------------------
# Photoionization cross sections [cm^2] at photon energy E [eV].
# --------------------------------------------------------------------------
def sigma_hydrogenic(E, E_th, Z):
    """Exact ground-state hydrogenic cross section, Osterbrock & Ferland Eq. 2.4."""
    s0 = 6.30e-18/(Z*Z)
    if E < E_th:
        return 0.0
    if E == E_th:
        return s0
    eps = np.sqrt(E/E_th - 1.0)
    return (s0*(E_th/E)**4*np.exp(4.0 - 4.0*np.arctan(eps)/eps)
            / (1.0 - np.exp(-2.0*PI/eps)))


def sigma_vfky96(E, E_th, E_0, s_0, y_a, P, y_w, y_0, y_1):
    """VFKY96 Eq. (1); s_0 in Mb, returns cm^2."""
    if E < E_th:
        return 0.0
    x = E/E_0 - y_0
    z = np.sqrt(x*x + y_1*y_1)
    Q = 5.5 - 0.5*P
    Fy = ((x - 1.0)**2 + y_w**2)*z**(-Q)*(1.0 + np.sqrt(z/y_a))**(-P)
    return s_0*Fy*MB2CM2


def sigma_HI(E):
    return sigma_hydrogenic(E, ETH_HI, 1.0)


def sigma_HeI(E):
    # PH2 parameters of phfit2.f for He I (2,2), as in src/photo_xsec.f90
    return sigma_vfky96(E, ETH_HEI, 1.361e1, 9.492e2, 1.469, 3.188,
                        2.039, 4.434e-1, 2.136)


def sigma_HeII(E):
    return sigma_hydrogenic(E, ETH_HEII, 2.0)


# channel -> (cross section, threshold [eV], g_n/(2 g_+))
CHANNELS = {
    "H I": (sigma_HI, ETH_HI, 1.0),          # H II   + e -> H I  (1s)
    "He I": (sigma_HeI, ETH_HEI, 0.25),      # He II  + e -> He I (1^1S)
    "He II": (sigma_HeII, ETH_HEII, 1.0),    # He III + e -> He II(1s)
}


# --------------------------------------------------------------------------
def milne_alpha1(channel, T_K):
    """alpha_1(T) [cm^3 s^-1] from the Milne relation, as milne_alpha1 does.

    The integral is taken in x = (E - E_th)/kT over [0, XMAX] and split at the
    near-threshold peak so the quadrature is limited by the cross section, not
    by the sampling of exp(-x).
    """
    sigma, E_th, gfac = CHANNELS[channel]
    kT_eV = KBOLTZ_CGS*T_K/EV2ERG
    kT_erg = KBOLTZ_CGS*T_K

    def integrand(x):
        E = E_th + x*kT_eV
        return E*E*sigma(E)*np.exp(-x)

    acc = 0.0
    edges = [0.0, 0.25, 1.0, 4.0, 12.0, 30.0, XMAX]
    for xa, xb in zip(edges[:-1], edges[1:]):
        val, _ = quad(integrand, xa, xb, epsabs=0.0, epsrel=1.0e-12, limit=400)
        acc += val
    # acc is Int E^2 sigma exp(-x) dx with E in eV; -> erg^3 cm^2
    acc *= kT_eV*EV2ERG**3
    return (gfac*(2.0/np.sqrt(PI))*kT_erg**(-1.5)
            * np.sqrt(2.0/ME_CGS)/(ME_CGS*CLIGHT_CGS**2)*acc)


# --------------------------------------------------------------------------
def rr_mao(T_K, c):
    """The Fortran rr_mao(T, a0, b0, c0, a1, b1, a2, b2); T [K], te [eV]."""
    a0, b0, c0, a1, b1, a2, b2 = c
    te = np.asarray(T_K)*(KBOLTZ_CGS/EV2ERG)
    return (a0*1.0e-10*te**(-b0 - c0*np.log(te))
            * (1.0 + a2*te**(-b2))/(1.0 + a1*te**(-b1)))


def fit_rr_mao(T, target, start):
    """Least squares in log space on the seven rr_mao coefficients."""
    logt = np.log(target)

    def resid(c):
        with np.errstate(over="ignore", invalid="ignore"):
            m = rr_mao(T, c)
        if not np.all(np.isfinite(m)) or np.any(m <= 0.0):
            return np.full_like(logt, 1.0e3)
        return np.log(m) - logt

    # a0 > 0 keeps the prefactor physical; a1, a2 >= 0 keep the rational
    # factor pole-free over the fitted range.
    lo = [1.0e-6, -5.0, -1.0, 0.0, -5.0, 0.0, -5.0]
    hi = [1.0e+3, 5.0, 1.0, 1.0e6, 5.0, 1.0e6, 5.0]
    best = None
    rng = np.random.default_rng(20260730)
    for k in range(NSTART):
        if k == 0:
            c0 = np.array(start, dtype=float)
        else:
            c0 = np.array(start, dtype=float)*np.exp(
                rng.normal(0.0, 0.45, size=7))
        c0 = np.clip(c0, lo, hi)
        try:
            sol = least_squares(resid, c0, bounds=(lo, hi), method="trf",
                                xtol=1.0e-14, ftol=1.0e-14, gtol=1.0e-14,
                                max_nfev=8000)
        except ValueError:
            continue
        err = np.max(np.abs(np.exp(resid(sol.x)) - 1.0))
        if best is None or err < best[1]:
            best = (sol.x, err)
        if best[1] < ERR_STOP:
            break
    return best


# --------------------------------------------------------------------------
def main():
    T = np.logspace(np.log10(T_FIT_LO), np.log10(T_FIT_HI), NT_FIT)
    results = {}

    print("Milne-relation alpha_1 fits (form: Mao & Kaastra 2016 rr_mao)")
    print(f"fit range {T_FIT_LO:.0e}-{T_FIT_HI:.0e} K, {NT_FIT} nodes\n")

    for name in ("H I", "He I", "He II"):
        target = np.array([milne_alpha1(name, t) for t in T])
        coeff, err = fit_rr_mao(T, target, MAO_COEFF[name])
        results[name] = (coeff, err, target)
        # residual over the production bisection range as well
        m = 3.0e3 <= T
        m &= T <= 3.0e4
        err_prod = np.max(np.abs(rr_mao(T[m], coeff)/target[m] - 1.0))
        print(f"{name:6s} max |fit/Milne - 1| = {100*err:.4f}%  "
              f"(3e3-3e4 K: {100*err_prod:.4f}%)")

    print("\nFortran coefficients (a0, b0, c0, a1, b1, a2, b2):")
    fname = {"H I": "alpha1_HII", "He I": "alpha1_HeII",
             "He II": "alpha1_HeIII"}
    for name in ("H I", "He I", "He II"):
        c = results[name][0]
        print(f"  {fname[name]:14s} " + ", ".join(f"{v:.6e}" for v in c))
    print()
    for name in ("H I", "He I", "He II"):
        c = results[name][0]
        print(f"  !--- {fname[name]}")
        print(f"    a = rr_mao(T, {c[0]:.5e}_wp, {c[1]:.5e}_wp, "
              f"{c[2]:.5e}_wp, {c[3]:.5e}_wp, &")
        print(f"               {c[4]:.5e}_wp, {c[5]:.5e}_wp, {c[6]:.5e}_wp)")

    print("\nalpha_1 at 1e4 K [cm^3 s^-1]:")
    print(f"{'channel':7s} {'new fit':>12s} {'Milne quad':>12s} "
          f"{'Cloudy c25':>12s} {'old Mao':>12s} {'fit/Cloudy':>11s} "
          f"{'Mao/Cloudy':>11s}")
    for name in ("H I", "He I", "He II"):
        c = results[name][0]
        new = float(rr_mao(1.0e4, c))
        exact = milne_alpha1(name, 1.0e4)
        cl = CLOUDY_1E4[name]
        old = float(rr_mao(1.0e4, MAO_COEFF[name]))
        print(f"{name:7s} {new:12.5e} {exact:12.5e} {cl:12.5e} "
              f"{old:12.5e} {new/cl:11.5f} {old/cl:11.5f}")

    print("\nalpha_1(T) from the new coefficients [cm^3 s^-1]:")
    print(f"{'T [K]':>9s} {'H I':>13s} {'He I':>13s} {'He II':>13s}")
    for t in (3.0e3, 5.0e3, 8.0e3, 1.0e4, 1.5e4, 2.0e4, 3.0e4):
        row = [float(rr_mao(t, results[n][0]))
               for n in ("H I", "He I", "He II")]
        print(f"{t:9.0f} " + " ".join(f"{v:13.5e}" for v in row))


if __name__ == "__main__":
    main()
