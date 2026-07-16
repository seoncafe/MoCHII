#!/usr/bin/env python3
"""Write star_llam.txt: a 4.0e4 K blackbody in PHYSICAL units (L_lambda vs
lambda), for par%spectrum_type = 'per_ang'.

The point of this example is that the ionizing luminosity is CARRIED BY THE
FILE, not set with par%luminosity.  So the file must be an absolute spectrum
whose ionizing-segment integral equals the stromgren_sphere value
L_TARGET = 3.177837e38 erg/s (13.598-100 eV).

Construction (same Planck shape as the code's planck_nu):
    B_E(T) = E^3 / (exp(E*ev2erg/(kB T)) - 1)          [arbitrary norm]
    K      = L_TARGET / int_{13.598}^{100} B_E(T) dE    [absolute scale]
    L_E    = K * B_E(T)                                 [erg/s/eV]
    L_lam  = L_E * |dE/dlambda| = L_E * E^2 / hc        [erg/s/A]   (E = hc/lambda)
so int L_lam dlambda over the ionizing segment = int L_E dE = L_TARGET exactly.

The table is a log grid over ~90-2100 A (wavelength-ascending, as the code
accepts), spanning the 13.598-100 eV ionizing segment (123.98-911.8 A) with a
margin on both ends.  A physical-type spectrum with par%luminosity UNSET makes
MoCHII DERIVE the ionizing luminosity from this file and print it in the log
(it recovers L_TARGET to interpolation accuracy).
"""
import numpy as np

EV2ERG   = 1.602176634e-12      # define.f90 ev2erg
KB_CGS   = 1.380649e-16         # define.f90 kboltz_cgs
HC_EVANG = 12398.42             # define.f90 hc_evAng [eV*Angstrom]

EION_MIN = 13.598               # ionizing-segment edges [eV]
EION_MAX = 100.0
T_STAR   = 4.0e4                # blackbody temperature [K]
L_TARGET = 3.177837e38          # target ionizing-band luminosity [erg/s]


def planck_nu(E, T):
    """Code's planck_nu shape (arbitrary norm), E [eV], T [K]."""
    x = E * EV2ERG / (KB_CGS * T)
    return E**3 / np.expm1(x)


def ion_integral(T):
    """Dense int_{13.598}^{100} B_E(T) dE for the absolute normalization K."""
    Efine = np.linspace(EION_MIN, EION_MAX, 20001)
    return np.trapz(planck_nu(Efine, T), Efine)


def main():
    K = L_TARGET / ion_integral(T_STAR)          # absolute scale [erg/s/eV]

    lam = np.logspace(np.log10(90.0), np.log10(2100.0), 300)   # [A] ascending
    E = HC_EVANG / lam                                          # [eV]
    Le = planck_nu(E, T_STAR) * K                               # L_E [erg/s/eV]
    Llam = Le * E**2 / HC_EVANG                                 # L_lam [erg/s/A]

    np.savetxt(
        "star_llam.txt", np.c_[lam, Llam],
        header="4.0e4 K blackbody, physical units for par%%spectrum_type='per_ang'\n"
               "col1 lambda [A] (ascending)   col2 L_lambda [erg/s/A]\n"
               "ionizing-segment (13.598-100 eV) integral = %.6e erg/s" % L_TARGET,
        fmt="%.8e")
    print("wrote star_llam.txt: %d points, %.1f-%.1f A" % (lam.size, lam.min(), lam.max()))
    print("K = %.6e erg/s/eV,  ionizing integral target = %.6e erg/s" % (K, L_TARGET))


if __name__ == "__main__":
    main()
