#!/usr/bin/env python3
"""Generate ABSOLUTE (physical-unit) spectrum files for the spectrum_type gates.

The single physical spectrum for source i is the Planck shape B_E(T_i) (same
shape as two_temp_spec.txt / the code's planck_nu) scaled by K_i so that its
ionizing-band integral int_{13.598}^{100} L_E dE = L_TARGET (1e37 erg/s).  So a
'per_ev' run with NO src_lum DERIVES src_lum ~ 1e37 and reproduces the shape run
(two_temp_file) per bin.  The SAME physical spectrum is written in every unit
system, so the 'per_ev'/'per_hz'/'per_ang'/'per_um' runs must agree:

  per_ev_2src.txt   cols  E [eV]     L_E1 L_E2   [erg/s/eV]   (E grid, ascending)
  per_ang_2src.txt  cols  lam [A]    L_l1 L_l2   [erg/s/A]    (lambda grid, asc.)
  per_um_2src.txt   cols  lam [um]   L_l1 L_l2   [erg/s/um]   (lambda grid, asc.)
  extje.txt         cols  E [eV]     J_E         [erg/s/cm^2/sr/eV]  band int = 1e-5

L_lambda = L_E * dE/dlambda = L_E * E^2 / hc.  The lambda tables use a uniform-
in-lambda grid (a natural wavelength table), so their NSUB interpolation differs
slightly from the E-grid 'per_ev' run - the residual the G-per_ang/per_um gates
bound (<0.5%).
"""
import os
import numpy as np

EV2ERG   = 1.602176634e-12      # define.f90 ev2erg
KB_CGS   = 1.380649e-16         # define.f90 kboltz_cgs
HC_EVANG = 12398.42             # define.f90 hc_evAng [eV*Angstrom]
HC_EVUM  = HC_EVANG * 1.0e-4    # [eV*micron] = 1.239842

EION_MIN = 13.598
EION_MAX = 100.0
L_TARGET = 1.0e37               # derived ionizing luminosity per source
J_TARGET = 1.0e-5               # external band-integrated J

HERE = os.path.dirname(os.path.abspath(__file__))


def planck_nu(E, T):
    """Code's planck_nu shape (arbitrary norm), E [eV], T [K]."""
    x = E * EV2ERG / (KB_CGS * T)
    return E**3 / np.expm1(x)


def ion_integral(T):
    """Dense int_{13.598}^{100} B_E(T) dE (for the absolute normalization K)."""
    Efine = np.linspace(EION_MIN, EION_MAX, 20001)
    return np.trapz(planck_nu(Efine, T), Efine)


def write_cols(path, x, cols, header):
    with open(path, "w") as f:
        f.write("# " + header + "\n")
        for i in range(x.size):
            f.write(f"{x[i]:.8e}" + "".join(f" {c[i]:.8e}" for c in cols) + "\n")


def main():
    T1, T2 = 3.0e4, 5.0e4
    K1 = L_TARGET / ion_integral(T1)
    K2 = L_TARGET / ion_integral(T2)

    # --- 'per_ev': E grid (same 400-pt log grid as two_temp_spec.txt) ---
    E = np.logspace(np.log10(10.0), np.log10(110.0), 400)
    Le1 = planck_nu(E, T1) * K1        # L_E [erg/s/eV], ionizing int = 1e37
    Le2 = planck_nu(E, T2) * K2
    write_cols(os.path.join(HERE, "per_ev_2src.txt"), E, [Le1, Le2],
               "spectrum_type=per_ev  col1 E [eV]  col2 L_E(3e4)  col3 L_E(5e4) "
               "[erg/s/eV]; ionizing int = 1e37 each")

    # --- 'per_ang' / 'per_um': uniform-in-lambda grid over the same E range ---
    lam_a = np.linspace(HC_EVANG / 110.0, HC_EVANG / 10.0, 400)   # [Angstrom] asc
    Ea = HC_EVANG / lam_a
    Lla1 = planck_nu(Ea, T1) * K1 * Ea**2 / HC_EVANG              # L_lambda [/A]
    Lla2 = planck_nu(Ea, T2) * K2 * Ea**2 / HC_EVANG
    write_cols(os.path.join(HERE, "per_ang_2src.txt"), lam_a, [Lla1, Lla2],
               "spectrum_type=per_ang  col1 lambda [A]  col2 L_lam(3e4)  "
               "col3 L_lam(5e4) [erg/s/A]")

    lam_um = np.linspace(HC_EVUM / 110.0, HC_EVUM / 10.0, 400)    # [micron] asc
    Eu = HC_EVUM / lam_um
    Llu1 = planck_nu(Eu, T1) * K1 * Eu**2 / HC_EVUM               # L_lambda [/um]
    Llu2 = planck_nu(Eu, T2) * K2 * Eu**2 / HC_EVUM
    write_cols(os.path.join(HERE, "per_um_2src.txt"), lam_um, [Llu1, Llu2],
               "spectrum_type=per_um  col1 lambda [um]  col2 L_lam(3e4)  "
               "col3 L_lam(5e4) [erg/s/um]")

    # --- external 'per_ev' J_E file: band-integrated int J_E dE = 1e-5 ---
    Kj = J_TARGET / ion_integral(4.0e4)
    Je = planck_nu(E, 4.0e4) * Kj
    write_cols(os.path.join(HERE, "extje.txt"), E, [Je],
               "spectrum_type=per_ev (external)  col1 E [eV]  col2 J_E "
               "[erg/s/cm^2/sr/eV]; ionizing int = 1e-5")

    print(f"K1={K1:.6e} K2={K2:.6e} Kj={Kj:.6e}")
    print("wrote per_ev_2src.txt, per_ang_2src.txt, per_um_2src.txt, extje.txt")


if __name__ == "__main__":
    main()
