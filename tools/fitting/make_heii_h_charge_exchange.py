#!/usr/bin/env python3
"""Write src/heii_h_charge_exchange_data.inc: the non-radiative rate
coefficient of

    He+(1s 2S) + H(1s 2S) -> He(1s2 1S) + H+ ,

from the reverse-direction calculation of Loreau, Ryabchenko & Vaeck (2014,
J. Phys. B 47, 135204; read), taken through microscopic reversibility.

Source
------
Loreau et al. compute their process (1), H+ + He(1s2) -> H(1s) + He+(1s),
with two fully quantal methods below 1 keV/u (which agree with each other and
with the measurements of Latypov & Shaporenko), and give a recommended fit,
their eq. (3) with the process (1) column of Table 2:

    log10(sigma / 1e-16 cm^2) = a0 + sum_{n=1}^{7} [a_n cos(n w log10 E)
                                                     + b_n sin(n w log10 E)],

E in eV/u.  They state its accuracy as "an uncertainty of less than 25%
between 100 and 400 eV/u ... and of less than 5% over the rest of the energy
range", the fit spanning 10 eV/u - 10 MeV/u; their quantal cross sections
(Figure 1) start at 15 eV/u.  The paper writes "log"; the base-10 reading is
the one that reproduces Figure 1 (2.3e-7 x 1e-16 cm^2 at 15 eV/u).

E in eV/u is taken as the laboratory energy per nucleon of the proton
projectile, the usual convention for p + He, so E_cm = E m_p m_He / (m_p +
m_He) in u.  The paper does not state the frame; reading E as E_cm instead
lowers the rate by 35-48% over 1e3-1e5 K (printed below as a sensitivity).

Reciprocity
-----------
Both channels are in their ground states.  With g_f = 2 x 2 = 4 (He+ 2S, H
2S) and g_r = 1 (He 1S, H+), the unpolarized forward cross section is

    sigma_f(E_f) = (g_r/g_f) (mu_r E_r)/(mu_f E_f) sigma_r(E_r),
    E_r = E_f + dE (center of mass),  dE = IP(He) - IP(H) = 10.989 eV.

The 1/4 is the ratio of statistical weights (the singlet approach
probability enters by construction).  k_f(T) = <sigma_f v> over a Maxwellian
in the forward relative energy.

Near threshold (2026-09-26)
--------------
Between the reverse threshold, x_th = dE (m_p + m_He)/(m_p m_He) = 13.66
eV/u, and 15 eV/u, where Loreau et al.'s quantal cross sections start
(their Figure 1), there are no data, and the Fourier fit does not vanish at
the threshold (2.0e-23 cm^2 there, against 2.3e-23 at 15 eV/u).  Carried
through the threshold, a finite sigma_r gives sigma_f ~ 1/E_f and a forward
rate that grows without bound as T falls (3.4e-16 at 1000 K); set to zero
it drops the open channel altogether (1.5e-21 at 1000 K).  Neither is the
physics.  For an exothermic ion-atom reaction at low collision energy the
product sigma_f v tends to a constant (capture on the ion-induced-dipole
potential, with a charge-transfer probability that does not vary with
energy); by reciprocity that is sigma_r proportional to
(E_r - dE)^(1/2), i.e. to (x - x_th)^(1/2), near threshold.  Below
15 eV/u sigma_r is therefore

    sigma_r(x) = sigma_r(15 eV/u) [(x - x_th)/(15 - x_th)]^(1/2),

joined to the fit at 15 eV/u.  The assumption is that the transfer
probability is constant over forward energies 0-1.08 eV, which no
calculation tests; the two limits above bracket it (at 1e4 K: 4.4e-17 with
zero, 1.27e-16 with the fit, 9.64e-17 with this form, beside Zygelman's
Figure 8 times 1/4, 9.3e-17).  EXHALE carries the same rate with the same
form since 2026-09-26 (its cooling_data/heii_h_charge_exchange_rate.py,
from this script; charge_exchange.f90, nonradiative_electron_capture_Hep_
from_H), reporting 9.644e-17 at 1e4 K, as here.

Validity (stated in the Fortran beside the table)
--------
 - 2e4 - 5e4 K: set by the quantal cross sections (the near-threshold form
   moves the rate by <= 13% here), within the frame reading (35-48%) and
   the fit accuracy.
 - below ~2e4 K: increasingly set by the near-threshold form above: at
   1e4 K the zero and fit limits are 0.46x and 1.32x the adopted value, at
   1000 K 1.5e-5x and 3.6x.  The adopted form gives a nearly constant
   rate, 9.6e-17 - 9.8e-17 over 1e3 - 2e4 K.  The table starts at 1000 K
   and is held below it.
 - above ~3e4 K: Loreau et al. disagree with Zygelman et al. (1989, PRA 40,
   2340), the calculation behind the Kingdon & Ferland (1996) fit used
   before, by up to 1e3 in the cross section above ~5 eV; they report the
   disagreement themselves (with Kimura et al. 1993, which used the Zygelman
   group's potential curves).  At 1e4 K the two agree (Zygelman's Figure 8,
   Maxwell-averaged and times the singlet 1/4, gives 9.3e-17 against
   9.6e-17 here; md/heii_h_charge_exchange/ANALYSIS_2026-09-25.md).  The
   Loreau value is adopted throughout because its two methods and the
   measurement agree; the high-temperature value is uncertain by that
   conflict.
The table ends at 1e6 K and is held above it.

Run:  python3 tools/fitting/make_heii_h_charge_exchange.py
"""
import math
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "src", "heii_h_charge_exchange_data.inc")

# Loreau et al. (2014) Table 2, process (1)
A = [-5.167, -5.200, -0.3265, -0.3786, -0.04399, -0.05318, -0.06957, 0.01072]
B = [0.0, 0.1131, 1.279, 0.4867, 0.1944, 0.1305, 0.03501, -0.04530]
W = 0.7800

U = 1.66053906660e-24          # g
EV = 1.602176634e-12           # erg
KB = 8.617333262e-5            # eV/K
ME = 5.485799090e-4            # u
MP = 1.007276467               # u
MHE_ATOM = 4.001506179 + 2*ME  # u (He nucleus + 2 electrons)
MU_R = MP*MHE_ATOM/(MP + MHE_ATOM)                  # H+ + He
MU_F = (MP + ME)*(MHE_ATOM - ME)/(MP + MHE_ATOM)    # He+ + H
DE = 24.587387 - 13.598434     # eV

LOGT = np.round(np.arange(3.0, 6.0 + 1e-9, 0.02), 10)


def sigma_rev(x):
    """Loreau process (1) [cm^2] at x eV/u (lab, per nucleon of the proton)."""
    L = np.log10(x)
    s = A[0] + sum(A[n]*np.cos(n*W*L) + B[n]*np.sin(n*W*L) for n in range(1, 8))
    return 10.0**s*1.0e-16


X_TH = DE*(MP + MHE_ATOM)/(MP*MHE_ATOM)     # reverse threshold [eV/u]
X_DATA = 15.0                               # first quantal energy [eV/u]


def k_forward(T, frame="lab", low="threshold", npts=200001):
    """<sigma_f v> [cm^3 s^-1]; below 15 eV/u, low = 'threshold' (the
    adopted (x - x_th)^(1/2) form), 'fit' (the fit as given), 'zero'."""
    kT = KB*T
    ef = np.exp(np.linspace(math.log(1e-7), math.log(60.0*kT + 60.0), npts))
    er = ef + DE
    x = er*(MP + MHE_ATOM)/(MP*MHE_ATOM) if frame == "lab" else er
    s = sigma_rev(x)
    if low == "zero":
        s = np.where(x < X_DATA, 0.0, s)
    elif low == "threshold":
        xth = X_TH if frame == "lab" else DE
        s15 = sigma_rev(np.array([X_DATA]))[0]
        s = np.where(x < X_DATA,
                     s15*np.sqrt(np.clip((x - xth)/(X_DATA - xth), 0.0, None)), s)
    sf = 0.25*(MU_R*er)/(MU_F*ef)*s
    f = ef*sf*np.exp(-ef/kT)*ef                     # E sigma e^{-E/kT} dE, dE = E dlnE
    integ = np.sum(0.5*(f[1:] + f[:-1])*np.diff(np.log(ef)))
    return math.sqrt(8.0/(math.pi*MU_F*U))*(kT*EV)**-1.5*integ*EV**2


def fortran_data(name, values, per_line=4):
    v = [f"{x:.10e}".replace("e", "E") + "_wp" for x in values]
    lines = [", ".join(v[i:i + per_line]) for i in range(0, len(v), per_line)]
    return f"  data {name} / &\n" + ", &\n".join("       " + ln for ln in lines) + " /\n"


def main():
    T = 10.0**LOGT
    k = np.array([k_forward(t) for t in T])
    # quadrature check: half the nodes
    k2 = np.array([k_forward(t, npts=100001) for t in T])
    quad = np.max(np.abs(k2/k - 1.0))
    with open(OUT, "w") as fh:
        fh.write(
            "! MoCHII: non-radiative He+ + H -> He + H+ [cm^3 s^-1], unpolarized,\n"
            "! from Loreau, Ryabchenko & Vaeck (2014, J. Phys. B 47, 135204) eq. (3)\n"
            "! and Table 2 (process 1, the reverse reaction) by microscopic\n"
            "! reversibility (statistical weights 1/4); below 15 eV/u, where the\n"
            "! quantal data stop, sigma_r ~ (x - x_th)^(1/2) joined at 15 eV/u\n"
            "! (sigma_f v constant at low energy).  GENERATED by\n"
            "! tools/fitting/make_heii_h_charge_exchange.py -- do not edit; rerun\n"
            "! the script.  Included by ion_balance_mod.\n"
            f"! Quadrature: {quad:.1e} relative between 200001 and 100001 nodes.\n"
            "!   heii_h_cx_logT   log10 T [K]\n"
            "!   heii_h_cx_logk   log10 k [cm^3 s^-1]\n")
        fh.write(f"  integer, parameter :: HEII_H_CX_NT = {len(LOGT)}\n")
        fh.write("  real(kind=wp), save :: heii_h_cx_logT(HEII_H_CX_NT)\n")
        fh.write("  real(kind=wp), save :: heii_h_cx_logk(HEII_H_CX_NT)\n")
        fh.write(fortran_data("heii_h_cx_logT", LOGT))
        fh.write(fortran_data("heii_h_cx_logk", np.log10(k)))
    print(f"wrote {os.path.normpath(OUT)}  (quadrature {quad:.1e})")

    print("\n   T [K]    k (adopted)  [sigma=0, fit as given below 15 eV/u]   E as E_cm   KF96 (before)")
    for t in (1e3, 3e3, 6e3, 1e4, 2e4, 3e4, 5e4, 1e5, 1e6):
        t4 = min(t, 1e5)/1e4
        kf96 = 7.47e-15*t4**2.06*(1 + 9.93*math.exp(-3.89*t4))
        print(f"{t:9.0f}  {k_forward(t):.3e}  [{k_forward(t, low='zero'):.3e}, "
              f"{k_forward(t, low='fit'):.3e}]  {k_forward(t, frame='cm'):.3e}  {kf96:.3e}")
    # check: the reverse rate averaged directly over sigma_r must equal
    # 4 (mu_f/mu_r)^1.5 exp(-dE/kT) k_f
    for t in (1e4, 1e5):
        kT = KB*t
        er = DE + np.exp(np.linspace(math.log(1e-7), math.log(60*kT + 60), 200001))
        s = sigma_rev(er*(MP + MHE_ATOM)/(MP*MHE_ATOM))
        f = er*s*np.exp(-er/kT)*(er - DE)
        integ = np.sum(0.5*(f[1:] + f[:-1])*np.diff(np.log(er - DE)))
        kr = math.sqrt(8/(math.pi*MU_R*U))*(kT*EV)**-1.5*integ*EV**2
        print(f"reciprocity check at {t:.0e} K: {kr/(4*(MU_F/MU_R)**1.5*math.exp(-DE/kT)*k_forward(t, low='fit')):.6f}"
              " (the fit as given on both sides)")


if __name__ == "__main__":
    main()
