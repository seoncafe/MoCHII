#!/usr/bin/env python3
"""He+ + H -> He + H+ (non-radiative, ground states) from Loreau, Ryabchenko
& Vaeck (2014, J. Phys. B 47, 135204), through microscopic reversibility.

Loreau et al. compute the reverse reaction, their process (1),
    H+ + He(1s2 1S) -> H(1s) + He+(1s),
and give a fit (their eq. 3, Table 2)
    log10(sigma / 1e-16 cm^2) = a0 + sum_n [a_n cos(n w log10 E) + b_n sin(n w log10 E)],
with E in eV/u, stated uncertainty < 25%, quantal data from 15 eV/u.

Reciprocity (spin-averaged, both channels in their ground states):
    g_f k_f^2 sigma_f(E_f) = g_r k_r^2 sigma_r(E_r),   E_r = E_f + dE (CM),
    g_f = 2 x 2 = 4 (He+ 2S, H 2S),  g_r = 1 x 1 = 1 (He 1S, H+),
so the unpolarized forward cross section is
    sigma_f = (1/4) (mu_r E_r)/(mu_f E_f) sigma_r(E_r).
The 1/4 is the ratio of statistical weights, not an assumption about a fit.

Energy convention: E [eV/u] is taken as the laboratory energy per nucleon
of the proton projectile (standard for p + He), E_cm = E m_p m_He/(m_p+m_He)
in u; the alternative reading E_cm = E [eV] is evaluated as a sensitivity.

Near threshold (E below the 15 eV/u where the quantal data start) three
variants bracket the rate: the fit as given, sigma = 0, and sigma held at its
15 eV/u value.  The +-25% fit uncertainty scales every value.

Also printed: the reverse rate He + H+ -> He+ + H averaged directly over
sigma_r (a check on the reciprocity arithmetic, which must reproduce
k_r = 4 (mu_f/mu_r)^(3/2) exp(-dE/kT) k_f), MoCHII's current rates
(KF96 non-radiative, Glover & Jappsen R26 radiative), and alpha_B(He I).

Standard library + numpy.  Run:  python3 loreau2014_rate.py
"""
import math
import numpy as np

A = [-5.167, -5.200, -0.3265, -0.3786, -0.04399, -0.05318, -0.06957, 0.01072]
B = [0.0, 0.1131, 1.279, 0.4867, 0.1944, 0.1305, 0.03501, -0.04530]
W = 0.7800

U = 1.66053906660e-24          # g
EV = 1.602176634e-12           # erg
KB = 8.617333262e-5            # eV/K
MP, MHE, ME = 1.007276467, 4.001506179 + 2*5.485799e-4, 5.485799e-4   # u (He atom)
M_HE_ION = MHE - ME
M_H = MP + ME
MU_R = MP*MHE/(MP + MHE)               # H+ + He
MU_F = M_H*M_HE_ION/(M_H + M_HE_ION)   # He+ + H
DE = 24.587387 - 13.598434             # eV, IP(He) - IP(H)


def sigma_rev(x):
    """Loreau process (1) [cm^2] at x eV/u."""
    L = np.log10(x)
    s = A[0] + sum(A[n]*np.cos(n*W*L) + B[n]*np.sin(n*W*L) for n in range(1, 8))
    return 10.0**s*1.0e-16


def x_of_ecm(ecm, convention):
    return ecm*(MP + MHE)/(MP*MHE) if convention == 'lab' else ecm


def maxwell(T, e, sig, mu):
    """<sigma v> over relative energy e [eV] (log grid) [cm^3 s^-1]."""
    kT = KB*T
    f = e*sig*np.exp(-e/kT)*e                      # E sigma e^{-E/kT} dE, dE = E dlnE
    integ = np.sum(0.5*(f[1:] + f[:-1])*np.diff(np.log(e)))
    return math.sqrt(8.0/(math.pi*mu*U))*(kT*EV)**-1.5*integ*EV**2


def k_forward(T, convention='lab', low='fit'):
    kT = KB*T
    ef = np.exp(np.linspace(math.log(1e-5), math.log(80*kT + 50), 60000))
    er = ef + DE
    x = x_of_ecm(er, convention)
    s = sigma_rev(x)
    x15 = 15.0
    if low == 'zero':
        s = np.where(x < x15, 0.0, s)
    elif low == 'flat':
        s = np.where(x < x15, sigma_rev(np.array([x15]))[0], s)
    sf = 0.25*(MU_R*er)/(MU_F*ef)*s
    return maxwell(T, ef, sf, MU_F)


def k_reverse_direct(T, convention='lab'):
    kT = KB*T
    er = DE + np.exp(np.linspace(math.log(1e-5), math.log(80*kT + 50), 60000))
    s = sigma_rev(x_of_ecm(er, convention))
    # Maxwellian over the reverse channel's relative energy, open above DE
    f = er*s*np.exp(-er/kT)*(er - DE)
    integ = np.sum(0.5*(f[1:] + f[:-1])*np.diff(np.log(er - DE)))
    return math.sqrt(8.0/(math.pi*MU_R*U))*(kT*EV)**-1.5*integ*EV**2


def kf96(T):
    t = min(T, 1e5)/1e4
    return 7.47e-15*t**2.06*(1 + 9.93*math.exp(-3.89*t))


def radiative(T):
    return 1.25e-15*(T/300.0)**0.25


def alphaB_HeI(T):
    """MoCHII 'badnell_milne' alpha_B(He II -> He I) at 1e4 K scale, for
    comparison only: Badnell alpha_A minus the Milne alpha_1 (recomb_mod)."""
    def rr(T, A_, B_, T0, T1, C, T2):
        bp = B_ + C*math.exp(-T2/T); tt = math.sqrt(T/T0)
        return A_/(tt*(1 + tt)**(1 - bp)*(1 + math.sqrt(T/T1))**(1 + bp))
    aA = rr(T, 5.235e-11, 0.6988, 7.301, 4.475e6, 0.0829, 1.682e5) + \
        T**-1.5*(1.417e-3*math.exp(-4.633e5/T) + 2.235e-4*math.exp(-5.532e5/T)
                 - 2.185e-5*math.exp(-8.887e5/T))
    te = T*KB
    a1 = 3.26862e-3*1e-10*te**(-6.97361e-1 - 3.14907e-4*math.log(te)) * \
        (1 + 2.58417e1*te**-9.90134e-1)/(1 + 5.76341e1*te**-1.18534)
    return aA - a1


def main():
    temps = [1e3, 3e3, 6e3, 1e4, 1.5e4, 2e4, 3e4, 5e4, 1e5]
    print('# forward He+ + H -> He + H+, non-radiative, unpolarized [cm^3 s^-1]')
    print('# T_K  k_fit  k_zero<15eV/u  k_flat<15eV/u  k_fit_CMconv  '
          'kr_direct/kr_reciprocity  KF96  radiative  alphaB_HeI  '
          'k_fit/radiative  KF96/k_fit')
    for T in temps:
        kf = k_forward(T); kz = k_forward(T, low='zero'); kfl = k_forward(T, low='flat')
        kc = k_forward(T, convention='cm')
        krd = k_reverse_direct(T)
        krr = 4*(MU_F/MU_R)**1.5*math.exp(-DE/(KB*T))*kf
        chk = krd/krr if krr > 0 else float('nan')
        print('%7.0f %.3e %.3e %.3e %.3e %.6f %.3e %.3e %.3e %.4f %.1f' % (
            T, kf, kz, kfl, kc, chk, kf96(T), radiative(T), alphaB_HeI(T),
            kf/radiative(T), kf96(T)/kf))


if __name__ == '__main__':
    main()
