#!/usr/bin/env python3
"""Write data/atomic/h_2s_2p_mixing.txt: rate coefficients for the collisional
transfer H(2s) -> H(2p) by electrons, protons, He+ and He2+, which quench the
H I two-photon continuum,

    f_2q = A_2s / (A_2s + n_e q_e + n_p q_p + n(He+) q_He+ + n(He2+) q_He2+) .

Physics: Pengelly & Seaton (1964, MNRAS 127, 165; read from the ADS scan).
2s and 2p are nearly degenerate, so a charged particle passing at impact
parameter R drives 2s -> 2p with probability P(R) = 4 Z^2 e^2 S/(3 hbar^2 R^2
v^2) (their eq. 34); P = 1/2 inside R_1 (eqs. 35-37), and the integral is cut
at the smallest of R_c = 1.12 hbar v/dE, 0.72 v tau and the Debye radius
(p. 170), giving Q = pi R_1^2 [1/2 + ln(R_c/R_1)] (eq. 36) and q = <Q v>
(eq. 42).  Z is the charge of the colliding particle, and v the relative
velocity for the reduced mass of the pair (eq. 43).  For 2s the line strength
is split between the two final levels, 2p_1/2 (weight 1/3, dE = the Lamb
shift, 0.0353 cm^-1) and 2p_3/2 (weight 2/3, dE = 0.330 cm^-1).  For those
separations the dE cut-off is the smallest of the three wherever
n_e < ~1e12 cm^-3, so the table is the low-density limit.  0.72 v tau, with
tau = 1/A_2s, is never the smallest.  Where R_c < R_1 the cross section is
pi R_c^2 / 2 (P = 1/2 out to R_c), which is the limit of eq. 36 that the
closed form (their eq. 43) does not carry: it matters for He2+ on the
2p_3/2 separation below ~1e4 K.

Normalization: Osterbrock & Ferland (2006) Table 4.10, which the book states
"comes from" Pengelly & Seaton (1964), prints the proton and electron rates
for both final levels at 1e4 and 2e4 K; they sum to the 5.31e-4 / 4.71e-4
cm^3 s^-1 that MOCASSIN (and MoCHII before this table) used as one lumped
coefficient.  The Maxwellian quadrature below reproduces the level split and
the temperature dependence of that table but sits above it by 13% (1e4 K) and
9.5% (2e4 K) for protons and 4.7% and 2.8% for electrons; the printed numbers are taken as the published
values, and the quadrature supplies (a) the temperature dependence outside
1e4-2e4 K and (b) He+ and He2+, which Table 4.10 does not give.  Each
collider is written as quadrature x c(T), with c = printed/quadrature for
electrons and for protons, interpolated in log T between 1e4 and 2e4 K and
held outside; He+ and He2+ take the proton c, i.e. the table keeps the
quadrature's ratios q_He+/q_p and q_He2+/q_p (1.07 and 3.21 at 1e4 K), in which
a normalization common to all heavy colliders cancels.

A_2s = 8.2246 s^-1: Drake (1986, PRA 34, 2871) eq. 27,
8.22938 Z^6 Z_r^4 (1 - m/M) [1 + 3.9448 (aZ)^2 - 2.040 (aZ)^4]/[1 + 4.6019 (aZ)^2]
with Z = Z_r = 1 and M the proton mass.
"""
import os

import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "..", "..",
                   "data", "atomic", "h_2s_2p_mixing.txt")

KB = 1.380649e-16          # erg/K
ME = 9.1093837015e-28      # g
MP = 1.67262192369e-24     # proton
MH = 1.6735575e-24         # H atom
MHE = 6.6464731e-24        # He atom
HBAR = 1.054571817e-27
C = 2.99792458e10
ALPHA = 7.2973525693e-3
A0 = 5.29177210903e-9
E_ESU = 4.80320471e-10
H = 2.0*np.pi*HBAR

# final levels of 2s_1/2: statistical-weight share of the line strength and
# the energy separation from 2s_1/2 [cm^-1]
LEVELS = {"2p1/2": (1.0/3.0, 0.0353), "2p3/2": (2.0/3.0, 0.330)}

COLLIDERS = {              # reduced mass with H, charge
    "e":    (ME*MH/(ME + MH), 1),
    "p":    (MP*MH/(MP + MH), 1),
    "He+":  ((MHE - ME)*MH/(MHE - ME + MH), 1),
    "He2+": ((MHE - 2*ME)*MH/(MHE - 2*ME + MH), 2),
}

# Osterbrock & Ferland (2006) Table 4.10 [cm^3 s^-1], summed over 2p_1/2 and
# 2p_3/2: protons 2.51 + 2.23 (1e4 K), 2.08 + 2.19 (2e4 K); electrons
# 0.22 + 0.35, 0.17 + 0.27 (all x 1e-4).
OF_T = np.array([1.0e4, 2.0e4])
OF_Q = {"p": np.array([4.74e-4, 4.27e-4]), "e": np.array([0.57e-4, 0.44e-4])}

A_2S = 8.22938*(1.0 - ME/MP)*(1.0 + 3.9448*ALPHA**2 - 2.040*ALPHA**4) \
       / (1.0 + 4.6019*ALPHA**2)


def q_quadrature(T, collider, ne=1.0):
    """<Q v> of Pengelly & Seaton's cross section over a Maxwellian of the
    relative velocity, summed over the two final levels."""
    mu, Z = COLLIDERS[collider]
    x = np.linspace(-14.0, 5.0, 20000)
    v = np.exp(x)*np.sqrt(2.0*KB*T/mu)
    f = 4.0*np.pi*v**2*(mu/(2.0*np.pi*KB*T))**1.5*np.exp(-mu*v**2/(2.0*KB*T))
    r_debye = np.sqrt(KB*T/(4.0*np.pi*ne*E_ESU**2))
    total = 0.0
    for w, de_cm in LEVELS.values():
        # R_1^2 = 8 Z^2 e^2 S/(3 hbar^2 v^2), S = w 27 e^2 a0^2 (eqs. 37, 40)
        r1sq = w*72.0*Z**2*(ALPHA*C*A0/v)**2
        rc = np.minimum(1.12*HBAR*v/(de_cm*H*C), r_debye)
        rc = np.minimum(rc, 0.72*v/A_2S)
        r1 = np.sqrt(r1sq)
        Q = np.where(rc > r1,
                     np.pi*r1sq*(0.5 + np.log(np.maximum(rc/r1, 1.0))),
                     0.5*np.pi*rc**2)
        total += np.trapz(Q*v*f*v, x)          # dv = v dx
    return total


def normalization(T, collider):
    key = "e" if collider == "e" else "p"
    c = OF_Q[key]/np.array([q_quadrature(t, key) for t in OF_T])
    lt = np.clip(np.log10(T), np.log10(OF_T[0]), np.log10(OF_T[1]))
    return np.interp(lt, np.log10(OF_T), c)


def main():
    # log grid plus the two temperatures Table 4.10 prints, exactly
    temps = np.unique(np.concatenate([np.logspace(3.0, 5.0, 41), OF_T]))
    rows = []
    for T in temps:
        rows.append([T] + [q_quadrature(T, k)*normalization(T, k)
                           for k in COLLIDERS])
    rows = np.array(rows)

    # checks: Table 4.10 reproduced at its two temperatures
    for t, qp, qe in zip(OF_T, OF_Q["p"], OF_Q["e"]):
        i = np.flatnonzero(rows[:, 0] == t)[0]
        assert abs(rows[i, 2]/qp - 1) < 1e-12 and abs(rows[i, 1]/qe - 1) < 1e-12
    # the density at which the Debye radius would take over from dE
    ne_debye = None
    for ne in np.logspace(8, 16, 81):
        if abs(q_quadrature(1.0e4, "p", ne)/q_quadrature(1.0e4, "p") - 1) > 1e-3:
            ne_debye = ne
            break

    hdr = [
        "# MoCHII: H(2s) -> H(2p) collisional rate coefficients [cm^3 s^-1] and",
        "# the 2s two-photon decay rate, for the quenching of the H I two-photon",
        "# continuum: f = A_2S/(A_2S + n_e q_e + n_p q_p + n(He+) q_He+ + n(He2+) q_He2+)",
        "# physics: Pengelly & Seaton (1964, MNRAS 127, 165) eqs. 34-43, integrated",
        "#   over the Maxwellian of the relative velocity, both final levels",
        "#   (2p1/2: 1/3, 0.0353 cm^-1; 2p3/2: 2/3, 0.330 cm^-1), dE cut-off",
        "# normalization: Osterbrock & Ferland (2006) Table 4.10 (from P&S 1964)",
        "#   reproduced exactly at 1e4 and 2e4 K for e and p; He+ and He2+ keep the",
        "#   quadrature's ratio to protons",
        "# low-density limit: the Debye cut-off changes q_p by 0.1%% at n_e = %.1e"
        % ne_debye,
        "# A_2S: Drake (1986, PRA 34, 2871) eq. 27 for Z = 1, M = m_p",
        "# written by tools/fitting/make_h2s_collisional_mixing.py",
        "# columns: T [K], q_e, q_p, q_He+, q_He2+; log-log interpolation,",
        "#   held at the end values outside 1e3-1e5 K",
        "A_2S  %.6f" % A_2S,
        "TABLE %d" % len(rows),
    ]
    with open(OUT, "w") as fh:
        fh.write("\n".join(hdr) + "\n")
        for r in rows:
            fh.write("  %.4e  %.5e  %.5e  %.5e  %.5e\n" % tuple(r))
    i4 = np.flatnonzero(rows[:, 0] == 1.0e4)[0]
    print("A_2S = %.6f s^-1" % A_2S)
    print("1e4 K: q_e %.3e q_p %.3e q_He+ %.3e q_He2+ %.3e" % tuple(rows[i4, 1:]))
    print("Debye cut-off matters from n_e ~ %.1e cm^-3" % ne_debye)
    print("wrote", os.path.normpath(OUT))


if __name__ == "__main__":
    main()
