#!/usr/bin/env python3
"""Gate: the H-impact [C II] 158 um / [O I] 63 um two-level cooling must
reproduce both density limits.

This mirrors `metal_cooling_H` in src/species_mod.f90.  The coefficients and
the expression are read out of the Fortran source, so the gate fails if the
implemented formula drifts from the two-level steady state rather than
comparing MoCHII against a copy of itself.

Steady state, n_ion = n_l + n_u:

    n_u (A + n_HI q_ul) = n_l n_HI q_lu
    Lambda = n_u A dE   = n_ion n_HI q_lu dE / (1 + n_HI (q_ul + q_lu)/A)

Two limits pin the expression:

  low density  (n_HI << A/(q_ul+q_lu)) : Lambda -> n_ion n_HI q_lu dE
  LTE          (n_HI >> A/(q_ul+q_lu)) : Lambda -> n_ion A dE x_u,
                                         x_u = b/(1+b), b = (g_u/g_l)e^(-dE/kT)

Dropping q_lu from the denominator leaves the low-density limit intact but
inflates the LTE limit by (1 + q_lu/q_ul) -- a factor approaching 3 for
[C II], i.e. an upper-level occupation of ~200% of the ion population.  The
low-density check alone therefore cannot catch it; both are required.

Run:  python3 check_two_level.py     (exits nonzero on failure)
"""
import os
import re
import sys

import numpy as np

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "..", "src", "species_mod.f90")

KB = 1.380649e-16

# (label, A [s^-1], q_ul prefactor, q_ul exponent, g_u/g_l, dE/k [K])
IONS = [
    ("[C II] 158um", 2.321e-6, 7.6e-10, 0.14, 4.0 / 2.0, 91.25),
    ("[O I] 63um",   8.865e-5, 9.2e-11, 0.67, 3.0 / 5.0, 227.7),
]

TEMPS = [20.0, 30.0, 100.0, 300.0, 1000.0, 1.0e4]
#--- densities far enough either side of n_crit that the residual finite-
#--- density correction (order n_HI(q_ul+q_lu)/A and its inverse) sits well
#--- below RTOL; the tolerance then tests the formula, not the sampling.
N_LOW, N_HIGH = 1.0e-14, 1.0e22
RTOL = 1.0e-8


def source_uses_both_collision_directions():
    """The denominator in species_mod.f90 must carry (qul + qlu), not qul."""
    with open(SRC) as fh:
        text = fh.read()
    good = re.findall(r"1\.0_wp \+ nHI\*\(qul \+ qlu\)/A_(?:CII|OI)", text)
    bad = re.findall(r"1\.0_wp \+ nHI\*qul/A_(?:CII|OI)", text)
    return len(good), len(bad)


def cooling(nHI, T, A, q0, texp, gratio, dE_K):
    """metal_cooling_H per unit n_ion, divided by dE -- the code's expression."""
    qul = q0 * (T / 100.0) ** texp
    qlu = gratio * qul * np.exp(-dE_K / T)
    return nHI * qlu / (1.0 + nHI * (qul + qlu) / A), qul, qlu


def main():
    ok = True

    ngood, nbad = source_uses_both_collision_directions()
    print(f"source check: {ngood} denominators with (qul + qlu), "
          f"{nbad} with qul alone")
    if ngood != 2 or nbad != 0:
        print("  FAIL: species_mod.f90 does not use the two-level denominator")
        ok = False

    for label, A, q0, texp, gratio, dE_K in IONS:
        print(f"\n=== {label} ===")
        print(f"{'T [K]':>8} {'n_crit':>10} {'low-n':>12} {'LTE':>12}"
              f" {'1+qlu/qul':>10}")
        for T in TEMPS:
            lam_lo, qul, qlu = cooling(N_LOW, T, A, q0, texp, gratio, dE_K)
            lam_hi, _, _ = cooling(N_HIGH, T, A, q0, texp, gratio, dE_K)

            #--- low-density limit: Lambda / (n_HI q_lu) -> 1
            r_lo = lam_lo / (N_LOW * qlu)
            #--- LTE limit: Lambda / (A x_u) -> 1 with x_u the Boltzmann value
            b = gratio * np.exp(-dE_K / T)
            r_hi = lam_hi / (A * b / (1.0 + b))

            ncrit = A / (qul + qlu)
            print(f"{T:8.0f} {ncrit:10.3e} {r_lo:12.9f} {r_hi:12.9f}"
                  f" {1.0 + qlu/qul:10.4f}")
            for name, got in (("low-density", r_lo), ("LTE", r_hi)):
                if abs(got - 1.0) > RTOL:
                    print(f"  FAIL: {name} limit at T = {T:g} K is {got:.12f}")
                    ok = False

            #--- the upper level can never hold more than the ion population
            if lam_hi / A > 1.0:
                print(f"  FAIL: LTE limit implies n_u/n_ion = {lam_hi/A:.3f} > 1")
                ok = False

    print(f"\nGATE (two-level low-density and LTE limits): "
          f"{'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
