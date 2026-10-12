#!/usr/bin/env python3
"""Reference values for tests/ne3_ground_term_photoionization, computed from
the published source alone, independently of the generator
(tools/fitting/ground_term_photoionization_table.py) and of the Fortran.

Source: Nahar (2019, New Astron. 67, 97), NORAD ne3.px.gd.txt, stored
unchanged as data/atomic/norad/ne3.px.gd.txt: the total photoionization
cross sections (with radiation damping) of Ne III 2s2 2p4 3P2, 3P1, 3P0.

(1) sigma(E) of the ground term as MoCHII defines it: sum_J (2J+1) sigma_J /
    9 on the photoelectron-energy scale eps = E - ETH, each level
    interpolated linearly on its OWN grid (the generator uses the 3P2 grid
    for all three), at 3000 energies drawn uniformly in eps below E_join
    (fixed seed) plus 1000 inside the first 2 eV, where the resonances are.
(2) The ground-term recombination rate the Milne relation gives on the
    LEVEL cross sections, each at its own threshold IP_J = ETH - E_J
    (E_J: the fine-structure energies of CHIANTI 11.0.2 ne_3.elvlc,
    observed), summed over J with weights (2J+1)/(2 g(Ne IV 4S) = 8):
        alpha_1 = sum_J (2J+1)/8 (2/sqrt(pi)) (kT)^(-3/2) sqrt(2/m_e)/(m_e c^2)
                  Int E^2 sigma_J(E) exp[-(E - IP_J)/kT] dE ,
    integrated exactly over each level's piecewise-linear table (Simpson on
    each interval) up to its photoelectron energy at E_join.  MoCHII places
    the three thresholds at ETH; the reference does not, so the comparison
    also measures that choice.  The tail above E_join (VFKY96 in MoCHII) is
    left out: exp(-(E_join - ETH)/kT) = exp(-26.5) at 5e4 K.
Nahar has not published Ne III recombination rate coefficients: NORAD
carries no ne3.rrc.txt (HTTP 404, 2026-10-11) and a NASA ADS search of her
neon papers finds none, so (2) -- the Milne relation on her published cross
sections, the low-n part of her own unified method -- is the reference.

Run from the repository root:  python3 tests/ne3_ground_term_photoionization/make_reference.py
"""
import math
import os
import random

RY = 13.605693122994                   # eV
KB = 1.380649e-16                      # erg/K
ME = 9.1093837015e-28                  # g
C = 2.99792458e10                      # cm/s
EV = 1.602176634e-12                   # erg/eV
ETH = 6.342331081121e+01               # NIST IP of Ne III 3P2 (element_ne.txt ETH, transition 3)
GPLUS = 4                              # Ne IV 2s2 2p3 4S3/2
TEMPS = [3.0e3, 5.0e3, 1.0e4, 2.0e4, 5.0e4]
SRC = os.path.join("data", "atomic", "norad", "ne3.px.gd.txt")
ELVLC = os.path.expanduser("~/RT_Codes/CHIANTI/dbase/ne/ne_3/ne_3.elvlc")
OUT = os.path.join("tests", "ne3_ground_term_photoionization", "reference.txt")


def read_levels(path):
    lines = open(path).read().split("\n")
    out, j = {}, 0
    while j < len(lines) - 1:
        t, u = lines[j].split(), lines[j + 1].split()
        if len(t) == 4 and t[0] == "0" and t[2] == "0" and t[3] == "1" and len(u) == 3:
            n = int(u[1])
            rows = [r.split() for r in lines[j + 3:j + 3 + n]]
            out[int(t[1])] = ([float(r[1]) * RY for r in rows], [float(r[3]) for r in rows])
            j += 3 + n
        else:
            j += 1
    return out


def interp(x, xs, ys):
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        m = (lo + hi) // 2
        if xs[m] <= x:
            lo = m
        else:
            hi = m
    return ys[lo] + (x - xs[lo]) / (xs[hi] - xs[lo]) * (ys[hi] - ys[lo])


def fine_structure_cm():
    """{2J: E [cm^-1]} of 2s2 2p4 3P from CHIANTI ne_3.elvlc (observed)"""
    e = {}
    for line in open(ELVLC):
        t = line.split()
        if len(t) >= 7 and t[0] in ("1", "2", "3"):
            e[int(round(2 * float(t[4])))] = float(t[5])
    return e


def e_join():
    for line in open(os.path.join("data", "atomic", "element_ne.txt")):
        if line.startswith("# PHTABLE: Nahar"):
            return float(line.split("E_join = ")[1].split()[0])
    raise SystemExit("element_ne.txt has no Ne III PHTABLE comment")


def main():
    lev = read_levels(SRC)
    g = {4: 5, 2: 3, 0: 1}
    ej = e_join()
    eps_max = ej - ETH
    rnd = random.Random(20261011)
    pts = sorted([rnd.uniform(0.0, eps_max) for _ in range(3000)] +
                 [rnd.uniform(0.0, 2.0) for _ in range(1000)])
    sig = [sum(g[k] * interp(e, *lev[k]) for k in g) / 9.0 for e in pts]
    efs = fine_structure_cm()
    pref = lambda T: (2 / math.sqrt(math.pi)) * (KB * T) ** -1.5 * math.sqrt(2 / ME) / (ME * C * C)
    alpha = []
    for T in TEMPS:
        kT = KB * T / EV
        tot = 0.0
        for k, w in g.items():
            ip = ETH - efs[k] / 8065.543937        # cm^-1 -> eV
            xs, ys = lev[k]
            acc = 0.0
            for i in range(len(xs) - 1):
                a, b = xs[i], min(xs[i + 1], eps_max)
                if b <= a:
                    break
                sa, sb = ys[i], ys[i] + (b - a) / (xs[i + 1] - xs[i]) * (ys[i + 1] - ys[i])
                f = lambda e, s: (ip + e) ** 2 * s * math.exp(-e / kT)
                m = 0.5 * (a + b)
                acc += (b - a) / 6 * (f(a, sa) + 4 * f(m, 0.5 * (sa + sb)) + f(b, sb))
            tot += w / (2.0 * GPLUS) * acc
        # [eV^3 Mb] -> erg^3 cm^2
        alpha.append(pref(T) * tot * EV ** 3 * 1e-18)
    with open(OUT, "w") as fh:
        fh.write("# reference for check_ne3_ground_term_photoionization.f90; written by "
                 "make_reference.py from data/atomic/norad/ne3.px.gd.txt (see its docstring)\n")
        fh.write(f"# E_join {ej:.10e}\n")
        fh.write(f"ALPHA1 {len(TEMPS)}\n")
        for T, a in zip(TEMPS, alpha):
            fh.write(f"{T:.6e} {a:.10e}\n")
        fh.write(f"SIGMA {len(pts)}\n")
        for e, s in zip(pts, sig):
            fh.write(f"{ETH + e:.14e} {s * 1e-18:.12e}\n")
    for T, a in zip(TEMPS, alpha):
        print(f"T = {T:8.0f} K  alpha_1(levels at their own thresholds) = {a:.5e} cm^3/s")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
