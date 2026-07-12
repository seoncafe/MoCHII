#!/usr/bin/env python3
"""Storey & Hummer (1995, MNRAS 272, 41) H I case-B line emissivities.

Parses the extended-format data file e1bx.d (the MOCASSIN-distributed
repackaging of the SH95 CDS data; identical reader as MOCASSIN
output_mod::hdatx and the reference reader in
~/RT_Codes/Storey_Hummer/read_data.f90):

  header: ntemp ndens
  blocks (ia=1..ntemp, ib=1..ndens):
    dens  temp  ntop ndum nlu nll
    e(j), j = 1..(2*ntop-nlu-nll)(nlu-nll+1)/2
  index of transition (nu -> nl):
    k = (2*ntop - nll - nl + 1)(nl - nll)/2 + ntop - nu + 1

Values are 4 pi j / (n_e n_p) [erg cm^3 s^-1].  Writes
data/atomic/sh95_hi_caseB.txt with the principal lines on the full
(T, n_e) grid for bilinear log interpolation in Fortran.
"""
import os
import numpy as np

SRC = os.path.expanduser("~/RT_Codes/Storey_Hummer/e1bx.d")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "..", "data", "atomic", "sh95_hi_caseB.txt")

# (n_upper, n_lower, label, vacuum wavelength [A])
LINES = [(3, 2, "Halpha", 6564.6), (4, 2, "Hbeta", 4862.7),
         (5, 2, "Hgamma", 4341.7), (6, 2, "Hdelta", 4102.9),
         (4, 3, "Palpha", 18756.1), (5, 3, "Pbeta", 12821.6),
         (7, 4, "Brgamma", 21661.2)]


def main():
    toks = open(SRC).read().split()
    p = 0
    ntemp, ndens = int(toks[p]), int(toks[p+1]); p += 2
    temps = np.zeros(ntemp); denss = np.zeros(ndens)
    E = {}
    for ia in range(ntemp):
        for ib in range(ndens):
            dens = float(toks[p]); flag = toks[p+1]
            temp = float(toks[p+2]); case = toks[p+3]
            ntop, ndum, nlu, nll = (int(t) for t in toks[p+4:p+8])
            p += 8
            ne_vals = (2*ntop - nlu - nll)*(nlu - nll + 1)//2
            vals = np.array(toks[p:p+ne_vals], dtype=float)
            p += ne_vals
            temps[ia] = temp;  denss[ib] = dens
            for (nu, nl, lab, wl) in LINES:
                k = (2*ntop - nll - nl + 1)*(nl - nll)//2 + ntop - nu + 1
                E.setdefault(lab, np.zeros((ntemp, ndens)))[ia, ib] = vals[k-1]

    with open(OUT, "w") as fh:
        fh.write("# Storey & Hummer (1995, MNRAS 272, 41) H I case-B line "
                 "emissivities 4 pi j / (n_e n_p) [erg cm^3 s^-1]\n")
        fh.write("# source file: Storey_Hummer/e1bx.d (MOCASSIN-format "
                 "repackaging of the CDS tables), parsed by "
                 "tools/fitting/make_sh95_lines.py (2026-07-12)\n")
        fh.write("# grid: NT temperatures [K], ND densities [cm^-3]; then "
                 "one block for each line: label nu nl lambda[A], NT x ND "
                 "values (T rows, n_e columns)\n")
        fh.write(f"GRID {ntemp} {ndens}\n")
        fh.write("T " + " ".join(f"{t:.4e}" for t in temps) + "\n")
        fh.write("NE " + " ".join(f"{d:.4e}" for d in denss) + "\n")
        for (nu, nl, lab, wl) in LINES:
            fh.write(f"LINE {lab} {nu} {nl} {wl:.1f}\n")
            for ia in range(ntemp):
                fh.write(" ".join(f"{E[lab][ia, ib]:.4e}"
                                  for ib in range(ndens)) + "\n")
    print(f"wrote {os.path.relpath(OUT)} ({ntemp} T x {ndens} ne)")
    # sanity: Hbeta at T=1e4, ne=1e2 should be ~1.24e-25 erg cm^3/s
    ia = np.argmin(abs(temps - 1e4)); ib = np.argmin(abs(denss - 1e2))
    print(f"check: 4pi j(Hbeta)/(ne np) at T={temps[ia]:.0f}, "
          f"ne={denss[ib]:.0f} = {E['Hbeta'][ia, ib]:.4e} (expect ~1.24e-25)")
    print(f"check: Halpha/Hbeta = {E['Halpha'][ia, ib]/E['Hbeta'][ia, ib]:.4f}"
          f" (expect ~2.87)")


if __name__ == "__main__":
    main()
