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


SRC2 = os.path.expanduser("~/RT_Codes/Storey_Hummer/data/e2b.d.gz")
OUT2 = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "data", "atomic", "sh95_heii_caseB.txt")

# He II lines (n_upper, n_lower, label, vacuum wavelength [A])
LINES2 = [(3, 2, "HeII1640", 1640.4), (4, 3, "HeII4686", 4687.0),
          (5, 3, "HeII3204", 3204.0), (5, 4, "HeII10126", 10126.4),
          (7, 4, "HeII5413", 5413.0), (9, 4, "HeII4543", 4542.9)]

# case-A source files (original CDS format, gzipped)
SRC_HI_A = os.path.expanduser("~/RT_Codes/Storey_Hummer/data/e1a.d.gz")
OUT_HI_A = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "data", "atomic", "sh95_hi_caseA.txt")
SRC_HEII_A = os.path.expanduser("~/RT_Codes/Storey_Hummer/data/e2a.d.gz")
OUT_HEII_A = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "..", "data", "atomic", "sh95_heii_caseA.txt")


def read_cds(src, out, lines, want_case, Z, ion_label, norm, provenance,
             checks="heii"):
    """Read the ORIGINAL CDS format (e1a.d / e2a.d style) and write the
    standard MoCHII (T, n_e) grid file.

    Format: header 'ntemp ndens'; then ntemp*ndens blocks, each a header
    'dens z temp case ntop ncut' (6 tokens; z = 1 for H I, 2 for He II;
    case is 'A' or 'B'; ncut = 25) followed by ncut(ncut-1)/2 emissivities
    indexed (intrat.f)  k = (ncut-nu)(ncut+nu-1)/2 + nl.  An alpha-tot table
    may follow the blocks; it is ignored because reading stops after
    ntemp*ndens blocks.  Values are 4 pi j / (n_e n_ion) [erg cm^3 s^-1].

    src         : source path (gzip if it ends with .gz, else plain)
    out         : output path
    lines       : list of (n_upper, n_lower, label, wavelength[A])
    want_case   : expected case token 'A' or 'B' (asserted per block)
    Z           : expected z token (asserted per block)
    ion_label   : ion label string (e.g. 'H I', 'He II')
    norm        : normalization string for the header comment
    provenance  : one-line source description for the header comment
    checks      : 'hi' or 'heii' selects the sanity prints
    """
    if src.endswith(".gz"):
        import gzip
        toks = gzip.open(src, "rt").read().split()
    else:
        toks = open(src).read().split()
    p = 0
    ntemp, ndens = int(toks[p]), int(toks[p+1]); p += 2
    temps = np.zeros(ntemp); denss = np.zeros(ndens)
    E = {}
    for ia in range(ntemp):
        for ib in range(ndens):
            dens = float(toks[p]); z = int(toks[p+1]); temp = float(toks[p+2])
            case = toks[p+3]; ncut = int(toks[p+5])
            assert case == want_case, f"unexpected case token {case}"
            assert z == Z, f"unexpected z token {z} (expected {Z})"
            p += 6
            ne_vals = ncut*(ncut-1)//2
            vals = np.array(toks[p:p+ne_vals], dtype=float)
            p += ne_vals
            temps[ia] = temp;  denss[ib] = dens
            for (nu, nl, lab, wl) in lines:
                k = (ncut-nu)*(ncut+nu-1)//2 + nl
                E.setdefault(lab, np.zeros((ntemp, ndens)))[ia, ib] = vals[k-1]

    with open(out, "w") as fh:
        fh.write("# Storey & Hummer (1995, MNRAS 272, 41) " + ion_label +
                 " case-" + want_case + " line emissivities " + norm +
                 " [erg cm^3 s^-1]\n")
        fh.write("# source file: " + provenance + ", parsed by "
                 "tools/fitting/make_sh95_lines.py (2026-07-13)\n")
        fh.write("# grid: NT temperatures [K], ND densities [cm^-3]; then "
                 "one block for each line: label nu nl lambda[A], NT x ND "
                 "values (T rows, n_e columns)\n")
        fh.write(f"GRID {ntemp} {ndens}\n")
        fh.write("T " + " ".join(f"{t:.4e}" for t in temps) + "\n")
        fh.write("NE " + " ".join(f"{d:.4e}" for d in denss) + "\n")
        for (nu, nl, lab, wl) in lines:
            fh.write(f"LINE {lab} {nu} {nl} {wl:.1f}\n")
            for ia in range(ntemp):
                fh.write(" ".join(f"{E[lab][ia, ib]:.4e}"
                                  for ib in range(ndens)) + "\n")
    print(f"wrote {os.path.relpath(out)} ({ntemp} T x {ndens} ne)")
    if checks == "hi":
        # sanity: Hbeta at T=1e4, ne=1e2 should be ~1.24e-25 erg cm^3/s
        ia = np.argmin(abs(temps - 1e4)); ib = np.argmin(abs(denss - 1e2))
        print(f"check: 4pi j(Hbeta)/(ne np) at T={temps[ia]:.0f}, "
              f"ne={denss[ib]:.0f} = {E['Hbeta'][ia, ib]:.4e} "
              f"(expect ~1.24e-25)")
        print(f"check: Halpha/Hbeta = "
              f"{E['Halpha'][ia, ib]/E['Hbeta'][ia, ib]:.4f} (expect ~2.8-2.9)")
    else:
        # sanity at T = 1e4 K: alpha_eff(4686) ~ 3.57e-13 cm^3/s (AGN3), so
        # 4 pi j(4686)/(ne nHeIII) = alpha_eff h nu ~ 1.51e-24 erg cm^3/s;
        # and 1640/4686 ~ 6.6-7.
        ia = np.argmin(abs(temps - 1e4)); ib = np.argmin(abs(denss - 1e4))
        print(f"check: 4pi j(4686)/(ne nHeIII) at T={temps[ia]:.0f}, "
              f"ne={denss[ib]:.0f} = {E['HeII4686'][ia, ib]:.4e} "
              f"(expect ~1.51e-24)")
        print(f"check: 1640/4686 = "
              f"{E['HeII1640'][ia, ib]/E['HeII4686'][ia, ib]:.3f} "
              f"(expect ~6.6-7)")


def main_heii():
    """He II (Z = 2) case B from the ORIGINAL CDS file e2b.d: header
    'ntemp ndens'; blocks 'dens z temp case ntop ncut' followed by
    ncut(ncut-1)/2 emissivities indexed (intrat.f)
        k = (ncut-nu)(ncut+nu-1)/2 + nl,  ncut = 25;
    an alpha-tot table follows the blocks (ignored).  Values are
    4 pi j / (n_e n_HeIII) [erg cm^3 s^-1]."""
    import gzip
    toks = gzip.open(SRC2, "rt").read().split()
    p = 0
    ntemp, ndens = int(toks[p]), int(toks[p+1]); p += 2
    temps = np.zeros(ntemp); denss = np.zeros(ndens)
    E = {}
    for ia in range(ntemp):
        for ib in range(ndens):
            dens = float(toks[p]); temp = float(toks[p+2])
            case = toks[p+3]; ncut = int(toks[p+5])
            assert case == "B", f"unexpected case token {case}"
            p += 6
            ne_vals = ncut*(ncut-1)//2
            vals = np.array(toks[p:p+ne_vals], dtype=float)
            p += ne_vals
            temps[ia] = temp;  denss[ib] = dens
            for (nu, nl, lab, wl) in LINES2:
                k = (ncut-nu)*(ncut+nu-1)//2 + nl
                E.setdefault(lab, np.zeros((ntemp, ndens)))[ia, ib] = vals[k-1]

    with open(OUT2, "w") as fh:
        fh.write("# Storey & Hummer (1995, MNRAS 272, 41) He II (Z=2) case-B "
                 "line emissivities 4 pi j / (n_e n_HeIII) [erg cm^3 s^-1]\n")
        fh.write("# source file: Storey_Hummer/data/e2b.d.gz (original CDS "
                 "table), parsed by tools/fitting/make_sh95_lines.py "
                 "(2026-07-12)\n")
        fh.write("# grid: NT temperatures [K], ND densities [cm^-3]; then "
                 "one block for each line: label nu nl lambda[A], NT x ND "
                 "values (T rows, n_e columns)\n")
        fh.write(f"GRID {ntemp} {ndens}\n")
        fh.write("T " + " ".join(f"{t:.4e}" for t in temps) + "\n")
        fh.write("NE " + " ".join(f"{d:.4e}" for d in denss) + "\n")
        for (nu, nl, lab, wl) in LINES2:
            fh.write(f"LINE {lab} {nu} {nl} {wl:.1f}\n")
            for ia in range(ntemp):
                fh.write(" ".join(f"{E[lab][ia, ib]:.4e}"
                                  for ib in range(ndens)) + "\n")
    print(f"wrote {os.path.relpath(OUT2)} ({ntemp} T x {ndens} ne)")
    # sanity at T = 1e4 K: alpha_eff(4686) ~ 3.57e-13 cm^3/s (AGN3), so
    # 4 pi j(4686)/(ne nHeIII) = alpha_eff h nu ~ 1.51e-24 erg cm^3/s;
    # and 1640/4686 ~ 6.6-7.
    ia = np.argmin(abs(temps - 1e4)); ib = np.argmin(abs(denss - 1e4))
    print(f"check: 4pi j(4686)/(ne nHeIII) at T={temps[ia]:.0f}, "
          f"ne={denss[ib]:.0f} = {E['HeII4686'][ia, ib]:.4e} "
          f"(expect ~1.51e-24)")
    print(f"check: 1640/4686 = "
          f"{E['HeII1640'][ia, ib]/E['HeII4686'][ia, ib]:.3f} (expect ~6.6-7)")


def main_hi_caseA():
    """H I (Z = 1) case A from the ORIGINAL CDS file e1a.d."""
    read_cds(SRC_HI_A, OUT_HI_A, LINES, "A", 1, "H I",
             "4 pi j / (n_e n_p)",
             "Storey_Hummer/data/e1a.d.gz (original CDS table)",
             checks="hi")


def main_heii_caseA():
    """He II (Z = 2) case A from the ORIGINAL CDS file e2a.d."""
    read_cds(SRC_HEII_A, OUT_HEII_A, LINES2, "A", 2, "He II (Z=2)",
             "4 pi j / (n_e n_HeIII)",
             "Storey_Hummer/data/e2a.d.gz (original CDS table)",
             checks="heii")


if __name__ == "__main__":
    main()
    main_heii()
    main_hi_caseA()
    main_heii_caseA()
