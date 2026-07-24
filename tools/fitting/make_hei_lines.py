#!/usr/bin/env python3
"""He I case-B collisional-radiative line emissivities (Porter et al. 2012, 2013).

These are full case-B collisional-radiative totals, not pure recombination
cascades: the tabulated 4 pi j / (n_e n_He+) is density-dependent because it
already includes collisional excitation out of the 2^3S metastable (the 10830
coefficient rises 6.6x from n_e=10 to n_e=1e4 cm^-3 at 1e4 K).

Parse Cloudy c23.01's data/he1_case_b.dat (the Porter grid: log10 of
4 pi j / (n_e n_He+) [erg cm^3 s^-1] on a 21-T x 14-n_e grid) and write the
principal H II-region diagnostic lines to data/atomic/hei_porter_caseB.txt in
the same GRID/T/NE/LINE layout the SH95 tables use (linear emissivities; the
Fortran reader interpolates in log space).  Copy he1_case_b.dat into
data/atomic once for a self-contained source, like the Badnell tables.
"""
import os

ATOM = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "data", "atomic")
# he1_case_b.dat copied into data/atomic (from Cloudy c23.01) for a
# self-contained source, like the Badnell tables.
SRC = os.path.join(ATOM, "he1_case_b.dat")
OUT = os.path.join(ATOM, "hei_porter_caseB.txt")

# principal He I diagnostics (air wavelengths, Cloudy convention); label <= 12
WANT = {
    "3888.63": ("HeI3889", 3, 2),
    "4026.20": ("HeI4026", 5, 2),
    "4471.49": ("HeI4471", 4, 2),
    "5875.64": ("HeI5876", 3, 2),
    "6678.15": ("HeI6678", 3, 2),
    "7065.22": ("HeI7065", 3, 2),
    "7281.35": ("HeI7281", 3, 2),
    "10830.25": ("HeI10830", 2, 2),
}


def main():
    lines = open(SRC).read().splitlines()
    # header: line 0 = "<magic>\t<nlines>"; density row and temperature row
    # follow their "#These are the ..." comment lines.
    dens = temps = None
    i = 0
    while i < len(lines):
        s = lines[i]
        if s.startswith("#These are the logs of densities"):
            dens = [float(x) for x in lines[i + 1].split()]
            i += 2; continue
        if s.startswith("#These are the temperatures"):
            temps = [float(x) for x in lines[i + 1].split()]
            i += 2; continue
        if dens is not None and temps is not None:
            break
        i += 1
    nD, nT = len(dens), len(temps)

    # collect the wanted line blocks: header "#He 1  <wl>A", then a
    # "wl QN_lo QN_up" row, then nD rows of (idx + nT values).
    blocks = {}
    j = 0
    while j < len(lines):
        if lines[j].startswith("#He 1"):
            hdr = lines[j + 1].split()
            wl = hdr[0]
            if wl in WANT:
                grid = []                      # [density][temperature] log10
                for d in range(nD):
                    row = lines[j + 2 + d].split()
                    grid.append([float(x) for x in row[1:1 + nT]])
                blocks[wl] = grid
            j += 2 + nD
        else:
            j += 1

    with open(OUT, "w") as fh:
        fh.write("# He I case-B collisional-radiative line emissivities "
                 "4 pi j / (n_e n_He+) [erg cm^3 s^-1] (density-dependent)\n")
        fh.write("# Porter et al. (2012 ApJL 756 L14; 2013 erratum), from "
                 "Cloudy c23.01 data/he1_case_b.dat; parsed by "
                 "tools/fitting/make_hei_lines.py\n")
        fh.write("# air wavelengths (Cloudy convention); grid: NT temperatures "
                 "[K], ND densities [cm^-3]; then LINE label nu nl lambda[A], "
                 "NT x ND values (T rows, n_e columns)\n")
        fh.write(f"GRID {nT} {nD}\n")
        fh.write("T  " + " ".join(f"{t:.4e}" for t in temps) + "\n")
        fh.write("NE " + " ".join(f"{10.0**d:.4e}" for d in dens) + "\n")
        # emit in ascending wavelength for readability
        for wl in sorted(WANT, key=float):
            if wl not in blocks:
                continue
            label, nu, nl = WANT[wl]
            fh.write(f"LINE {label} {nu} {nl} {wl}\n")
            g = blocks[wl]                     # [density][temperature]
            for it in range(nT):               # transpose -> T rows, n_e cols
                fh.write(" ".join(f"{10.0**g[idn][it]:.4e}" for idn in range(nD))
                         + "\n")
    print(f"wrote {os.path.relpath(OUT)} ({len(blocks)} lines, {nT} T x {nD} ne)")


if __name__ == "__main__":
    main()
