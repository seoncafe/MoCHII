#!/usr/bin/env python3
"""He I case-B collisional-radiative line emissivities (Porter et al. 2012, 2013).

These are full case-B collisional-radiative totals, not pure recombination
cascades: the tabulated 4 pi j / (n_e n_He+) is density-dependent because it
already includes collisional excitation out of the 2^3S metastable (the 10833
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

# Principal He I diagnostics.  The keys are the AIR wavelengths that index
# the Cloudy source file; the written wavelength and the label are VACUUM, the
# convention of every other block of the MoCHII line table (the SH95 H I/He II
# tables carry vacuum values and the n-level wavelengths come from level energy
# differences, which are vacuum by construction).  Label <= 12 characters.
WANT = {
    "3888.63": ("HeI3890", 3, 2),
    "4026.20": ("HeI4027", 5, 2),
    "4471.49": ("HeI4473", 4, 2),
    "5875.64": ("HeI5877", 3, 2),
    "6678.15": ("HeI6680", 3, 2),
    "7065.22": ("HeI7067", 3, 2),
    "7281.35": ("HeI7283", 3, 2),
    "10830.25": ("HeI10833", 2, 2),
}


def air_to_vacuum(lam_air):
    """Vacuum wavelength [A] of an air wavelength [A].

    Inverse of the IAU standard dispersion of air (Morton 2000, ApJS 130,
    403); accurate well below the 0.01 A the tabulated values carry.  Check:
    10830.25 -> 10833.217, matching the NIST fine-structure components stored
    in data/atomic/hei_metastable.txt (10833.22 for the J = 1 line).
    """
    s2 = (1.0e4/lam_air)**2
    n = (1.0 + 0.00008336624212083
         + 0.02408926869968/(130.1065924800 - s2)
         + 0.0001599740894897/(38.92568793293 - s2))
    return lam_air*n


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
        fh.write("# VACUUM wavelengths (the source file is indexed by air "
                 "wavelength, Cloudy convention; converted here so this block "
                 "matches the rest of the MoCHII line table).  Air -> vacuum: "
                 "3888.63->3889.73, 4026.20->4027.34, 4471.49->4472.75, "
                 "5875.64->5877.27, 6678.15->6679.99, 7065.22->7067.17, "
                 "7281.35->7283.36, 10830.25->10833.22\n")
        fh.write("# grid: NT temperatures [K], ND densities [cm^-3]; then "
                 "LINE label nu nl lambda[A], NT x ND values "
                 "(T rows, n_e columns)\n")
        fh.write(f"GRID {nT} {nD}\n")
        fh.write("T  " + " ".join(f"{t:.4e}" for t in temps) + "\n")
        fh.write("NE " + " ".join(f"{10.0**d:.4e}" for d in dens) + "\n")
        # emit in ascending wavelength for readability
        for wl in sorted(WANT, key=float):
            if wl not in blocks:
                continue
            label, nu, nl = WANT[wl]
            fh.write(f"LINE {label} {nu} {nl} "
                     f"{air_to_vacuum(float(wl)):.2f}\n")
            g = blocks[wl]                     # [density][temperature]
            for it in range(nT):               # transpose -> T rows, n_e cols
                fh.write(" ".join(f"{10.0**g[idn][it]:.4e}" for idn in range(nD))
                         + "\n")
    print(f"wrote {os.path.relpath(OUT)} ({len(blocks)} lines, {nT} T x {nD} ne)")


if __name__ == "__main__":
    main()
