#!/usr/bin/env python3
"""Generate tabulated Planck spectrum files for the source-spectrum gates.

Both files use the SAME shape the code's planck_nu evaluates,
    B_E ~ E^3 / (exp(E * ev2erg / (kB * T)) - 1)     [arbitrary normalization],
on 400 log-spaced points over 10-110 eV (the ionizing band 13.598-100 eV sits
strictly inside, so no zero-outside effect enters).  Each column is normalized
independently inside the code (source i to src_lum(i); the external field to
pi*J*A), so the absolute scale here is immaterial - only the shape matters.

Writes:
  tests/multi_src/two_temp_spec.txt  3 cols: E, B_E(3e4 K), B_E(5e4 K)
                                     (feeds src_spectrum_file of two_temp_file.in)
  tests/ext_field/ext_spec_4e4.txt   2 cols: E, B_E(4e4 K)
                                     (feeds ext_spectrum of ext_planckfile.in)
"""
import os
import numpy as np

EV2ERG = 1.602176634e-12      # define.f90 ev2erg
KB_CGS = 1.380649e-16         # define.f90 kboltz_cgs


def planck_nu(E, T):
    """Code's planck_nu shape (arbitrary normalization), E [eV], T [K]."""
    x = E * EV2ERG / (KB_CGS * T)
    return E**3 / np.expm1(x)


def write_table(path, E, cols, labels):
    with open(path, "w") as f:
        f.write("# tabulated Planck B_E (arbitrary normalization); shape = "
                "E^3/(exp(E*ev2erg/(kB*T))-1)\n")
        f.write("# col1 = E [eV]; " + "; ".join(labels) + "\n")
        for i in range(E.size):
            f.write(f"{E[i]:.8e}" + "".join(f" {c[i]:.8e}" for c in cols) + "\n")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    E = np.logspace(np.log10(10.0), np.log10(110.0), 400)

    # two-temperature multi-column file (sources 1, 2 of two_temp_file.in)
    write_table(os.path.join(here, "two_temp_spec.txt"), E,
                [planck_nu(E, 3.0e4), planck_nu(E, 5.0e4)],
                ["col2 = B_E(3e4 K)", "col3 = B_E(5e4 K)"])

    # external 4e4 K 2-column file (ext_planckfile.in)
    extdir = os.path.join(here, "..", "ext_field")
    write_table(os.path.join(extdir, "ext_spec_4e4.txt"), E,
                [planck_nu(E, 4.0e4)], ["col2 = B_E(4e4 K)"])

    print("wrote two_temp_spec.txt and ../ext_field/ext_spec_4e4.txt "
          f"({E.size} points, {E[0]:.3f}-{E[-1]:.3f} eV)")


if __name__ == "__main__":
    main()
