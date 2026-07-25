#!/usr/bin/env python3
"""Collect the fitted line-cooling coefficients into one text table.

Reads every ``data/atomic/cooling_<ion>.txt`` and writes
``data/atomic/cooling_fit_parameters.txt``: a single, human-readable listing
of the coefficients of

    Lambda_ion(T) = T^(-1/2) * sum_i A_i * exp(-T_i / T)   [erg cm^3 s^-1]

so the whole fitted set can be read, cited, or re-implemented without opening
the individual files.  The individual files remain the primary source; this table
is generated from them and carries no independent numbers.
"""
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ATOMIC = os.path.join(HERE, "..", "..", "data", "atomic")
OUT = os.path.join(ATOMIC, "cooling_fit_parameters.txt")

ROMAN = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII"]


def ion_label(tag):
    """'o_3' -> 'O III'."""
    el, stage = tag.rsplit("_", 1)
    return "%s %s" % (el.capitalize(), ROMAN[int(stage)])


def read_fit(path):
    """Return (n_terms, [(A_i, T_i)], fit range, max error) of one ion file."""
    lo = hi = err = None
    terms = []
    nterm = None
    for line in open(path):
        if line.startswith("#"):
            m = re.search(r"fit range:\s*([\d.eE+-]+)\s*-\s*([\d.eE+-]+)\s*K", line)
            if m:
                lo, hi = float(m.group(1)), float(m.group(2))
            m = re.search(r"max fit error:\s*([\d.]+)%", line)
            if m:
                err = float(m.group(1))
            continue
        tok = line.split()
        if not tok:
            continue
        if nterm is None and len(tok) == 1:
            nterm = int(tok[0])
        elif len(tok) >= 2:
            terms.append((float(tok[0]), float(tok[1])))
    return nterm, terms, (lo, hi), err


def main():
    files = sorted(glob.glob(os.path.join(ATOMIC, "cooling_*.txt")))
    files = [f for f in files if "fit_parameters" not in f]
    rows = []
    for f in files:
        tag = os.path.basename(f)[len("cooling_"):-len(".txt")]
        nterm, terms, rng, err = read_fit(f)
        rows.append((ion_label(tag), tag, nterm, terms, rng, err))
    rows.sort(key=lambda r: (r[1].rsplit("_", 1)[0], int(r[1].rsplit("_", 1)[1])))

    with open(OUT, "w") as fh:
        w = fh.write
        w("# MoCHII: fitted collisional line-cooling coefficients (all ions)\n")
        w("#\n")
        w("#   Lambda_ion(T) = T^(-1/2) * sum_{i=1}^{N} A_i * exp(-T_i / T)\n")
        w("#\n")
        w("# Lambda_ion is a rate coefficient in erg cm^3 s^-1; the volumetric\n")
        w("# cooling rate of the ion is n_e * n_ion * Lambda_ion(T) [erg cm^-3 s^-1].\n")
        w("# A_i is in erg cm^3 s^-1 K^(1/2), T_i in K, T in K.\n")
        w("#\n")
        w("# Source: CHIANTI v11.0.2 level energies and effective collision strengths\n")
        w("# (Burgess-Tully descaling) in the optically thin n_e -> 0 limit of\n")
        w("# statistical equilibrium: all ions in the lowest level, excitation out of\n")
        w("# it only, each excitation radiating the full excitation energy through the\n")
        w("# cascade; transition energies are the observed level differences.  No\n")
        w("# density suppression above n_crit.  Fitted by\n")
        w("# tools/fitting/fit_cooling.py; collected here by\n")
        w("# tools/fitting/make_cooling_table.py from the cooling_<ion>.txt files,\n")
        w("# which remain the primary source (they carry the full provenance header).\n")
        w("#\n")
        w("# Columns:  ion  tag  N  fit_range_K  max_err_%  then N lines of  A_i  T_i\n")
        w("#\n")
        for label, tag, nterm, terms, rng, err in rows:
            w("\n%-8s %-8s N=%d  range= %.1e - %.1e K  max_err= %.2f%%\n"
              % (label, tag, len(terms), rng[0], rng[1], err))
            for A, T in terms:
                w("    %16.8e  %14.6e\n" % (A, T))
    print("wrote", os.path.normpath(OUT))
    print("ions:", len(rows), " terms:", sum(len(r[3]) for r in rows))


if __name__ == "__main__":
    main()
