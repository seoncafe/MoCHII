#!/usr/bin/env python3
"""Collect the fitted line-cooling coefficients into one text table.

Reads each selected model's primary cooling files and writes
``data/atomic/excitation/<model>/cooling_fit_parameters.txt``: a readable listing
of the coefficients of

    Lambda_ion(T) = T^(-1/2) * sum_i A_i * exp(-T_i / T)   [erg cm^3 s^-1]

so the whole fitted set can be read, cited, or re-implemented without opening
the individual files.  The individual files remain the primary source; this table
is generated from them and carries no independent numbers.
"""
import argparse
import glob
import io
import os
import re

from atomic_output import atomic_write_text

HERE = os.path.dirname(os.path.abspath(__file__))
ATOMIC = os.path.join(HERE, "..", "..", "data", "atomic")

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
            m_source = re.search(r"Numerical fit interval:\s*([\d.eE+-]+)\.\.([\d.eE+-]+)\s*K", line)
            if m_source:
                lo, hi = float(m_source.group(1)), float(m_source.group(2))
            m_error = re.search(r"Maximum relative error on \d+ temperatures:\s*([\d.eE+-]+)", line)
            if m_error:
                err = 100.0 * float(m_error.group(1).rstrip('.'))
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=('chianti', 'mocassin_v2'), default='chianti')
    args = parser.parse_args()
    atomic = os.path.join(ATOMIC, 'excitation', args.model)
    output = os.path.join(atomic, 'cooling_fit_parameters.txt')
    files = sorted(glob.glob(os.path.join(atomic, "cooling_*.txt")))
    files = [f for f in files if "fit_parameters" not in f]
    rows = []
    for f in files:
        tag = os.path.basename(f)[len("cooling_"):-len(".txt")]
        nterm, terms, rng, err = read_fit(f)
        if (nterm != len(terms) or not terms or None in rng or err is None):
            raise ValueError(f"Incomplete cooling coefficients or fit metadata: {f}")
        rows.append((ion_label(tag), tag, nterm, terms, rng, err))
    rows.sort(key=lambda r: (r[1].rsplit("_", 1)[0], int(r[1].rsplit("_", 1)[1])))

    with io.StringIO() as fh:
        w = fh.write
        w("# MoCHII: fitted collisional line-cooling coefficients (all ions)\n")
        w("# Selected metal excitation model: %s\n" % args.model)
        w("#\n")
        w("#   Lambda_ion(T) = T^(-1/2) * sum_{i=1}^{N} A_i * exp(-T_i / T)\n")
        w("#\n")
        w("# Lambda_ion is a rate coefficient in erg cm^3 s^-1; the volumetric\n")
        w("# cooling rate of the ion is n_e * n_ion * Lambda_ion(T) [erg cm^-3 s^-1].\n")
        w("# A_i is in erg cm^3 s^-1 K^(1/2), T_i in K, T in K.\n")
        w("#\n")
        w("# Sources: identified excitation data recorded in each primary cooling file,\n")
        w("# including CHIANTI, Stout, and corrected supplied MOCASSIN models, in the optically thin n_e -> 0 limit of\n")
        w("# stationary collisions and radiative cascades among stable states.\n")
        w("# Bound collision-connected states use the checked physical energy selection.\n")
        w("# Nonnegative nodes use cubic interpolation unless it overshoots below zero;\n")
        w("# those channels use PCHIP through the same source nodes.\n")
        w("# Identified Stout temperature tables use log-PCHIP in T and Upsilon,\n")
        w("# with the adopted constant first-node Upsilon below the source interval.\n")
        w("# Coefficients describe the total fit and can be signed. No finite-density\n")
        w("# suppression is included. Fit provenance is recorded in each primary file; collected by\n")
        w("# tools/fitting/make_cooling_table.py from the cooling_<ion>.txt files,\n")
        w("# which remain the primary source (they carry the full provenance header).\n")
        w("#\n")
        w("# Columns:  ion  tag  N  fit_range_K  max_err_%  then N lines of  A_i  T_i\n")
        w("#\n")
        for label, tag, nterm, terms, rng, err in rows:
            w("\n%-8s %-8s N=%d  range= %.1e - %.1e K  max_err= %.2f%%\n"
              % (label, tag, len(terms), rng[0], rng[1], err))
            for A, T in terms:
                # Retain the source file's double-precision coefficients;
                # signed sums can amplify rounding of individual terms.
                w("    %.16e  %.16e\n" % (A, T))
        atomic_write_text(output, fh.getvalue())
    print("wrote", os.path.normpath(output))
    print("ions:", len(rows), " terms:", sum(len(r[3]) for r in rows))


if __name__ == "__main__":
    main()
