#!/usr/bin/env python3
"""Atomic-hydrogen impact rates for the ground-term fine-structure levels of
C I, C II and O I: data/atomic/hcoll_<el>_<stage>.txt.

nlevel_mod reads hcoll_<el>_<stage>.txt beside nlevel_<el>_<stage>.txt and
adds n_HI times these rates to the electron-impact rates of the same
statistical equilibrium, so the level populations -- hence both the line
emissivities written at output time and the cooling increment the thermal
balance adds (species_mod::metal_cooling_H) -- carry both colliders in one
solve, and the two saturate together.

Sources
-------
C I and O I: Abrahamsson, Krems & Dalgarno (2007), ApJ 654, 1171 --
  quantum close-coupling on the Kalemos et al. (1999) CH and Parlant &
  Yarkony (1999) OH potentials, read from the journal PDF (2026-09-24).
  Their eq. (1) gives the EXCITATION rate coefficient j -> j',
      ln(k / 1e-11 cm^3 s^-1) = sum_{i=0}^{N} a_i T^(-i b),
  with a_i and b from Table 2 (C + H; valid 5-1000 K) and Table 4 (O + H;
  valid 30-1000 K).  In the text layer of the PDF the minus signs of the
  exponent and of every other coefficient are glyphs that plain text
  extraction drops; the signs below were recovered from the glyph stream
  and are confirmed by the self-check, which reproduces the paper's own
  Table 1 ("Accurate" column) and Table 3 ("This work") to the stated fit
  accuracy.  These replace, for O I, the 9.2e-11 (T/100)^0.67 [O I] 63 um
  de-excitation rate of the earlier two-level channel (no source recorded;
  0.20-0.46x the present rate over 50-1000 K), and add C I, which had no
  H-impact channel at all, and the [O I] 145 um / [C I] 609 and 370 um
  ground-term transitions the two-level treatment could not carry.
C II: Barinovs, van Hemert, Krems & Dalgarno (2005), ApJ 620, 537, as the
  power law q_ul = 7.6e-10 (T/100)^0.14 cm^3 s^-1 already used for
  [C II] 158 um + H (the paper itself was not re-read here; the law agrees
  with an independent fit to the LAMDA c+.dat tabulation of the same
  calculation to 2.5% over 60-1000 K).

Evaluation (nlevel_mod::hcoll_rate_matrix): the DE-excitation rate is
evaluated at T clamped to [Tmin, Tmax], the range of the calculation --
held at the endpoint rather than extrapolated, since a polynomial in
T^(-b) or a power law is not a statement about temperatures outside it --
and the excitation rate follows by detailed balance at the actual T with
the level energies of the n-level atom, so detailed balance holds at every
temperature.  Above 1000 K (2000 K for C II) the H-impact rate is therefore
held fixed; where gas that warm is also partly ionized, electrons dominate
the excitation anyway.

Run:  python3 make_hcoll_fine_structure.py [--output-dir DIR]
The level indices are checked against nlevel_<el>_<stage>.txt in the output
directory (statistical weight and energy of each addressed level) before a
file is written, so a rate cannot be attached to the wrong level.
"""
import argparse
from datetime import datetime
import os
import math
from zoneinfo import ZoneInfo

DATE = datetime.now(ZoneInfo("Asia/Seoul")).date().isoformat()
# Ground-term levels the rows address: index -> (g, energy [cm^-1], label),
# energies from NIST ASD (C I 3P 0/16.40/43.40, C II 2P 0/63.42,
# O I 3P 0/158.265/226.977).
GROUND_TERM_LEVELS = {
    "c_1": {1: (1.0, 0.0, "3P0"), 2: (3.0, 16.40, "3P1"), 3: (5.0, 43.40, "3P2")},
    "c_2": {1: (2.0, 0.0, "2P1/2"), 2: (4.0, 63.42, "2P3/2")},
    "o_1": {1: (5.0, 0.0, "3P2"), 2: (3.0, 158.265, "3P1"), 3: (1.0, 226.977, "3P0")},
}
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "..", "data", "atomic")

# Abrahamsson et al. (2007) eq. (1) coefficients a_0..a_N and b.
ABR_C = {  # Table 2, C(3P) + H; levels of nlevel_c_1: 1=3P0 2=3P1 3=3P2
    (1, 2): ([3.6593, 56.6023, -802.9765, 5025.1882, -17874.4255,
              38343.6655, -49249.4895, 34789.3941, -10390.9809], 1.0/4.0),
    (1, 3): ([10.8377, -173.4153, 2024.0272, -13391.6549, 52198.5522,
              -124518.3586, 178182.5823, -140970.6106, 47504.5861], 1.0/3.0),
    (2, 3): ([15.8996, -201.3030, 1533.6164, -6491.0083, 15921.9239,
              -22691.1632, 17334.7529, -5517.9360, 0.0], 1.0/4.0),
}
ABR_O = {  # Table 4, O(3P) + H; levels of nlevel_o_1: 1=3P2 2=3P1 3=3P0
    (1, 2): ([4.581, -156.118, 2679.979, -78996.962, 1308323.468,
              -13011761.861, 71010784.971, -162826621.855], 3.0/4.0),
    (1, 3): ([3.297, -168.382, 1844.099, -68362.889, 1376864.737,
              -17964610.169, 134374927.808, -430107587.886], 3.0/4.0),
    (2, 3): ([3.437, 17.443, -618.761, 3757.156, -12736.468, 22785.266,
              -22759.228, 12668.261], 1.0/2.0),
}
# The paper's tabulated rates [1e-11 cm^3 s^-1] for the self-check.
TAB1_C = {  # Table 1 "Accurate": T -> (k0->1, k0->2, k1->2)
    10: (3.68, 0.065, 0.678), 20: (11.08, 1.39, 4.61),
    40: (19.44, 7.035, 13.02), 100: (29.75, 23.83, 29.79),
    200: (37.56, 41.19, 45.20), 300: (42.90, 51.70, 54.79),
    500: (50.68, 65.80, 68.12), 1000: (64.15, 88.64, 90.61)}
TAB3_O = {  # Table 3 "This work": T -> (k2->1, k2->0, k1->0)
    50: (0.189, 0.0080, 1.10), 70: (0.740, 0.0548, 2.70),
    100: (2.20, 0.243, 5.38), 150: (5.49, 0.827, 9.89),
    200: (9.06, 1.59, 13.80), 300: (15.83, 3.24, 19.65),
    500: (26.79, 6.18, 26.21), 700: (34.91, 8.47, 29.49),
    1000: (43.75, 11.01, 32.03)}
# The published fit reproduces its own tables to <= 1% everywhere except at
# the cold edge of the O + H range: at 50 K the eq. (1) form gives k1->0
# 9.1% above Table 3 (1.20 vs 1.10e-11) and k2->1 1.7% below.  That is the
# paper's fit, not a transcription error (70 K already agrees to 0.1%), so
# 50 K is printed but held out of the tolerance.
FIT_EDGE_O = 50


def k_abr(coef, b, T):
    """Abrahamsson et al. (2007) eq. (1): excitation rate [cm^3 s^-1]."""
    return 1.0e-11 * math.exp(sum(a * T ** (-i * b) for i, a in enumerate(coef)))


def self_check():
    worst = 0.0
    for lab, table, fits, edge in (("C I", TAB1_C, ABR_C, None),
                                   ("O I", TAB3_O, ABR_O, FIT_EDGE_O)):
        for T, ref in table.items():
            devs = [k_abr(*fits[key], T) / 1e-11 / r - 1.0
                    for key, r in zip([(1, 2), (1, 3), (2, 3)], ref)]
            if T != edge:
                worst = max(worst, max(abs(d) for d in devs))
            print(f"  {lab}  T={T:5d}  fit/table - 1 = "
                  + "  ".join(f"{d:+.4f}" for d in devs)
                  + ("   (fit edge, not tested)" if T == edge else ""))
    print(f"  worst |fit/table - 1| = {worst:.4f}")
    assert worst < 0.012, worst


def verify_ground_term_indices(outdir, tag):
    """Check the addressed levels of nlevel_<tag>.txt in outdir.

    Each index must carry the expected statistical weight and an energy
    within 0.5 cm^-1 of the NIST ground-term value.
    """
    path = os.path.join(outdir, f"nlevel_{tag}.txt")
    records = [line.split() for line in open(path)
               if line.strip() and not line.lstrip().startswith("#")]
    start = next(i for i, r in enumerate(records) if r[0] == "NLEV")
    levels = {int(r[0]): (float(r[1]), float(r[2]))
              for r in records[start+1:start+1+int(records[start][1])]}
    for index, (g, energy, label) in GROUND_TERM_LEVELS[tag].items():
        found = levels.get(index)
        if found is None or found[0] != g or abs(found[1] - energy) > 0.5:
            raise ValueError(f"{path}: level {index} is {found}, expected {label} "
                             f"g={g:g} E={energy} cm^-1")
        print(f"  {tag} level {index}: g={found[0]:g} E={found[1]:.3f} cm^-1 = {label}")


def write(name, header, rows, outdir=OUTDIR):
    verify_ground_term_indices(outdir, name.removeprefix("hcoll_").removesuffix(".txt"))
    path = os.path.join(outdir, name)
    with open(path, "w") as fh:
        for line in header:
            fh.write("# " + line + "\n")
        fh.write(f"# generated by tools/fitting/make_hcoll_fine_structure.py "
                 f"on {DATE}\n")
        fh.write(f"NHCOL {len(rows)}\n")
        for r in rows:
            fh.write(" ".join(r) + "\n")
    print(f"wrote {os.path.relpath(path)}")


FORMAT = [
    "Row: lo hi form Tmin Tmax c_1 .. c_10  (levels as in nlevel_<el>_<st>.txt,",
    "  lo the lower level).  form 1: Abrahamsson et al. (2007) eq. (1),",
    "  EXCITATION lo->hi: ln(k/1e-11) = sum_i c_(i+2) T^(-i c_1), i = 0..8.",
    "  form 2: DE-excitation hi->lo, q = c_1 (T/100)^c_2.  The de-excitation",
    "  rate is held at T clamped to [Tmin, Tmax]; excitation by detailed",
    "  balance at the actual T.  Rates in cm^3 s^-1.",
]


def rows_abr(table, tmin, tmax):
    out = []
    for (lo, hi), (coef, b) in table.items():
        c = [b] + coef + [0.0] * (9 - len(coef))
        out.append([f"{lo:d}", f"{hi:d}", "1", f"{tmin:.1f}", f"{tmax:.1f}"]
                   + [f"{x:.10e}" for x in c])
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=OUTDIR,
                        help="directory holding the nlevel atoms the rows address")
    outdir = parser.parse_args(argv).output_dir
    self_check()
    write("hcoll_o_1.txt",
          ["O I ground-term fine-structure excitation by atomic hydrogen:",
           "Abrahamsson, Krems & Dalgarno (2007), ApJ 654, 1171, Table 4",
           "(eq. 1; valid 30-1000 K).  3P2->3P1 is [O I] 63 um, 3P1->3P0",
           "[O I] 145 um."] + FORMAT,
          rows_abr(ABR_O, 30.0, 1000.0), outdir)
    write("hcoll_c_1.txt",
          ["C I ground-term fine-structure excitation by atomic hydrogen:",
           "Abrahamsson, Krems & Dalgarno (2007), ApJ 654, 1171, Table 2",
           "(eq. 1; valid 5-1000 K).  3P0->3P1 is [C I] 609 um, 3P1->3P2",
           "[C I] 370 um."] + FORMAT,
          rows_abr(ABR_C, 5.0, 1000.0), outdir)
    write("hcoll_c_2.txt",
          ["C II 2P1/2-2P3/2 ([C II] 158 um) de-excitation by atomic",
           "hydrogen: Barinovs et al. (2005), ApJ 620, 537, as the power law",
           "7.6e-10 (T/100)^0.14, held outside 20-2000 K."] + FORMAT,
          [["1", "2", "2", "20.0", "2000.0",
            f"{7.6e-10:.10e}", f"{0.14:.10e}"] + [f"{0.0:.10e}"] * 8], outdir)


if __name__ == "__main__":
    main()
