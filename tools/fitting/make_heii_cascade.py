#!/usr/bin/env python3
"""Write data/atomic/heii_caseB_n2.txt: what a He III case-B recombination
does when it reaches n = 2 of He II.

In case B every recombination that is not to the ground level cascades down to
n = 2 and leaves through 2^2S (two-photon, 40.813 eV shared by a pair) or
2^2P (Ly-alpha, 40.813 eV in one photon).  Both channels photoionize hydrogen,
so the diffuse field needs the branching

    f_2s(T) = alpha_eff(2^2S) / alpha_B ,      f_2p = 1 - f_2s .

Both numerator and denominator are taken from Storey & Hummer (1995), so the
ratio is a closure of ONE cascade calculation and not a difference of
determinations from different sources -- the mistake that put a 2.4% He I
error into the diffuse packet count when alpha_B was formed as alpha_A minus a
Mao level-resolved alpha_1.

  alpha_B    : the alpha-tot table that follows the emissivity blocks of the
               original CDS file e2b.d (Z = 2, case B), read here at
               n_e = 100 cm^-3.  intrat.f calls it 'alpha-tot'.
  alpha(2S)  : the SH95 n_e = 100 values for the 2^2S level, the same four
               anchors MOCASSIN's emission_mod carries and MoCHII has been
               using for the written He II two-photon continuum.  They cover
               5000-30000 K; the Fortran clamps outside that range, so f_2s is
               held at its endpoint values there.

Cross-check reported at the end: MOCASSIN's independent He II Ly-alpha
emissivity fit, 10^(-0.792 log10 Te + 6.01) x 1e-25 erg cm^3 s^-1, which is a
different determination of the same 2^2P branch.
"""
import gzip
import os

import numpy as np

SRC = os.path.expanduser("~/RT_Codes/Storey_Hummer/data/e2b.d.gz")
OUT = os.path.join(os.path.dirname(__file__), "..", "..",
                   "data", "atomic", "heii_caseB_n2.txt")

NE_WANT = 1.0e2          # cm^-3

# SH95 alpha_eff(2^2S) for He II at n_e = 100 cm^-3, in 1e-13 cm^3 s^-1.
ALPHA_2S_T = np.array([5.0e3, 1.0e4, 1.5e4, 3.0e4])
ALPHA_2S_V = np.array([6.161, 4.091, 3.189, 2.035]) * 1.0e-13

EV = 1.602176634e-12      # erg
E_LYA = 40.813            # eV, He II 2 -> 1


def read_alpha_tot(src):
    """alpha-tot [cm^3 s^-1] on the (T, n_e) grid of the CDS file e2b.d.

    The emissivity blocks come first: ntemp*ndens of them, each a 6-token
    header (dens z temp case ntop ncut) followed by ncut(ncut-1)/2 values.
    What follows is the alpha-tot table, ntemp*ndens values laid out with
    temperature varying slowest.
    """
    toks = gzip.open(src, "rt").read().split()
    p = 0
    ntemp, ndens = int(toks[p]), int(toks[p + 1])
    p += 2
    temps = np.zeros(ntemp)
    denss = np.zeros(ndens)
    for ia in range(ntemp):
        for ib in range(ndens):
            denss[ib] = float(toks[p])
            assert int(toks[p + 1]) == 2, "expected Z = 2"
            temps[ia] = float(toks[p + 2])
            assert toks[p + 3] == "B", "expected case B"
            ncut = int(toks[p + 5])
            p += 6 + ncut * (ncut - 1) // 2
    rest = np.array(toks[p:], dtype=float)
    assert rest.size == ntemp * ndens, (
        f"alpha-tot block has {rest.size} values, expected {ntemp*ndens}")
    return temps, denss, rest.reshape(ntemp, ndens)


def loglog(T, xt, xv):
    """Log-log interpolation with endpoint clamping -- the rule the Fortran
    applies, so the file and the code agree by construction."""
    lt = np.log10(np.clip(T, xt[0], xt[-1]))
    return 10.0 ** np.interp(lt, np.log10(xt), np.log10(xv))


def main():
    temps, denss, atot = read_alpha_tot(SRC)
    ib = int(np.argmin(np.abs(denss - NE_WANT)))
    assert abs(denss[ib] - NE_WANT) / NE_WANT < 1e-6, "n_e = 100 column absent"
    aB = atot[:, ib]

    # sanity: alpha_B must fall with T and sit near the textbook 1.5e-12 at 1e4
    assert np.all(np.diff(aB) < 0.0), "alpha_B is not monotonic in T"
    i4 = int(np.argmin(np.abs(temps - 1.0e4)))
    print(f"alpha_B(He II, 1e4 K, n_e = 100) = {aB[i4]:.4e} cm^3/s "
          f"(AGN3 quotes ~1.55e-12)")

    with open(OUT, "w") as fh:
        fh.write("# MoCHII: He II n = 2 exit branching for the He III case-B "
                 "cascade\n")
        fh.write("# alpha_B: Storey & Hummer (1995, MNRAS 272, 41) alpha-tot, "
                 "Z = 2, case B, n_e = 100 cm^-3,\n")
        fh.write("#          from the block following the emissivities of the "
                 "CDS file e2b.d\n")
        fh.write("# alpha_2s: SH95 alpha_eff(2^2S), Z = 2, case B, "
                 "n_e = 100 cm^-3 (4 anchors, 5000-30000 K)\n")
        fh.write("# f_2s = alpha_2s/alpha_B is a closure of ONE cascade "
                 "calculation; f_2p = 1 - f_2s\n")
        fh.write("# written by tools/fitting/make_heii_cascade.py\n")
        fh.write("# columns: T [K], alpha [cm^3 s^-1]; log-log interpolation, "
                 "clamped outside\n")
        fh.write(f"ALPHA_B {len(temps)}\n")
        for t, a in zip(temps, aB):
            fh.write(f"{t:12.4e} {a:12.4e}\n")
        fh.write(f"ALPHA_2S {len(ALPHA_2S_T)}\n")
        for t, a in zip(ALPHA_2S_T, ALPHA_2S_V):
            fh.write(f"{t:12.4e} {a:12.4e}\n")
    print(f"wrote {os.path.normpath(OUT)}")

    # --- report the branching and cross-check f_2p against MOCASSIN's
    # --- independent He II Ly-alpha emissivity fit
    print(f"\n{'T [K]':>8} {'alpha_B':>11} {'alpha_2s':>11} {'f_2s':>7} "
          f"{'f_2p':>7} {'alpha_2p':>11} {'MOCASSIN':>11} {'ratio':>7}")
    for T in (5.0e3, 1.0e4, 1.5e4, 2.0e4, 3.0e4):
        b = loglog(T, temps, aB)
        s = loglog(T, ALPHA_2S_T, ALPHA_2S_V)
        f2s = s / b
        a2p = b - s
        j_moc = 10.0 ** (-0.792 * np.log10(T) + 6.01) * 1.0e-25
        a2p_moc = j_moc / (E_LYA * EV)
        print(f"{T:8.0f} {b:11.4e} {s:11.4e} {f2s:7.4f} {1-f2s:7.4f} "
              f"{a2p:11.4e} {a2p_moc:11.4e} {a2p/a2p_moc:7.4f}")


if __name__ == "__main__":
    main()
