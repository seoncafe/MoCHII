#!/usr/bin/env python3
"""Validate the file-driven source spectra (par%src_spectrum_file + the
external-field fast-path spectrum fix).

Two gates on the optically thin uniform box (band-integration convention
identical to check_multi.py / check_ext.py):

  G-file-equivalence   par%src_spectrum_file (two_temp_file) reproduces the
                       per-source Planck run (two_temp).  The two internal
                       sources take their 3e4/5e4 K shapes from columns of one
                       multi-column file instead of src_tstar; the RNG stream is
                       unchanged, so the volume-summed J_nu(bin) must match
                       twotemp_rates.h5 to <1% (populated bins) and the
                       band-integrated total to <0.5%.

  G-ext-spectrum       external-only fast path honors ext_spectrum AND ext_tstar
                       (the pre-fix path ignored both, always using the global
                       spectrum).  ext_planckfile (ext_spectrum = 4e4 K file,
                       decoy par%tstar = 2e4 K) must match ext_tstar (ext_tstar =
                       4e4 K, same decoy) to <1% per bin and <0.5% band-total,
                       and both must tally interior <J>/J_ext ~ 1.00.  A match to
                       the 4e4 K twin (not the 2e4 K decoy) proves the fix.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..",
                                "tools", "python"))
from mochii_output import read_sections

# cgs, identical to define.f90
EV2ERG = 1.602176634e-12
H_CGS = 6.62607015e-27
EV_TO_HZ = EV2ERG / H_CGS

HERE = os.path.dirname(__file__)
EION_MIN = 13.598


def load(rates_file):
    """Return (Jleaf[nleaf], Jnu[nbin,nleaf], xyz[3,nleaf], ebin, debin)."""
    sec = read_sections(rates_file)
    jnu = np.asarray(sec["J_nu"]["data"], float)
    ebin = np.asarray(sec["E_bin"]["data"], float).ravel()
    debin = np.asarray(sec["dE_bin"]["data"], float).ravel()
    xyz = np.asarray(sec["LeafXYZ"]["data"], float)
    if xyz.shape[0] != 3:
        xyz = xyz.T
    nleaf = xyz.shape[1]
    nbin = ebin.size
    if jnu.shape == (nleaf, nbin):
        jnu = jnu.T
    elif jnu.shape != (nbin, nleaf):
        raise ValueError(f"J_nu shape {jnu.shape} matches neither orientation")
    dnu = debin * EV_TO_HZ
    ion = ebin >= EION_MIN
    Jleaf = (jnu[ion, :] * dnu[ion, None]).sum(axis=0)
    return Jleaf, jnu, xyz, ebin, debin


def compare_spectra(tag, test_file, ref_file, ref_label):
    """Per-bin and band-total comparison of two volume-summed spectra."""
    _, jt, _, ebin, _ = load(test_file)
    _, jr, _, _, _ = load(ref_file)
    spec_t = jt.sum(axis=1)              # volume-summed J_nu(bin)
    spec_r = jr.sum(axis=1)
    ion = ebin >= EION_MIN
    peak = spec_r[ion].max()
    pop = ion & (spec_r > 0.02 * peak)   # well-populated bins only
    rel = np.abs(spec_t[pop] - spec_r[pop]) / spec_r[pop]
    med, mx = np.median(rel), rel.max()
    tot_ratio = spec_t[ion].sum() / spec_r[ion].sum()
    print(f"\n=== {tag} ===")
    print(f"  reference                           = {ref_label}")
    print(f"  populated ionizing bins             = {pop.sum()} of {ion.sum()}")
    print(f"  median |test - ref| / ref           = {med*100:.3f}%")
    print(f"  max    per-bin rel diff             = {mx*100:.3f}%")
    print(f"  band-integrated ratio test/ref      = {tot_ratio:.5f}")
    ok = (mx < 0.01) and (abs(tot_ratio - 1.0) < 0.005)
    print(f"  --> {'PASS' if ok else 'FAIL'} "
          f"(max per-bin < 1% and |band ratio - 1| < 0.5%)")
    return ok


def interior_J(rates_file, jext, half):
    Jleaf, _, xyz, _, _ = load(rates_file)
    m = ((np.abs(xyz[0]) < half) & (np.abs(xyz[1]) < half)
         & (np.abs(xyz[2]) < half))
    Jmean = Jleaf[m].mean()
    return Jmean, Jmean / jext, m.sum()


def main():
    passed = True

    # G-file-equivalence
    ft = os.path.join(HERE, "twotemp_file_rates.h5")
    fr = os.path.join(HERE, "twotemp_rates.h5")
    if os.path.exists(ft) and os.path.exists(fr):
        passed &= compare_spectra("G-file-equivalence", ft, fr,
                                  "twotemp_rates.h5 (per-source src_tstar)")
    else:
        print("G-file-equivalence: rates files missing "
              "(run two_temp.in and two_temp_file.in)")
        passed = False

    # G-ext-spectrum
    ef = os.path.join(HERE, "..", "ext_field", "ext_planckfile_rates.h5")
    et = os.path.join(HERE, "..", "ext_field", "ext_tstar_rates.h5")
    if os.path.exists(ef) and os.path.exists(et):
        ok = compare_spectra("G-ext-spectrum", ef, et,
                             "ext_tstar_rates.h5 (ext_tstar = 4e4 K)")
        # interior <J>/J_ext ~ 1.00 for both runs (J_ext = 1e-5, half = 0.5*xmax)
        jext, half = 1.0e-5, 1.0
        for nm, f in (("ext_planckfile (ext_spectrum file)", ef),
                      ("ext_tstar     (ext_tstar Planck)", et)):
            Jmean, ratio, n = interior_J(f, jext, half)
            print(f"  interior {nm}: <J> = {Jmean:.6e}  "
                  f"<J>/J_ext = {ratio:.4f}  ({n} leaves)")
            ok &= abs(ratio - 1.0) < 0.05
        print(f"  --> {'PASS' if ok else 'FAIL'} "
              f"(spectra match AND both interior <J>/J_ext within 5%)")
        passed &= ok
    else:
        print("G-ext-spectrum: rates files missing "
              "(run ext_planckfile.in and ext_tstar.in)")
        passed = False

    print(f"\nOVERALL: {'PASS' if passed else 'FAIL'}")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
