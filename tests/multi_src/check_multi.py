#!/usr/bin/env python3
"""Validate the multi-component ionizing source model.

Three gates, all on an optically thin uniform box (tau_ion << 1), where the
band-integrated mean intensity J follows an analytic superposition:

  G-two-points  two internal point sources.  Leaf J(r) vs
                J(r) = sum_i L_i / (16 pi^2 |r - r_i|^2)  (r, r_i in cm).
  G-mixed       one internal point source + isotropic external field.  Leaf J(r)
                vs J_ext + L / (16 pi^2 r^2).
  G-spectra     two point sources of equal L but different tstar.  The
                volume-summed spectrum J_nu(bin) of two_temp must equal the sum
                of the two single-source runs (two_temp_s1 + two_temp_s2).

J is read from '<base>_rates.h5' and band-integrated over the ionizing bins
exactly as tests/ext_field/check_ext.py does (same cgs eV_to_Hz).
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
PC2CM = 3.0856776e18                    # define.f90 pc2cm (distance_unit='pc')

HERE = os.path.dirname(__file__)


def load(rates_file):
    """Return (Jleaf[nleaf], Jnu[nbin,nleaf], xyz[3,nleaf] code units, ebin, debin)."""
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
    ion = ebin >= 13.598
    Jleaf = (jnu[ion, :] * dnu[ion, None]).sum(axis=0)
    return Jleaf, jnu, xyz, ebin, debin


def near_source(xyz, srcs, ncell):
    """Mask leaves within ncell cell widths of ANY source (code units)."""
    cell = 4.0 / 32.0                   # 2*xmax/nx = 0.125 pc
    m = np.zeros(xyz.shape[1], bool)
    for (sx, sy, sz) in srcs:
        r = np.sqrt((xyz[0]-sx)**2 + (xyz[1]-sy)**2 + (xyz[2]-sz)**2)
        m |= (r < ncell*cell)
    return m


def report_superpos(tag, rates_file, srcs, lums, Jext=0.0):
    """srcs: list of (x,y,z) code units; lums: ionizing L [erg/s]; Jext added."""
    Jleaf, _, xyz, _, _ = load(rates_file)
    exclude = near_source(xyz, srcs, 3)
    keep = ~exclude
    Jana = np.full(xyz.shape[1], Jext, float)
    for (sx, sy, sz), L in zip(srcs, lums):
        rcm = np.sqrt((xyz[0]-sx)**2 + (xyz[1]-sy)**2 + (xyz[2]-sz)**2) * PC2CM
        Jana += L / (16.0 * np.pi**2 * rcm**2)
    ratio = Jleaf[keep] / Jana[keep]
    med = np.median(np.abs(ratio - 1.0))
    mean_ratio = ratio.mean()
    print(f"\n=== {tag} ===")
    print(f"  leaves kept (>3 cells from sources) = {keep.sum()} of {keep.size}")
    print(f"  median |J_meas/J_analytic - 1|      = {med*100:.2f}%")
    print(f"  mean   J_meas/J_analytic            = {mean_ratio:.4f}")
    print(f"  ratio scatter (std)                 = {ratio.std():.4f}")
    ok = (med < 0.03) and (abs(mean_ratio - 1.0) < 0.02)
    print(f"  --> {'PASS' if ok else 'FAIL'} "
          f"(median<3% and |mean-1|<0.02)")
    return ok


def report_spectra(tag, two_file, s1_file, s2_file):
    _, jnu2, _, ebin, _ = load(two_file)
    _, js1, _, _, _ = load(s1_file)
    _, js2, _, _, _ = load(s2_file)
    spec2 = jnu2.sum(axis=1)             # volume-summed J_nu(bin), two-source
    specs = (js1 + js2).sum(axis=1)      # sum of the two single-source runs
    ion = ebin >= 13.598
    peak = spec2[ion].max()
    pop = ion & (spec2 > 0.02 * peak)    # well-populated bins only
    rel = np.abs(spec2[pop] - specs[pop]) / specs[pop]
    med = np.median(rel)
    mx = rel.max()
    tot_ratio = spec2[ion].sum() / specs[ion].sum()
    print(f"\n=== {tag} ===")
    print(f"  populated ionizing bins            = {pop.sum()} of {ion.sum()}")
    print(f"  median |two - (s1+s2)| / (s1+s2)    = {med*100:.2f}%")
    print(f"  max    per-bin rel diff            = {mx*100:.2f}%")
    print(f"  band-integrated ratio two/(s1+s2)  = {tot_ratio:.4f}")
    ok = (med < 0.03) and (abs(tot_ratio - 1.0) < 0.02)
    print(f"  --> {'PASS' if ok else 'FAIL'} "
          f"(median<3% and |total-1|<0.02)")
    return ok


def main():
    passed = True
    # G-two-points: src1 (-1,0,0) L=3e37, src2 (+1,+0.5,0) L=1e37
    f = os.path.join(HERE, "twop_rates.h5")
    if os.path.exists(f):
        passed &= report_superpos(
            "G-two-points", f,
            srcs=[(-1.0, 0.0, 0.0), (1.0, 0.5, 0.0)],
            lums=[3.0e37, 1.0e37])
    else:
        print("twop_rates.h5 not found"); passed = False
    # G-mixed: point origin L=1e37 + external J=1e-5
    f = os.path.join(HERE, "mixed_rates.h5")
    if os.path.exists(f):
        passed &= report_superpos(
            "G-mixed", f,
            srcs=[(0.0, 0.0, 0.0)], lums=[1.0e37], Jext=1.0e-5)
    else:
        print("mixed_rates.h5 not found"); passed = False
    # G-spectra: two_temp vs two single-source runs
    ft = os.path.join(HERE, "twotemp_rates.h5")
    f1 = os.path.join(HERE, "twotemp_s1_rates.h5")
    f2 = os.path.join(HERE, "twotemp_s2_rates.h5")
    if os.path.exists(ft) and os.path.exists(f1) and os.path.exists(f2):
        passed &= report_spectra("G-spectra", ft, f1, f2)
    else:
        print("two_temp rates files not found"); passed = False

    print(f"\nOVERALL: {'PASS' if passed else 'FAIL'}")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
