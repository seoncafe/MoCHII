#!/usr/bin/env python3
"""Validate the isotropic external radiation field (source_geometry=external).

Physics: an optically thin box/sphere (tau_ion << 1) uniformly illuminated by
an isotropic external field of band-integrated mean intensity J must tally an
interior band-integrated mean intensity equal to J.  This tests the
normalization L = pi*J*A_surface and the cosine-weighted inward sampling.

Reads '<base>_rates.h5': J_nu(bin,leaf) [erg/s/cm^2/Hz/sr], E_bin/dE_bin [eV],
LeafXYZ [code units].  Band-integrates over the ionizing bins (E >= eion_min):
    J = sum_bin J_nu(bin,leaf) * dnu(bin),   dnu[Hz] = dE_bin[eV] * eV_to_Hz.
eV_to_Hz uses the SAME cgs constants the code uses to build ion_dnu
(ev2erg / h_planck_cgs), so this integral reproduces the code's own band
convention exactly.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..",
                                "tools", "python"))
from mochii_output import read_sections

# cgs, identical to define.f90 (ev2erg / h_planck_cgs) -> Hz per eV
EV2ERG = 1.602176634e-12
H_CGS = 6.62607015e-27
EV_TO_HZ = EV2ERG / H_CGS           # = 2.417989...e14 Hz/eV


def parse_in(fname):
    """Pull par%ext_intensity, par%eion_min, par%rmax from a namelist .in."""
    vals = {}
    with open(fname) as fh:
        for ln in fh:
            s = ln.split("!", 1)[0].strip()
            if not s.lower().startswith("par%") or "=" not in s:
                continue
            key, rhs = s.split("=", 1)
            key = key.strip().lower().replace("par%", "")
            try:
                vals[key] = float(rhs.strip().split()[0])
            except ValueError:
                pass
    return vals


def band_integrated_J(rates_file, eion_min):
    sec = read_sections(rates_file)
    jnu = np.asarray(sec["J_nu"]["data"], float)
    ebin = np.asarray(sec["E_bin"]["data"], float).ravel()
    debin = np.asarray(sec["dE_bin"]["data"], float).ravel()
    xyz = np.asarray(sec["LeafXYZ"]["data"], float)
    if xyz.shape[0] != 3:
        xyz = xyz.T                         # -> (3, nleaf)
    nleaf = xyz.shape[1]
    nbin = ebin.size
    # orient J_nu to (nbin, nleaf)
    if jnu.shape == (nbin, nleaf):
        pass
    elif jnu.shape == (nleaf, nbin):
        jnu = jnu.T
    else:
        raise ValueError(f"J_nu shape {jnu.shape} matches neither "
                         f"(nbin={nbin}, nleaf={nleaf})")
    dnu = debin * EV_TO_HZ                  # Hz
    ion = ebin >= eion_min                  # ionizing-band bins only
    # J(leaf) = sum over ionizing bins of J_nu * dnu
    Jleaf = (jnu[ion, :] * dnu[ion, None]).sum(axis=0)
    return Jleaf, xyz


def report(tag, rates_file, in_file, interior_mask_fn):
    p = parse_in(in_file)
    Jext = p["ext_intensity"]
    eion_min = p.get("eion_min", 13.598)
    Jleaf, xyz = band_integrated_J(rates_file, eion_min)
    mask = interior_mask_fn(xyz, p)
    Ji = Jleaf[mask]
    Jmean = Ji.mean()
    ratio = Jmean / Jext
    print(f"\n=== {tag} ===")
    print(f"  input  J (ext_intensity) = {Jext:.6e}")
    print(f"  interior leaves sampled  = {mask.sum()} of {xyz.shape[1]}")
    print(f"  measured interior <J>    = {Jmean:.6e}")
    print(f"  ratio <J>/J_ext          = {ratio:.4f}")
    print(f"  interior J scatter       = {Ji.std()/Jmean*100:.2f}% (leaf-to-leaf)")
    # per-axis uniformity: split interior into low/high halves on each axis
    for a, nm in enumerate("xyz"):
        c = xyz[a, mask]
        lo = Ji[c < np.median(c)].mean()
        hi = Ji[c >= np.median(c)].mean()
        print(f"  {nm}: <J>(low half)={lo:.4e}  <J>(high half)={hi:.4e}  "
              f"asym={ (hi-lo)/Jmean*100:+.2f}%")
    ok = abs(ratio - 1.0) < 0.05
    print(f"  --> {'PASS' if ok else 'FAIL'} (|ratio-1| "
          f"{'<' if ok else '>='} 5%)")
    return ok


def main():
    here = os.path.dirname(__file__)
    passed = True
    # external_rec: interior = |x|,|y|,|z| < 0.5*half-extent (half-extent=xmax)
    f = os.path.join(here, "ext_rec_rates.h5")
    if os.path.exists(f):
        def rec_mask(xyz, p):
            h = 0.5 * p.get("xmax", 2.0)          # 0.5 * box half-extent
            return ((np.abs(xyz[0]) < h) & (np.abs(xyz[1]) < h)
                    & (np.abs(xyz[2]) < h))
        passed &= report("external_rec (box faces)", f,
                         os.path.join(here, "ext_rec.in"), rec_mask)
    else:
        print("ext_rec_rates.h5 not found (run ext_rec.in first)")
        passed = False
    # external_sph: interior = r < 0.6*rmax
    f = os.path.join(here, "ext_sph_rates.h5")
    if os.path.exists(f):
        def sph_mask(xyz, p):
            r = np.sqrt(xyz[0]**2 + xyz[1]**2 + xyz[2]**2)
            return r < 0.6 * p["rmax"]
        passed &= report("external_sph (bounding sphere)", f,
                         os.path.join(here, "ext_sph.in"), sph_mask)
    else:
        print("ext_sph_rates.h5 not found (run ext_sph.in first)")
        passed = False
    print(f"\nOVERALL: {'PASS' if passed else 'FAIL'}")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
