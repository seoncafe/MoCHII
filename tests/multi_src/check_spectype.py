#!/usr/bin/env python3
"""Validate par%spectrum_type (physical spectrum units), the derive/rescale
normalization semantics, absolute external fields, and the ISRF presets.

All runs are on the optically thin uniform box (car 32^3, xmax=2 pc,
nH_const=1e-4, gas_niter=0), band-integrating J_nu with the code's own cgs
eV->Hz factor.  Gates:

  G-per_ev-equiv     spectrum_type='per_ev' (absolute L_E), NO src_lum -> src_lum
                     DERIVED ~1e37; volume-summed J_nu(bin) reproduces
                     twotemp_file (shape, src_lum=1e37) per populated bin.
  G-lambda-conv      the SAME spectra as L_lambda vs Angstrom / micron match
                     the 'per_ev' run to <0.5% per bin (grid difference only).
  G-rescale          'per_ev' file WITH explicit src_lum = 3e37/1e37 reproduces the
                     'shape' twin with the same values; the rescale NOTE is
                     logged.
  G-ext-je           external absolute J_E file (band int 1e-5): interior <J> ~
                     0.998e-5; adding ext_intensity=2e-5 doubles the interior J.
  G-presets          draine/habing/mathis: interior u_FUV = (4pi/c)<J_FUV>
                     matches the analytic integral; ionizing bins carry ZERO;
                     habing = draine/1.71; mathis within a factor ~2 of draine;
                     the no-add_fuv run aborts.
"""
import os
import re
import sys
import subprocess
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..",
                                "tools", "python"))
from mochii_output import read_sections

# cgs, identical to define.f90
EV2ERG   = 1.602176634e-12
H_CGS    = 6.62607015e-27
CLIGHT   = 2.99792458e10
EV_TO_HZ = EV2ERG / H_CGS
HC_EVANG = 12398.42
EION_MIN = 13.598
FOURPI   = 4.0 * np.pi

HERE = os.path.dirname(__file__)


# ---------------------------------------------------------------- io helpers
def load(rates_file):
    """Return jnu[nbin,nleaf], ebin, debin, xyz[3,nleaf]."""
    sec = read_sections(rates_file)
    jnu = np.asarray(sec["J_nu"]["data"], float)
    ebin = np.asarray(sec["E_bin"]["data"], float).ravel()
    debin = np.asarray(sec["dE_bin"]["data"], float).ravel()
    xyz = np.asarray(sec["LeafXYZ"]["data"], float)
    if xyz.shape[0] != 3:
        xyz = xyz.T
    nleaf, nbin = xyz.shape[1], ebin.size
    if jnu.shape == (nleaf, nbin):
        jnu = jnu.T
    elif jnu.shape != (nbin, nleaf):
        raise ValueError(f"J_nu shape {jnu.shape}")
    return jnu, ebin, debin, xyz


def J_band(jnu, ebin, debin, lo, hi=None):
    """Per-leaf band-integrated J over lo<=E(<hi) [erg/s/cm^2/sr]."""
    dnu = debin * EV_TO_HZ
    m = ebin >= lo
    if hi is not None:
        m &= ebin < hi
    return (jnu[m, :] * dnu[m, None]).sum(axis=0)


def interior_mask(xyz, half=1.0):
    return ((np.abs(xyz[0]) < half) & (np.abs(xyz[1]) < half)
            & (np.abs(xyz[2]) < half))


def vol_spectrum(jnu):
    return jnu.sum(axis=1)          # volume-summed J_nu(bin)


def compare_bins(spec_t, spec_r, ebin, lo=EION_MIN):
    ion = ebin >= lo
    peak = spec_r[ion].max()
    pop = ion & (spec_r > 0.02 * peak)
    rel = np.abs(spec_t[pop] - spec_r[pop]) / spec_r[pop]
    tot = spec_t[ion].sum() / spec_r[ion].sum()
    return np.median(rel), rel.max(), tot, pop.sum(), ion.sum()


# ---------------------------------------------------------------- analytic ISRF
def draine_je(E):
    F = 1.658e6*E**2 - 2.152e5*E**3 + 6.919e3*E**4    # F = E*N(E), energy-weighted
    return np.where((E >= 5.0) & (E <= 13.6), EV2ERG*F, 0.0)


def bbody_si(T, lam_um):
    c, h, kB = 2.99792458e8, 6.62606957e-34, 1.3806488e-23
    hc2, hckB = 2.0*h*c**2, h*c/kB
    lam_m = lam_um*1.0e-6
    x = hckB/(T*lam_m)
    return np.where(x >= 700.0, 0.0, hc2/lam_m**5/np.expm1(np.minimum(x, 700.0)))


def mathis_jl_si(lam_um):
    Jl = np.zeros_like(lam_um)
    m1 = (lam_um >= 0.0912) & (lam_um < 0.110); Jl[m1] = 3069.0*lam_um[m1]**3.4172
    m2 = (lam_um >= 0.110) & (lam_um < 0.134);  Jl[m2] = 1.627
    m3 = (lam_um >= 0.134) & (lam_um < 0.250);  Jl[m3] = 0.0566*lam_um[m3]**(-1.6678)
    m4 = lam_um >= 0.250
    Jl[m4] = (1.0e-14*bbody_si(7500.0, lam_um[m4])
              + 1.65e-13*bbody_si(4000.0, lam_um[m4])
              + 4.0e-13*bbody_si(3000.0, lam_um[m4]))
    return Jl


def mathis_je(E):
    lam_um = HC_EVANG*1.0e-4/E
    return mathis_jl_si(lam_um) * (HC_EVANG*1.0e-7)/E**2


def preset_je(E, name):
    if name == "draine":  return draine_je(E)
    if name == "habing":  return draine_je(E)/1.71
    if name == "mathis":  return mathis_je(E)
    raise ValueError(name)


def u_fuv_analytic(name, elo=6.0, ehi=EION_MIN):
    E = np.linspace(elo, ehi, 40001)
    return FOURPI/CLIGHT * np.trapz(preset_je(E, name), E)


# ---------------------------------------------------------------- log parsing
def grep(logfile, pattern):
    if not os.path.exists(logfile):
        return []
    with open(logfile) as f:
        return [ln.rstrip() for ln in f if re.search(pattern, ln)]


# ---------------------------------------------------------------- gates
def g_le_equivalence():
    print("\n=== G-per_ev-equiv (spectrum_type='per_ev' derive vs shape) ===")
    fle = os.path.join(HERE, "st_per_ev_rates.h5")
    ftt = os.path.join(HERE, "twotemp_file_rates.h5")
    if not (os.path.exists(fle) and os.path.exists(ftt)):
        print("  rates missing"); return False
    jle, eb, _, _ = load(fle)
    jtt, _, _, _ = load(ftt)
    med, mx, tot, npop, nion = compare_bins(vol_spectrum(jle), vol_spectrum(jtt), eb)
    ratio = vol_spectrum(jle)[eb >= EION_MIN] / vol_spectrum(jtt)[eb >= EION_MIN]
    ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
    shape_scatter = ratio.std() / ratio.mean()
    print(f"  populated ionizing bins       = {npop} of {nion}")
    print(f"  median |le - shape| / shape   = {med*100:.3f}%")
    print(f"  max    per-bin rel diff       = {mx*100:.3f}%")
    print(f"  band-integrated ratio pe/tt   = {tot:.5f}")
    print(f"  per-bin ratio scatter (shape) = {shape_scatter*100:.4f}%")
    dl = grep(os.path.join(HERE, "st_per_ev.log"), r"L_ion \(derived\)")
    for ln in dl:
        print("  log:", ln.strip())
    derived = [float(x) for ln in dl for x in re.findall(r"[-+]?\d\.\d+E[-+]\d+", ln)]
    dok = len(derived) == 2 and all(abs(d/1.0e37 - 1.0) < 0.02 for d in derived)
    ok = (mx < 0.5e-2) and (abs(tot-1.0) < 5e-3) and (shape_scatter < 1e-3) and dok
    print(f"  --> {'PASS' if ok else 'FAIL'} "
          f"(per-bin<0.5%, |band-1|<0.5%, shape-scatter<0.1%, derived~1e37)")
    return ok


def g_lambda_conv(tag, testfile):
    print(f"\n=== G-lambda-conv ({tag} vs per_ev) ===")
    ft = os.path.join(HERE, testfile)
    fle = os.path.join(HERE, "st_per_ev_rates.h5")
    if not (os.path.exists(ft) and os.path.exists(fle)):
        print("  rates missing"); return False
    jt, eb, _, _ = load(ft)
    jle, _, _, _ = load(fle)
    med, mx, tot, npop, nion = compare_bins(vol_spectrum(jt), vol_spectrum(jle), eb)
    print(f"  populated ionizing bins       = {npop} of {nion}")
    print(f"  median |{tag} - per_ev| / per_ev = {med*100:.3f}%")
    print(f"  max    per-bin rel diff       = {mx*100:.3f}%")
    print(f"  band-integrated ratio         = {tot:.5f}")
    ok = (mx < 0.5e-2) and (abs(tot-1.0) < 5e-3)
    print(f"  --> {'PASS' if ok else 'FAIL'} (max per-bin < 0.5%, |band-1| < 0.5%)")
    return ok


def g_rescale():
    print("\n=== G-rescale ('per_ev'+explicit src_lum vs 'shape'+same) ===")
    fl = os.path.join(HERE, "st_per_ev_resc_rates.h5")
    fs = os.path.join(HERE, "st_sresc_rates.h5")
    if not (os.path.exists(fl) and os.path.exists(fs)):
        print("  rates missing"); return False
    jl, eb, _, _ = load(fl)
    js, _, _, _ = load(fs)
    med, mx, tot, npop, nion = compare_bins(vol_spectrum(jl), vol_spectrum(js), eb)
    print(f"  populated ionizing bins       = {npop} of {nion}")
    print(f"  median |per_ev - shape| / shape = {med*100:.4f}%")
    print(f"  max    per-bin rel diff       = {mx*100:.4f}%")
    print(f"  band-integrated ratio         = {tot:.6f}")
    note = grep(os.path.join(HERE, "st_per_ev_rescale.log"),
                r"rescaled from file integral")
    for ln in note:
        print("  log:", ln.strip())
    nok = len(note) >= 2
    ok = (mx < 1e-3) and (abs(tot-1.0) < 5e-4) and nok
    print(f"  --> {'PASS' if ok else 'FAIL'} "
          f"(per-bin < 0.1%, |band-1| < 0.05%, rescale NOTE x2 logged)")
    return ok


def g_ext_je():
    print("\n=== G-ext-je (external absolute J_E) ===")
    f1 = os.path.join(HERE, "st_extje_rates.h5")
    f2 = os.path.join(HERE, "st_extje2_rates.h5")
    if not (os.path.exists(f1) and os.path.exists(f2)):
        print("  rates missing"); return False
    j1, eb, db, xyz = load(f1)
    m = interior_mask(xyz, 1.0)
    Jm1 = J_band(j1, eb, db, EION_MIN)[m].mean()
    j2, eb2, db2, xyz2 = load(f2)
    m2 = interior_mask(xyz2, 1.0)
    Jm2 = J_band(j2, eb2, db2, EION_MIN)[m2].mean()
    print(f"  st_extje  interior <J>        = {Jm1:.6e}   <J>/1e-5 = {Jm1/1e-5:.4f}")
    print(f"  st_extje2 interior <J>        = {Jm2:.6e}   <J>/2e-5 = {Jm2/2e-5:.4f}")
    print(f"  doubling ratio <J2>/<J1>      = {Jm2/Jm1:.4f}")
    ok = (abs(Jm1/1e-5 - 1.0) < 0.01) and (abs(Jm2/Jm1 - 2.0) < 0.02)
    print(f"  --> {'PASS' if ok else 'FAIL'} "
          f"(<J1>/1e-5 within 1%, <J2>/<J1> within 1% of 2)")
    return ok


def g_presets():
    print("\n=== G-presets (draine / habing / mathis) ===")
    u = {}
    ok = True
    for name in ("draine", "habing", "mathis"):
        f = os.path.join(HERE, f"st_{name}_rates.h5")
        if not os.path.exists(f):
            print(f"  {name}: rates missing"); ok = False; continue
        jn, eb, db, xyz = load(f)
        mask = interior_mask(xyz, 1.0)
        Jfuv = J_band(jn, eb, db, 0.0, EION_MIN)[mask].mean()      # FUV interior J
        Jion = J_band(jn, eb, db, EION_MIN)[mask].mean()           # ionizing interior J
        um = FOURPI/CLIGHT * Jfuv
        ua = u_fuv_analytic(name)
        u[name] = um
        ionfrac = Jion / max(Jfuv, 1e-99)
        print(f"  {name:6s}: u_FUV(meas) = {um:.4e}   u_FUV(analytic) = {ua:.4e}"
              f"   meas/ana = {um/ua:.4f}   ion/FUV = {ionfrac:.2e}")
        gate = (abs(um/ua - 1.0) < 0.03) and (ionfrac < 1e-3)
        ok &= gate
    if "draine" in u and "habing" in u:
        r = u["habing"]/u["draine"]
        print(f"  habing/draine = {r:.4f}  (expect 1/1.71 = {1/1.71:.4f})")
        ok &= abs(r - 1/1.71) < 0.01
    if "draine" in u and "mathis" in u:
        r = u["mathis"]/u["draine"]
        print(f"  mathis/draine = {r:.4f}  (sanity: within a factor ~2)")
        ok &= (0.4 < r < 2.5)
    print(f"  u_FUV summary [erg/cm^3]: "
          + "  ".join(f"{k}={v:.3e}" for k, v in u.items()))
    # no-add_fuv abort
    logf = os.path.join(HERE, "st_draine_noadd.log")
    msg = grep(logf, r"FUV-only|add_fuv")
    aborted = len(msg) > 0
    print(f"  no-add_fuv abort message present: {aborted}")
    for ln in msg[:2]:
        print("  log:", ln.strip())
    ok &= aborted
    print(f"  --> {'PASS' if ok else 'FAIL'} "
          f"(u_FUV within 3%, ionizing~0, habing=draine/1.71, mathis~O(draine), abort)")
    return ok


def main():
    passed = True
    passed &= g_le_equivalence()
    passed &= g_lambda_conv("per_ang", "st_per_ang_rates.h5")
    passed &= g_lambda_conv("per_um", "st_per_um_rates.h5")
    passed &= g_rescale()
    passed &= g_ext_je()
    passed &= g_presets()
    print(f"\nOVERALL: {'PASS' if passed else 'FAIL'}")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
