#!/usr/bin/env python3
"""G0 gate: MoCHII Gamma(r) vs the analytic point-source attenuation.

Setup (g0_gamma.in): uniform sphere nH = 1 cm^-3, R = 1 pc (voxelized on a
level-5 octree, boxlen 2 pc), central point source, blackbody 1e5 K over
13.598-100 eV carrying L = 1e38 erg/s, fixed state x_HI = 0.1, x_HeI = 1,
y_He = 0.1, no dust, no scattering.

Analytic reference, evaluated with the SAME frequency bins and bin-center
cross sections as the code (this isolates the transport + tally + rate
bookkeeping; bin-count convergence is a separate check below).  The cross
sections come from tools/python/mochii_rates.py, which mirrors
src/photo_xsec.f90: H I and He II from the exact hydrogenic expression, He I
from VFKY96.  This check used to evaluate H I and He II from the VFKY96 fit
instead, which is 0.775% high at the H I threshold -- where a 1e5 K blackbody
puts most of the weight of the rate integral -- so the "analytic" reference was
not the integral the code computed.  tests/rates/check_python_rates.py holds the
module to the Fortran.

    Gamma_i(r) = sum_b L_b sigma_i(E_b) exp(-kappa(E_b) r_cm)
                 / (4 pi r_cm^2 h nu_b)
    kappa(E)   = nH [ x_HI sigma_HI + y_He x_HeI sigma_HeI ] (x_HeII = 0)

Pass criteria (r in [0.2, 0.85] pc; inside: uniform medium, away from the
center cell-size gradient and the voxelized sphere edge):
    median |Gamma_MC/Gamma_ana - 1| < 1%,  max < 3%.

Also: (a) cross-section values vs the independent EXHALE implementations
(hydrogenic H formula; VFKY96 He I) -- now a check of one hydrogenic
implementation against another, not of a fit against the exact form,
(b) the 16-bin rate integral vs a 2048-bin quasi-continuous integral
(PLAN: 10-20 bins -> percent level).
"""
import os as _os
import sys

import numpy as np
import h5py

sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                 "..", "..", "tools", "python"))
from mochii_rates import sigma_vfky96, sigma_HI, sigma_HeI, sigma_HeII

# --- constants (match define.f90) ---
EV2ERG = 1.602176634e-12
HP     = 6.62607015e-27
KB     = 1.380649e-16
PC2CM  = 3.0856776e18
ETH_HI, ETH_HeI, ETH_HeII = 13.598, 24.587, 54.416

# --- run parameters (match g0_gamma.in) ---
NNU, EMIN, EMAX = 16, 13.598, 100.0
TSTAR, LBAND    = 1.0e5, 1.0e38
NH, XHI, XHEI, YHE = 1.0, 0.1, 1.0, 0.1
RSPH_CM = 1.0 * PC2CM


def band_bins(nnu, emin, emax, nsub=32):
    """Replicate ion_band_mod: log bins, geometric centers, Planck weights
    via nsub-point midpoint integration, normalized to LBAND."""
    edge = np.exp(np.linspace(np.log(emin), np.log(emax), nnu + 1))
    ec = np.sqrt(edge[:-1]*edge[1:])
    lum = np.zeros(nnu)
    for i in range(nnu):
        es = edge[i] + (edge[i+1]-edge[i])*(np.arange(nsub)+0.5)/nsub
        x = es*EV2ERG/(KB*TSTAR)
        lum[i] = np.sum(es**3/np.expm1(x))*(edge[i+1]-edge[i])/nsub
    lum *= LBAND/lum.sum()
    return ec, lum


def gamma_ana(r_cm, ec, lum, sig):
    """Analytic bin-sum Gamma at radius r (uniform medium to RSPH_CM)."""
    kap = NH*(XHI*sigma_HI(ec) + YHE*XHEI*sigma_HeI(ec))     # cm^-1 per bin
    path = np.minimum(r_cm, RSPH_CM)
    tau = np.outer(path, kap)
    hnu = ec*EV2ERG
    return np.sum(lum*sig(ec)/hnu*np.exp(-tau), axis=1)/(4*np.pi*r_cm**2)


def gamma_ana_cellavg(xyz_pc, half_pc, ec, lum, sig, nsub=4):
    """Cell-volume average of gamma_ana on an nsub^3 sub-grid: the MC tally
    is a volume average over each leaf, so the reference must be too (the
    1/r^2 curvature alone is a ~(Delta/2r)^2 systematic at small r)."""
    off = (np.arange(nsub) + 0.5)/nsub - 0.5           # cell-centered offsets
    ox, oy, oz = np.meshgrid(off, off, off, indexing="ij")
    pts = np.stack([o.ravel() for o in (ox, oy, oz)])  # (3, nsub^3)
    nl = xyz_pc.shape[1]
    acc = np.zeros(nl)
    for k in range(pts.shape[1]):
        p = xyz_pc + 2.0*half_pc*pts[:, k:k+1]
        r = np.sqrt((p**2).sum(axis=0))*PC2CM
        acc += gamma_ana(r, ec, lum, sig)
    return acc/pts.shape[1]


# ---------------------------------------------------------------- main ---
f = h5py.File("g0_gamma_rates.h5", "r")
g_mc  = {s: f[f"Gamma_{s}"]["data"][:] for s in ("HI", "HeI", "HeII")}
xyz   = f["LeafXYZ"]["data"][:]          # (3, nleaf) as stored
e_bin = f["E_bin"]["data"][:]
r_pc  = np.sqrt((xyz**2).sum(axis=0))
r_cm  = r_pc*PC2CM

ec, lum = band_bins(NNU, EMIN, EMAX)
assert np.allclose(ec, e_bin, rtol=1e-10), "bin centers disagree with the code"

print("=== G0 gate: Gamma(r) vs analytic attenuation (cell-averaged) ===")
HALF_PC = 1.0/32.0            # level-5 half cell [pc]
sel = (r_pc > 0.15) & (r_pc < 0.85)
ok = True
for name, sig in (("HI", sigma_HI), ("HeI", sigma_HeI)):
    ana = gamma_ana_cellavg(xyz[:, sel], HALF_PC, ec, lum, sig)
    ratio = g_mc[name][sel]/ana
    med  = np.median(np.abs(ratio - 1.0))
    bias = ratio.mean() - 1.0
    mx   = np.abs(ratio - 1.0).max()
    stat = "PASS" if (med < 0.01 and mx < 0.03) else "FAIL"
    ok = ok and stat == "PASS"
    print(f"Gamma_{name:4s}: median |dev| = {med*100:6.3f}%   mean bias = "
          f"{bias*100:+6.3f}%   max |dev| = {mx*100:6.3f}%   [{stat}]")
print(f"Gamma_HeII : max = {g_mc['HeII'].max():.3e} s^-1 "
      f"(x_HeII = 0 -> rates finite but no absorber; sanity only)")

# --- cross sections vs independent EXHALE implementations ---
print("\n=== cross sections vs EXHALE (independent formulas) ===")
E = np.logspace(np.log10(13.7), np.log10(500.0), 200)


def sig_H_hydrogenic(E, Z=1.0):
    """EXHALE cross_sec.f90 sigma(E,Z): exact nonrelativistic hydrogenic."""
    E0 = 13.6*Z*Z
    eps = np.sqrt(np.maximum(E/E0 - 1.0, 1e-12))
    s = 6.3/(Z*Z)*(E0/E)**4*np.exp(4.0 - 4.0*np.arctan(eps)/eps) \
        / (1.0 - np.exp(-2.0*np.pi/eps))*1e-18
    return np.where(E >= 0.99999*E0, s, 0.0)


dev_H = np.abs(sigma_HI(E)/sig_H_hydrogenic(E) - 1.0).max()
Ehe = E[E >= 24.6]
dev_He = np.abs(sigma_HeI(Ehe)/sigma_vfky96(Ehe, 24.59, 1.361e1, 9.492e2,
                1.469, 3.188, 2.039, 4.434e-1, 2.136) - 1.0).max()
print(f"sigma_HI  vs EXHALE hydrogenic         : max dev = {dev_H*100:5.2f}% "
      f"(same formula; E_th 13.598 vs EXHALE's 13.6 eV)")
print(f"sigma_HeI vs EXHALE VFKY96 parameters  : max dev = {dev_He*100:.2e}%")
print(f"threshold values: sigma_HI(13.6+) = {sigma_HI(13.599)/1e-18:6.3f} Mb "
      f"(6.30), sigma_HeI(24.6+) = {sigma_HeI(24.6)/1e-18:6.3f} Mb (~7.4), "
      f"sigma_HeII(54.4+) = {sigma_HeII(54.42)/1e-18:6.3f} Mb (~1.58)")

# --- bin-count convergence of the rate integral ---
print("\n=== bin-count convergence (r = 0.5 pc) ===")
r0 = np.array([0.5*PC2CM])
g_ref = None
for nnu in (2048, 64, 32, 16, 8):
    ecb, lumb = band_bins(nnu, EMIN, EMAX)
    kap = NH*(XHI*sigma_HI(ecb) + YHE*XHEI*sigma_HeI(ecb))
    g = np.sum(lumb*sigma_HI(ecb)/(ecb*EV2ERG)*np.exp(-kap*r0[0]))/(4*np.pi*r0[0]**2)
    if g_ref is None:
        g_ref = g
        print(f"  nnu = {nnu:5d}: Gamma_HI = {g:.6e} s^-1 (reference)")
    else:
        print(f"  nnu = {nnu:5d}: Gamma_HI = {g:.6e} s^-1  dev = "
              f"{(g/g_ref-1)*100:+6.3f}%")

# --- figure ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

order = np.argsort(r_pc)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True,
                               gridspec_kw=dict(height_ratios=[3, 1]))
for name, sig, c in (("HI", sigma_HI, "C0"), ("HeI", sigma_HeI, "C1")):
    ana_c = gamma_ana_cellavg(xyz, HALF_PC, ec, lum, sig)
    ax1.plot(r_pc[order], g_mc[name][order], ".", ms=1.5, color=c, alpha=0.3)
    ax1.plot(r_pc[order], ana_c[order], "k-", lw=0.8)
    ax2.plot(r_pc[order], g_mc[name][order]/ana_c[order], ".", ms=1.5,
             color=c, alpha=0.3, label=name)
ax1.set_yscale("log")
ax1.set_ylabel(r"$\Gamma$ [s$^{-1}$]")
ax1.set_title(r"G0 gate: $\Gamma(r)$, MC (points) vs analytic (line)")
ax2.axhline(1.0, color="k", lw=0.8)
ax2.axvspan(0.0, 0.2, color="0.9")
ax2.axvspan(0.85, 1.05, color="0.9")
ax2.set_ylim(0.9, 1.1)
ax2.set_xlim(0.0, 1.05)
ax2.set_xlabel(r"$r$ [pc]")
ax2.set_ylabel("MC / analytic")
ax2.legend(frameon=False, markerscale=8)
fig.tight_layout()
fig.savefig("g0_gamma_check.png", dpi=140)
print("\nwrote g0_gamma_check.png")
print("GATE:", "PASS" if ok else "FAIL")
#--- The gate must fail the PROCESS, not just print: a CI runner (and a
#--- human skimming a log) otherwise reads a FAIL as a successful run.
sys.exit(0 if ok else 1)

