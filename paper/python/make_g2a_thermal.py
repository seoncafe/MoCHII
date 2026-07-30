#!/usr/bin/env python3
"""G2a gate: H/He nebula with thermal balance vs a 1D exact reference.

Same Stromgren setup as tests/g1_stromgren (nH = 100, Q_H = 1e49,
tstar = 4e4 K, case B), but T_e solved from heating = cooling per cell.
The 1D reference is the Gauss-Seidel radial marching of the G1 check with
the SAME thermal balance as thermal_mod/cooling_mod:
  - photoheating from the attenuated bin fluxes;
  - recombination cooling (Hui & Gnedin 1997 case B; He II as kT alpha;
    He III hydrogenic x2);
  - free-free on the code's own gbar(T, Z), Z = 1 for H II and He II and
    Z = 2 for He III;
  - collisional-ionization cooling;
  - H I line cooling from the SAME Tier-1 file data/atomic/cooling_tier1_h_1.txt.

Expected physics: without metal cooling Te ~ 13-20 kK (the known pure-H/He
overshoot, docs/PLAN.md section 8).  Gate: median |Te_MC/Te_1D - 1| < 1%
over the ionized interior (0.3 < r < 2.8 pc).

Rate coefficients come from tools/python/mochii_rates.py, so that the reference
recombines and radiates at the rates the run used; otherwise the deviation the
figure reports is partly a difference between two atomic datasets rather than a
property of the transport.  Three such differences were removed:

  - recombination, which mirrors src/recomb_mod.f90 for the default
    par%recomb_model = 'badnell_milne' (Badnell total minus the Milne
    ground-level rate); the run plotted is case B (par%case_ab = 'B').  This
    generator used to carry its own Hui & Gnedin (1997) case-B coefficients,
    which over the 19-24 kK of this test differ by 1.2% and made 0.20 of the
    0.97 percentage points of median deviation atomic data rather than
    transport;
  - the free-free Gaunt factor, which cooling_mod resolves by net charge.  The
    1.1 + 0.34 exp[-(5.5-logT)^2/3] approximation this generator used to carry
    is within 0.06-0.6% of the code at Z = 1, so its T dependence was never the
    problem; using it for He III as well was, since gbar(Z=2) is 6.5% below
    gbar(Z=1) at 10^4 K;
  - the H I and He II photoionization cross sections, which photo_xsec takes
    from the exact hydrogenic expression where this generator used the VFKY96 fit.
    The fit is 0.775% high at the H I threshold, and a hot star's rate integral
    carries most of its weight there.  He I is VFKY96 on both sides.

tests/rates/check_python_rates.py holds the module to the Fortran.  The
recombination *cooling* coefficients below stay Hui & Gnedin, as cooling_mod
has them.
"""
import sys as _sys

import numpy as np
import os as _os
# Figures are written to paper/figures/ regardless of the working directory.
FIGDIR = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "figures")

_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  "..", "..", "tools", "python"))
from mochii_rates import (alphaB_HII, alphaB_HeII, alphaB_HeIII, gbar_ff,
                          sigma_HI, sigma_HeI, sigma_HeII,
                          betaB_HII, beta_HeII, beta_HeIII,
                          ci_HI, ci_HeI, ci_HeII)

# Plot 1 in THIN Monte Carlo cells (random sample); all numbers use every cell.
THIN = 5
import h5py

EV2ERG = 1.602176634e-12
KB     = 1.380649e-16
PC2CM  = 3.0856776e18
ETH_HI, ETH_HeI, ETH_HeII = 13.598, 24.587, 54.416

NNU, EMIN, EMAX = 32, 13.598, 100.0
TSTAR, LBAND    = 4.0e4, 3.177837e38
NH, YHE         = 100.0, 0.1
RSPH_PC         = 4.0
TE_MIN, TE_MAX  = 3.0e3, 5.0e4












def load_tier1(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            rows.append(s)
    n = int(rows[0])
    A, Ti = [], []
    for k in range(1, n + 1):
        a, t = (float(v) for v in rows[k].split())
        A.append(a); Ti.append(t)
    return np.array(A), np.array(Ti)


T1_A, T1_T = load_tier1("../../data/atomic/cooling_h_1.txt")


def lam_HI(T):
    return np.sum(T1_A*np.exp(-T1_T/T))/np.sqrt(T)


def band_bins(nnu=NNU, emin=EMIN, emax=EMAX, nsub=32):
    edge = np.exp(np.linspace(np.log(emin), np.log(emax), nnu + 1))
    ec = np.sqrt(edge[:-1]*edge[1:])
    lum = np.zeros(nnu)
    for i in range(nnu):
        es = edge[i] + (edge[i+1]-edge[i])*(np.arange(nsub)+0.5)/nsub
        lum[i] = np.sum(es**3/np.expm1(es*EV2ERG/(KB*TSTAR)))*(edge[i+1]-edge[i])/nsub
    lum *= LBAND/lum.sum()
    return ec, lum


def solve_ion_scalar(gH, gHe1, gHe2, T, xHI, xHeI, xHeII):
    """Mirror of ion_balance_mod::solve_ion_cell (case B)."""
    aH, aHe2, aHe3 = alphaB_HII(T), alphaB_HeII(T), alphaB_HeIII(T)
    cH, cHe1, cHe2 = ci_HI(T), ci_HeI(T), ci_HeII(T)
    xHeIII = max(0.0, 1.0 - xHeI - xHeII)
    ne = max(NH*((1.0 - xHI) + YHE*(xHeII + 2.0*xHeIII)), 1e-12*NH)
    for _ in range(200):
        ne_old = ne
        den = gH + (cH + aH)*ne
        xHI = aH*ne/den if den > 0 else 1.0
        r1 = (gHe1 + cHe1*ne)/(aHe2*ne)
        r2 = (gHe2 + cHe2*ne)/(aHe3*ne)
        xHeI = 1.0/(1.0 + r1 + r1*r2)
        xHeII = xHeI*r1
        xHeIII = xHeII*r2
        ne = NH*((1.0 - xHI) + YHE*(xHeII + 2.0*xHeIII))
        ne = max(0.5*(ne + ne_old), 1e-12*NH)
        if abs(ne - ne_old) <= 1e-10*ne:
            break
    return xHI, xHeI, xHeII, ne


def cooling_total(T, ne, xHI, xHeI, xHeII):
    xHeIII = max(0.0, 1.0 - xHeI - xHeII)
    nHI, nHII = NH*xHI, NH*(1.0 - xHI)
    nHeI, nHeII_n, nHeIII = NH*YHE*xHeI, NH*YHE*xHeII, NH*YHE*xHeIII
    cool = ne*(nHII*betaB_HII(T) + nHeII_n*beta_HeII(T)
               + nHeIII*beta_HeIII(T))
    #--- free-free: Z = 1 for H II and He II, Z = 2 (so Z^2 = 4) for He III.
    cool += 1.42554e-27*np.sqrt(T)*ne*((nHII + nHeII_n)*gbar_ff(T, 1)
                                       + 4.0*nHeIII*gbar_ff(T, 2))
    cool += ne*(nHI*ci_HI(T)*ETH_HI + nHeI*ci_HeI(T)*ETH_HeI
                + nHeII_n*ci_HeII(T)*ETH_HeII)*EV2ERG
    cool += ne*nHI*lam_HI(T)
    return cool


def net_rate(gH, gHe1, gHe2, hH, hHe1, hHe2, T, xs):
    xHI, xHeI, xHeII, ne = solve_ion_scalar(gH, gHe1, gHe2, T, *xs)
    heat = NH*(xHI*hH + YHE*(xHeI*hHe1 + xHeII*hHe2))
    return heat - cooling_total(T, ne, xHI, xHeI, xHeII), (xHI, xHeI, xHeII), ne


def reference_1d(nr=8000, nsweep=4):
    ec, lum = band_bins()
    sH, sHe1, sHe2 = sigma_HI(ec), sigma_HeI(ec), sigma_HeII(ec)
    hnu = ec*EV2ERG
    dr_pc = RSPH_PC/nr
    dr_cm = dr_pc*PC2CM
    r = (np.arange(nr) + 0.5)*dr_pc
    r_cm = r*PC2CM
    xHI = np.full(nr, 1e-3); xHeI = np.zeros(nr); xHeII = np.ones(nr)
    Te = np.full(nr, 1.0e4)
    for sweep in range(nsweep):
        te_prev = Te.copy()
        tau = np.zeros(len(ec))
        for i in range(nr):
            xs = (xHI[i], xHeI[i], xHeII[i])
            for _ in range(2):     # self-shell attenuation consistency
                kap_i = NH*(xs[0]*sH + YHE*(xs[1]*sHe1 + xs[2]*sHe2))
                att = np.exp(-(tau + 0.5*kap_i*dr_cm))
                flux = lum*att/(4*np.pi*r_cm[i]**2)
                gH   = float((flux*sH/hnu).sum())
                gHe1 = float((flux*sHe1/hnu).sum())
                gHe2 = float((flux*sHe2/hnu).sum())
                hH   = float((flux*sH*(1.0 - ETH_HI/ec)).sum())
                hHe1 = float((flux*sHe1*(1.0 - ETH_HeI/ec)).sum())
                hHe2 = float((flux*sHe2*(1.0 - ETH_HeII/ec)).sum())
                #--- bisection on log T (mirror of thermal_mod)
                tlo, thi = TE_MIN, TE_MAX
                nlo, xs_lo, _ = net_rate(gH, gHe1, gHe2, hH, hHe1, hHe2, tlo, xs)
                nhi, xs_hi, _ = net_rate(gH, gHe1, gHe2, hH, hHe1, hHe2, thi, xs)
                if nlo <= 0.0:
                    te, xs = tlo, xs_lo
                elif nhi >= 0.0:
                    te, xs = thi, xs_hi
                else:
                    for _ in range(60):
                        tm = np.sqrt(tlo*thi)
                        nm, xs_m, _ = net_rate(gH, gHe1, gHe2, hH, hHe1, hHe2, tm, xs)
                        if nm > 0:
                            tlo = tm
                        else:
                            thi = tm
                        if thi/tlo - 1.0 < 1e-5:
                            break
                    te = np.sqrt(tlo*thi)
                    _, xs, _ = net_rate(gH, gHe1, gHe2, hH, hHe1, hHe2, te, xs)
            xHI[i], xHeI[i], xHeII[i] = xs
            Te[i] = te
            kap_i = NH*(xs[0]*sH + YHE*(xs[1]*sHe1 + xs[2]*sHe2))
            tau = tau + kap_i*dr_cm
        if np.abs(Te/te_prev - 1.0).max() < 1e-6:
            break
    return r, xHI, xHeI, xHeII, Te


# ------------------------------------------------------------------ main --
print("=== 1D thermal reference (rates from tools/python/mochii_rates.py) ===")
r1d, xHI1d, xHeI1d, xHeII1d, Te1d = reference_1d()
i_in = r1d < 2.5
print(f"1D: Te at r=0.5 pc = {Te1d[np.argmin(abs(r1d-0.5))]:.1f} K, "
      f"at 2.0 pc = {Te1d[np.argmin(abs(r1d-2.0))]:.1f} K")
print(f"1D: R(x_HI=0.5) = {r1d[np.argmin(abs(xHI1d-0.5))]:.4f} pc")

f = h5py.File("../../tests/g2a_thermal/g2a_thermal_rates.h5", "r")
xyz = f["LeafXYZ"]["data"][:]
r_mc = np.sqrt((xyz**2).sum(axis=0))
te_mc = f["T_e"]["data"][:]
xHI_mc = f["x_HI"]["data"][:]
ne_mc = f["n_e"]["data"][:]

sel = (r_mc > 0.3) & (r_mc < 2.8) & (ne_mc > 0)
te_ref = np.interp(r_mc[sel], r1d, Te1d)
dev = te_mc[sel]/te_ref - 1.0
print(f"\n=== G2a gate: Te(r), MC vs 1D over 0.3 < r < 2.8 pc ===")
print(f"median |dev| = {np.median(np.abs(dev))*100:.3f}%   "
      f"mean bias = {dev.mean()*100:+.3f}%   max |dev| = {np.abs(dev).max()*100:.3f}%")
stat = "PASS" if np.median(np.abs(dev)) < 0.01 else "FAIL"
print(f"GATE: {stat}")

# --- plot ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Paper figures: axis, tick, and title fonts sized to ~85% of the body text.
_S = 1.58
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "figure.titlesize": _b * 1.2 * _S,
})

order = np.argsort(r_mc)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
# Random 1-in-THIN sample of the Monte Carlo cells for plotting; the gate
# numbers above are computed from the full cell set.  The same sample is used
# in both panels.
rng = np.random.default_rng(20260725)
_gas = np.nonzero(ne_mc > 0)[0]
sel_gas = rng.choice(_gas, size=max(1, _gas.size // THIN), replace=False)
MC_STYLE = dict(marker="o", ms=1.8, mfc="tab:red", mec="none", ls="none",
                alpha=0.5, rasterized=True)
ax1.plot(r_mc[sel_gas], te_mc[sel_gas], label="MoCHII", **MC_STYLE)
ax1.plot(r1d, Te1d, "k-", lw=1.5, label="1D exact")
ax1.set_xlabel(r"$r$ [pc]");  ax1.set_ylabel(r"$T_e$ [K]")
ax1.set_xlim(0, 4);  ax1.legend(frameon=False, markerscale=3)
ax1.set_title(r"$T_e(r)$, H/He only (no metal cooling)")

ax2.plot(r_mc[sel_gas], xHI_mc[sel_gas], **MC_STYLE)
ax2.plot(r1d, xHI1d, "k-", lw=1.5)
ax2.set_xlabel(r"$r$ [pc]");  ax2.set_ylabel(r"$x_{\rm HI}$")
ax2.set_xlim(0, 4)
ax2.set_title(r"H ionization with thermal balance")
fig.tight_layout()
fig.savefig(_os.path.join(FIGDIR, "g2a_thermal_check.pdf"), dpi=200, bbox_inches="tight")
print("wrote g2a_thermal_check.pdf")
