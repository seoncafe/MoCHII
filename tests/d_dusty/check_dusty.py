#!/usr/bin/env python3
"""Dusty Stromgren gate (PLAN section 7 item 2).

Same G1 configuration (nH = 100, Q_H = 1e49, case B, Te = 1e4 K) plus MW
dust (WD01/D03 R_V=3.1 EUV extinction, absorption only,
C_abs(13.6 eV)/H = 1.83e-21): dust competes for Lyman-continuum photons
and shrinks the Stromgren sphere.  The exact reference is the 1D
Gauss-Seidel marching with the identical per-bin dust absorption; the
classic sharp-edge reduction (photon balance
y^3 = 1 - int_0^y 3 y'^2 tau_d e^{-...}) is quoted for context via the
implicit solution of Petrosian-type photon accounting.

Cases: full dust (global_dgr) — gate |R_eff(MC)/R_eff(1D) - 1| < 1%;
laursen09_live (dust ~ n_HI + 0.01 n_HII) — must land near the
dust-free radius (interior dust vanishes with the computed x_HII).
"""
import numpy as np
import h5py
import os

EV2ERG = 1.602176634e-12
KB     = 1.380649e-16
PC2CM  = 3.0856776e18
ETH_HI, ETH_HeI, ETH_HeII = 13.598, 24.587, 54.416
NNU, EMIN, EMAX = 32, 13.598, 100.0
TSTAR, LBAND    = 4.0e4, 3.177837e38
NH, YHE, TE     = 100.0, 0.1, 1.0e4
RSPH_PC         = 4.0
CEXT_REF        = 4.868e-22       # D03 C_ext/H at 0.55 um (par%cext_dust)
KEXT = "../../data/kext_albedo_WD_MW_3.1_60_D03.all_2009"
R_EFF_FREE_1D = 3.0365            # G1 1D reference, no dust


def sigma_vfky96(E, Eth, E0, s0, ya, P, yw, y0, y1):
    E = np.asarray(E, dtype=float)
    x = E/E0 - y0
    z = np.sqrt(x*x + y1*y1)
    Q = 5.5 - 0.5*P
    s = s0*((x-1.0)**2 + yw**2)*z**(-Q)*(1.0 + np.sqrt(z/ya))**(-P)*1e-18
    return np.where(E >= Eth, s, 0.0)


def sig_HI(E):   return sigma_vfky96(E, ETH_HI, 4.298e-1, 5.475e4, 3.288e1, 2.963, 0, 0, 0)
def sig_HeI(E):  return sigma_vfky96(E, ETH_HeI, 1.361e1, 9.492e2, 1.469, 3.188, 2.039, 4.434e-1, 2.136)
def sig_HeII(E): return sigma_vfky96(E, ETH_HeII, 1.720, 1.369e4, 3.288e1, 2.963, 0, 0, 0)


def alphaB_HII(T):
    lam = 2.0*157807.0/T
    return 2.753e-14*lam**1.500/(1.0 + (lam/2.740)**0.407)**2.242


def alphaB_HeII(T):  return 1.260e-14*(2.0*285335.0/T)**0.750


def alphaB_HeIII(T):
    lam = 2.0*631515.0/T
    return 2.0*2.753e-14*lam**1.500/(1.0 + (lam/2.740)**0.407)**2.242


def voronov(T, dE, P, A, X, K):
    U = dE*EV2ERG/(KB*T)
    return A*(1.0 + P*np.sqrt(U))*U**K*np.exp(-U)/(X + U)


def ci_HI(T):   return voronov(T, 13.6, 0.0, 2.91e-8, 0.232, 0.39)
def ci_HeI(T):  return voronov(T, 24.6, 0.0, 1.75e-8, 0.180, 0.35)
def ci_HeII(T): return voronov(T, 54.4, 1.0, 2.05e-9, 0.265, 0.25)


def band_bins():
    edge = np.exp(np.linspace(np.log(EMIN), np.log(EMAX), NNU + 1))
    ec = np.sqrt(edge[:-1]*edge[1:])
    lum = np.zeros(NNU)
    for i in range(NNU):
        es = edge[i] + (edge[i+1]-edge[i])*(np.arange(32)+0.5)/32
        lum[i] = np.sum(es**3/np.expm1(es*EV2ERG/(KB*TSTAR)))*(edge[i+1]-edge[i])/32
    lum *= LBAND/lum.sum()
    return ec, lum


def dust_sabs(ec):
    """(1-albedo) C_ext(E)/C_ext(0.55um) from the same D03 table, same
    interpolation as gas_opacity_mod (log-log C_ext, linear albedo)."""
    rows = []
    for ln in open(KEXT):
        t = ln.split()
        if len(t) >= 6:
            try:
                rows.append((float(t[0]), float(t[1]), float(t[3])))
            except ValueError:
                pass
    a = np.array(rows)
    a = a[np.argsort(a[:, 0])]
    lam, alb, cext = a.T

    def at(l):
        i = np.searchsorted(lam, l)
        w = (np.log(l) - np.log(lam[i-1]))/(np.log(lam[i]) - np.log(lam[i-1]))
        ce = np.exp(np.log(cext[i-1])*(1-w) + np.log(cext[i])*w)
        ab = alb[i-1]*(1-w) + alb[i]*w
        return (1.0 - ab)*ce
    cref_i = np.searchsorted(lam, 0.55)
    w = (np.log(0.55)-np.log(lam[cref_i-1]))/(np.log(lam[cref_i])-np.log(lam[cref_i-1]))
    cref = np.exp(np.log(cext[cref_i-1])*(1-w) + np.log(cext[cref_i])*w)
    return np.array([at(1.23984/e) for e in ec])/cref


def solve_ion(gH, gHe1, gHe2, xs):
    aH, aHe2, aHe3 = alphaB_HII(TE), alphaB_HeII(TE), alphaB_HeIII(TE)
    cH, cHe1, cHe2 = ci_HI(TE), ci_HeI(TE), ci_HeII(TE)
    xHI, xHeI, xHeII = xs
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
        ne = NH*((1.0 - xHI) + YHE*(xHeII + 2.0*xHeII*r2))
        ne = max(0.5*(ne + ne_old), 1e-12*NH)
        if abs(ne - ne_old) <= 1e-12*ne:
            break
    return xHI, xHeI, xHeII


def reference_1d(dust="full", nr=8000, nsweep=3):
    ec, lum = band_bins()
    sH, sHe1, sHe2 = sig_HI(ec), sig_HeI(ec), sig_HeII(ec)
    sd = dust_sabs(ec)*CEXT_REF          # C_abs,dust(E)/H [cm^2]
    hnu = ec*EV2ERG
    dr_cm = (RSPH_PC/nr)*PC2CM
    r = (np.arange(nr) + 0.5)*(RSPH_PC/nr)
    r_cm = r*PC2CM
    xHI = np.full(nr, 1e-3); xHeI = np.zeros(nr); xHeII = np.ones(nr)
    for sweep in range(nsweep):
        x_prev = xHI.copy()
        tau = np.zeros(NNU)
        for i in range(nr):
            xs = (xHI[i], xHeI[i], xHeII[i])
            for _ in range(2):
                if dust == "full":
                    fd = 1.0
                elif dust == "l09":
                    fd = xs[0] + 0.01*(1.0 - xs[0])
                else:
                    fd = 0.0
                kap_i = NH*(xs[0]*sH + YHE*(xs[1]*sHe1 + xs[2]*sHe2)) + NH*fd*sd
                att = np.exp(-(tau + 0.5*kap_i*dr_cm))
                flux = lum*att/(4*np.pi*r_cm[i]**2)
                gH = float((flux*sH/hnu).sum())
                gHe1 = float((flux*sHe1/hnu).sum())
                gHe2 = float((flux*sHe2/hnu).sum())
                xs = solve_ion(gH, gHe1, gHe2, xs)
            xHI[i], xHeI[i], xHeII[i] = xs
            if dust == "full":
                fd = 1.0
            elif dust == "l09":
                fd = xs[0] + 0.01*(1.0 - xs[0])
            else:
                fd = 0.0
            kap_i = NH*(xs[0]*sH + YHE*(xs[1]*sHe1 + xs[2]*sHe2)) + NH*fd*sd
            tau = tau + kap_i*dr_cm
        if np.abs(xHI - x_prev).max() < 1e-8:
            break
    vion = np.sum((1.0 - xHI)*4*np.pi*(r*PC2CM)**2)*(r[1]-r[0])*PC2CM
    return r, xHI, (3.0*vion/(4*np.pi))**(1/3)/PC2CM


def r_eff_mc(fname):
    f = h5py.File(fname, "r")
    xyz = f["LeafXYZ"]["data"][:]
    xhi = f["x_HI"]["data"][:]
    ne = f["n_e"]["data"][:]
    half = 4.0/64
    vion = np.sum((1.0 - xhi[ne > 0])*(2*half)**3)
    return (3.0*vion/(4.0*np.pi))**(1.0/3.0)


print("=== 1D references (same bins, rates, and D03 dust absorption) ===")
r1, x1, Re_full = reference_1d("full")
r2, x2, Re_l09 = reference_1d("l09")
print(f"dust-free : R_eff = {R_EFF_FREE_1D:.4f} pc (G1)")
print(f"full dust : R_eff = {Re_full:.4f} pc  (reduction {Re_full/R_EFF_FREE_1D:.3f})")
print(f"l09-live  : R_eff = {Re_l09:.4f} pc  (reduction {Re_l09/R_EFF_FREE_1D:.3f})")

print("\n=== MoCHII runs ===")
for tag, fn, ref in (("full dust", "d_dust_full_rates.h5", Re_full),
                     ("l09-live", "d_dust_l09_rates.h5", Re_l09)):
    if not os.path.exists(fn):
        print(f"[{tag}] {fn} not found; skipped")
        continue
    Re = r_eff_mc(fn)
    print(f"[{tag:10s}] R_eff = {Re:.4f} pc   vs 1D {Re/ref-1:+.2%}   "
          f"{'PASS' if abs(Re/ref-1) < 0.01 else 'FAIL'}")
