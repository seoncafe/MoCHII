#!/usr/bin/env python3
"""G2a gate: H/He nebula with thermal balance vs a 1D exact reference.

Same Stromgren setup as tests/g1_stromgren (nH = 100, Q_H = 1e49,
tstar = 4e4 K, case B), but T_e solved from heating = cooling per cell.
The 1D reference is the Gauss-Seidel radial marching of the G1 check with
the SAME thermal balance as thermal_mod/cooling_mod:
  - photoheating from the attenuated bin fluxes;
  - recombination cooling (Hui & Gnedin 1997 case B; He II as kT alpha;
    He III hydrogenic x2);
  - free-free with gbar approximated by 1.1 + 0.34 exp(-(5.5-logT)^2/3)
    (the Fortran integrates the ported Hummer getGauntFF; the difference
    is ~1-2% of a ~10% cooling term — noted in the report);
  - collisional-ionization cooling;
  - H I line cooling from the SAME Tier-1 file data/atomic/cooling_tier1_h_1.txt.

Expected physics: without metal cooling Te ~ 13-20 kK (the known pure-H/He
overshoot, docs/PLAN.md section 8).  Gate: median |Te_MC/Te_1D - 1| < 1%
over the ionized interior (0.3 < r < 2.8 pc).
"""
import numpy as np
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


def alphaA_HII(T):
    lam = 2.0*157807.0/T
    return 1.269e-13*lam**1.503/(1.0 + (lam/0.522)**0.470)**1.923


def alphaB_HII(T):
    lam = 2.0*157807.0/T
    return 2.753e-14*lam**1.500/(1.0 + (lam/2.740)**0.407)**2.242


def alphaA_HeII(T):  return 3.000e-14*(2.0*285335.0/T)**0.654
def alphaB_HeII(T):  return 1.260e-14*(2.0*285335.0/T)**0.750


def alphaB_HeIII(T):
    lam = 2.0*631515.0/T
    return 2.0*2.753e-14*lam**1.500/(1.0 + (lam/2.740)**0.407)**2.242


def betaB_HII(T):
    lam = 2.0*157807.0/T
    return 3.435e-30*T*lam**1.970/(1.0 + (lam/2.250)**0.376)**3.720


def beta_HeII_B(T):  return KB*T*alphaB_HeII(T)


def beta_HeIII_B(T):
    lam = 2.0*631515.0/T
    return 2.0*3.435e-30*T*lam**1.970/(1.0 + (lam/2.250)**0.376)**3.720


def voronov(T, dE, P, A, X, K):
    U = dE*EV2ERG/(KB*T)
    return A*(1.0 + P*np.sqrt(U))*U**K*np.exp(-U)/(X + U)


def ci_HI(T):   return voronov(T, 13.6, 0.0, 2.91e-8, 0.232, 0.39)
def ci_HeI(T):  return voronov(T, 24.6, 0.0, 1.75e-8, 0.180, 0.35)
def ci_HeII(T): return voronov(T, 54.4, 1.0, 2.05e-9, 0.265, 0.25)


def gbar_ff(T):
    return 1.1 + 0.34*np.exp(-(5.5 - np.log10(T))**2/3.0)


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


T1_A, T1_T = load_tier1("../../data/atomic/cooling_tier1_h_1.txt")


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
    cool = ne*(nHII*betaB_HII(T) + nHeII_n*beta_HeII_B(T)
               + nHeIII*beta_HeIII_B(T))
    cool += 1.42554e-27*np.sqrt(T)*ne*((nHII + nHeII_n)*gbar_ff(T)
                                       + 4.0*nHeIII*gbar_ff(T))
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
    sH, sHe1, sHe2 = sig_HI(ec), sig_HeI(ec), sig_HeII(ec)
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
print("=== 1D thermal reference (same rates; gbar approximation noted) ===")
r1d, xHI1d, xHeI1d, xHeII1d, Te1d = reference_1d()
i_in = r1d < 2.5
print(f"1D: Te at r=0.5 pc = {Te1d[np.argmin(abs(r1d-0.5))]:.1f} K, "
      f"at 2.0 pc = {Te1d[np.argmin(abs(r1d-2.0))]:.1f} K")
print(f"1D: R(x_HI=0.5) = {r1d[np.argmin(abs(xHI1d-0.5))]:.4f} pc")

f = h5py.File("g2a_thermal_rates.h5", "r")
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

order = np.argsort(r_mc)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
sel_gas = ne_mc > 0
ax1.plot(r_mc[sel_gas], te_mc[sel_gas], ".", ms=1.0, color="C0", alpha=0.25,
         label="MoCHII")
ax1.plot(r1d, Te1d, "k-", lw=1.5, label="1D exact")
ax1.set_xlabel(r"$r$ [pc]");  ax1.set_ylabel(r"$T_e$ [K]")
ax1.set_xlim(0, 4);  ax1.legend(frameon=False, markerscale=10)
ax1.set_title(r"$T_e(r)$, H/He only (no metal cooling yet)")

ax2.plot(r_mc[sel_gas], xHI_mc[sel_gas], ".", ms=1.0, color="C0", alpha=0.25)
ax2.plot(r1d, xHI1d, "k-", lw=1.5)
ax2.set_xlabel(r"$r$ [pc]");  ax2.set_ylabel(r"$x_{\rm HI}$")
ax2.set_xlim(0, 4)
ax2.set_title(r"H ionization with thermal balance")
fig.tight_layout()
fig.savefig("g2a_thermal_check.png", dpi=140)
print("wrote g2a_thermal_check.png")
