#!/usr/bin/env python3
"""G1 gate: Stromgren sphere vs the analytic R_S and a 1D exact reference.

Reference 1 (sharp edge):  R_S = (3 Q_H / (4 pi nH^2 alpha_B))^(1/3).
Reference 2 (exact): 1D radial marching with the SAME frequency bins,
cross sections, and rate coefficients as the Fortran — at each radius the
local equilibrium x(r) from the attenuated bin fluxes, accumulating the
optical depth of each shell; iterated to self-consistency (the marching
uses x from the previous sweep; converges in a few sweeps).

Gate criteria:
  (1) R_eff = (3 V_ion / 4 pi)^(1/3) from the MC leaves (V_ion = sum of
      x_HII * V_cell) within 1% of the 1D reference R_eff;
  (2) x_HI(r) profile tracks the 1D reference (shell medians);
  (3) refined-grid run (levels 4-6) matches the uniform level-6 run
      (AMR <-> uniform methodology, PLAN section 5).

Run after:  mpirun -np 8 ../../MoCHII.x g1_stromgren.in
            mpirun -np 8 ../../MoCHII.x g1_stromgren_ref.in
"""
import numpy as np
import h5py

EV2ERG = 1.602176634e-12
KB     = 1.380649e-16
PC2CM  = 3.0856776e18
ETH_HI, ETH_HeI, ETH_HeII = 13.598, 24.587, 54.416

NNU, EMIN, EMAX = 32, 13.598, 100.0
TSTAR, LBAND    = 4.0e4, 3.177837e38
NH, YHE, TE     = 100.0, 0.1, 1.0e4
RSPH_PC         = 4.0


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


def hui_gnedin_B(T, TTR, pref, lam0, p1, p2):
    lam = 2.0*TTR/T
    return pref*lam**p1/(1.0 + (lam/lam0)**0.407)**p2


def alphaB_HII(T):   return hui_gnedin_B(T, 157807.0, 2.753e-14, 2.740, 1.500, 2.242)
def alphaB_HeII(T):  return 1.260e-14*(2.0*285335.0/T)**0.750
def alphaB_HeIII(T): return 2.0*hui_gnedin_B(T, 631515.0, 2.753e-14, 2.740, 1.500, 2.242)


def voronov(T, dE, P, A, X, K):
    U = dE*EV2ERG/(KB*T)
    return A*(1.0 + P*np.sqrt(U))*U**K*np.exp(-U)/(X + U)


def ci_HI(T):   return voronov(T, 13.6, 0.0, 2.91e-8, 0.232, 0.39)
def ci_HeI(T):  return voronov(T, 24.6, 0.0, 1.75e-8, 0.180, 0.35)
def ci_HeII(T): return voronov(T, 54.4, 1.0, 2.05e-9, 0.265, 0.25)


def band_bins(nnu=NNU, emin=EMIN, emax=EMAX, nsub=32):
    edge = np.exp(np.linspace(np.log(emin), np.log(emax), nnu + 1))
    ec = np.sqrt(edge[:-1]*edge[1:])
    lum = np.zeros(nnu)
    for i in range(nnu):
        es = edge[i] + (edge[i+1]-edge[i])*(np.arange(nsub)+0.5)/nsub
        lum[i] = np.sum(es**3/np.expm1(es*EV2ERG/(KB*TSTAR)))*(edge[i+1]-edge[i])/nsub
    lum *= LBAND/lum.sum()
    return ec, lum


def equilibrium(gH, gHe1, gHe2, T=TE, nH=NH, yHe=YHE):
    """Same damped fixed point as ion_balance_mod (vectorized over r)."""
    aH, aHe2, aHe3 = alphaB_HII(T), alphaB_HeII(T), alphaB_HeIII(T)
    cH, cHe1, cHe2 = ci_HI(T), ci_HeI(T), ci_HeII(T)
    xHI = np.ones_like(gH)*0.5
    ne = nH*(0.5 + yHe*0.5)*np.ones_like(gH)
    for _ in range(300):
        ne_old = ne.copy()
        xHI = aH*ne/(gH + (cH + aH)*ne)
        r1 = (gHe1 + cHe1*ne)/(aHe2*ne)
        r2 = (gHe2 + cHe2*ne)/(aHe3*ne)
        xHeI = 1.0/(1.0 + r1 + r1*r2)
        xHeII = xHeI*r1
        xHeIII = xHeII*r2
        ne = nH*((1.0 - xHI) + yHe*(xHeII + 2.0*xHeIII))
        ne = np.maximum(0.5*(ne + ne_old), 1e-12*nH)
        if np.abs(ne/ne_old - 1.0).max() < 1e-12:
            break
    return xHI, xHeI, xHeII


def reference_1d(nr=20000, nsweep=4):
    """Gauss-Seidel radial marching with the same bins/rates.

    Radiation flows outward only (case B on-the-spot: no diffuse field), so
    solving each shell with the ALREADY-UPDATED optical depth of the inner
    shells makes one sweep essentially exact; a Jacobi sweep (previous-sweep
    opacities) advances the I-front only ~one absorption length per sweep
    and never converges in practice.  The shell's own attenuation to its
    center uses the current-shell solution iterated to consistency.
    Returns r_pc, x profiles.
    """
    ec, lum = band_bins()
    sH, sHe1, sHe2 = sig_HI(ec), sig_HeI(ec), sig_HeII(ec)
    hnu = ec*EV2ERG
    dr_pc = RSPH_PC/nr
    dr_cm = dr_pc*PC2CM
    r = (np.arange(nr) + 0.5)*dr_pc
    r_cm = r*PC2CM
    xHI = np.ones(nr); xHeI = np.ones(nr); xHeII = np.zeros(nr)
    for sweep in range(nsweep):
        x_prev = xHI.copy()
        tau = np.zeros(len(ec))
        for i in range(nr):
            xh, xh1, xh2 = xHI[i], xHeI[i], xHeII[i]
            for _ in range(3):     # self-shell attenuation consistency
                kap_i = NH*(xh*sH + YHE*(xh1*sHe1 + xh2*sHe2))
                att = np.exp(-(tau + 0.5*kap_i*dr_cm))
                flux = lum*att/(4*np.pi*r_cm[i]**2)
                gH   = float((flux*sH/hnu).sum())
                gHe1 = float((flux*sHe1/hnu).sum())
                gHe2 = float((flux*sHe2/hnu).sum())
                xh, xh1, xh2 = (float(v[0]) for v in equilibrium(
                    np.array([gH]), np.array([gHe1]), np.array([gHe2])))
            xHI[i], xHeI[i], xHeII[i] = xh, xh1, xh2
            kap_i = NH*(xh*sH + YHE*(xh1*sHe1 + xh2*sHe2))
            tau = tau + kap_i*dr_cm
        if np.abs(xHI - x_prev).max() < 1e-8:
            break
    return r, xHI, xHeI, xHeII


def load_run(fname):
    f = h5py.File(fname, "r")
    xyz = f["LeafXYZ"]["data"][:]
    return dict(
        r=np.sqrt((xyz**2).sum(axis=0)),
        xHI=f["x_HI"]["data"][:], xHeI=f["x_HeI"]["data"][:],
        xHeII=f["x_HeII"]["data"][:], ne=f["n_e"]["data"][:],
        nleaf=xyz.shape[1])


def r_eff(run, half_cell_arr):
    """Effective Stromgren radius from the ionized volume of leaves with
    gas (nH > 0 within the sphere)."""
    vol = (2.0*half_cell_arr)**3
    inside = run["ne"] > 0            # leaves with gas
    vion = np.sum((1.0 - run["xHI"][inside])*vol[inside])
    return (3.0*vion/(4.0*np.pi))**(1.0/3.0)


# ------------------------------------------------------------------ main --
print("=== 1D exact reference (same bins/rates) ===")
r1d, xHI1d, xHeI1d, xHeII1d = reference_1d()
RS_sharp = (3.0e49/(4*np.pi*NH**2*alphaB_HII(TE)))**(1/3)/PC2CM
i05 = np.argmin(np.abs(xHI1d - 0.5))
vion_1d = np.sum((1.0 - xHI1d)*4*np.pi*(r1d*PC2CM)**2)*(r1d[1]-r1d[0])*PC2CM
R_eff_1d = (3.0*vion_1d/(4*np.pi))**(1/3)/PC2CM
print(f"sharp-edge analytic R_S        = {RS_sharp:.4f} pc")
print(f"1D reference: R(x_HI=0.5)      = {r1d[i05]:.4f} pc")
print(f"1D reference: R_eff (V_ion)    = {R_eff_1d:.4f} pc")

runs = {}
for tag, fname, hcell in (("uniform L6", "g1_stromgren_rates.h5", None),
                          ("refined 4-6", "g1_stromgren_ref_rates.h5", None),
                          ("front 4-6", "g1_stromgren_front_rates.h5", None)):
    try:
        run = load_run(fname)
    except FileNotFoundError:
        print(f"[{tag}] {fname} not found; skipped")
        continue
    # cell half sizes from leaf spacing: uniform grid -> constant; refined ->
    # recover from the level structure via nearest-neighbor spacing is
    # fragile, so store half size from the run resolution instead:
    if run["nleaf"] == 262144:
        half = np.full(run["nleaf"], 4.0/64)      # level 6, boxlen 8
    else:
        # refined grid: infer level from cell volume via nH>0 leaf spacing —
        # use the fact that make_amr_sphere emits levels 4-6 on boxlen 8:
        # match each leaf to the smallest grid step consistent with center
        # coordinates: centers are odd multiples of (4/2^L)/ ... simpler:
        # centers at level L sit at (k+0.5)*8/2^L - 4; test divisibility.
        half = np.empty(run["nleaf"])
        xyz_r = None
        f = h5py.File(fname, "r")
        xyz = f["LeafXYZ"]["data"][:]
        for L in (6, 5, 4):
            step = 8.0/2**L
            on = np.all(np.abs(((xyz + 4.0)/step - 0.5) -
                               np.round((xyz + 4.0)/step - 0.5)) < 1e-6, axis=0)
            half[on] = step/2
    runs[tag] = (run, half)
    Re = r_eff(run, half)
    dev = Re/R_eff_1d - 1.0
    print(f"[{tag:12s}] nleaf = {run['nleaf']:7d}   R_eff = {Re:.4f} pc   "
          f"dev vs 1D = {dev*100:+.2f}%")

# --- profile comparison plot ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.plot(r1d, xHI1d, "k-", lw=1.5, label="1D exact")
ax1.plot(r1d, 1.0 - xHI1d, "k--", lw=1.0)
colors = {"uniform L6": "C0", "refined 4-6": "C1", "front 4-6": "C2"}
for tag, (run, half) in runs.items():
    sel = run["ne"] > 0
    ax1.plot(run["r"][sel], run["xHI"][sel], ".", ms=1.0,
             color=colors[tag], alpha=0.25, label=tag)
ax1.axvline(RS_sharp, color="0.5", ls=":", lw=1, label=r"$R_S$ (sharp)")
ax1.set_xlabel(r"$r$ [pc]");  ax1.set_ylabel(r"$x_{\rm HI}$, $x_{\rm HII}$")
ax1.set_xlim(0, 4);  ax1.legend(frameon=False, markerscale=10)
ax1.set_title(r"H ionization structure")

ax2.plot(r1d, xHeI1d, "k-", lw=1.5, label=r"1D $x_{\rm HeI}$")
ax2.plot(r1d, xHeII1d, "k--", lw=1.5, label=r"1D $x_{\rm HeII}$")
for tag, (run, half) in runs.items():
    sel = run["ne"] > 0
    ax2.plot(run["r"][sel], run["xHeI"][sel], ".", ms=1.0,
             color=colors[tag], alpha=0.25)
    ax2.plot(run["r"][sel], run["xHeII"][sel], ".", ms=1.0,
             color=colors[tag], alpha=0.25)
ax2.set_xlabel(r"$r$ [pc]");  ax2.set_ylabel(r"$x_{\rm HeI}$, $x_{\rm HeII}$")
ax2.set_xlim(0, 4);  ax2.legend(frameon=False)
ax2.set_title(r"He ionization structure")
fig.tight_layout()
fig.savefig("g1_stromgren_check.png", dpi=140)
print("wrote g1_stromgren_check.png")
