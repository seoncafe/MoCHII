#!/usr/bin/env python3
"""Q1 / V3 variance study: ordinary MC (vary iseed) vs RQMC (vary qmc_seed) on
the fixed-state G0 attenuation gate.  Cell-wise Gamma_HI relative error vs the
analytic cell-averaged reference (reused from tests/g0_gamma/check_g0_gamma.py),
plus a radial-profile error where launch RQMC is expected to help most.
"""
import glob, re
import numpy as np
import h5py

EV2ERG = 1.602176634e-12
HP = 6.62607015e-27
KB = 1.380649e-16
PC2CM = 3.0856776e18
ETH_HI, ETH_HeI = 13.598, 24.587
NNU, EMIN, EMAX = 16, 13.598, 100.0
TSTAR, LBAND = 1.0e5, 1.0e38
NH, XHI, XHEI, YHE = 1.0, 0.1, 1.0, 0.1
RSPH_CM = 1.0 * PC2CM
HALF_PC = 1.0 / 32.0


def sigma_vfky96(E, Eth, E0, s0, ya, P, yw, y0, y1):
    E = np.asarray(E, float)
    x = E / E0 - y0
    z = np.sqrt(x * x + y1 * y1)
    Q = 5.5 - 0.5 * P
    s = s0 * ((x - 1.0)**2 + yw**2) * z**(-Q) * (1.0 + np.sqrt(z / ya))**(-P) * 1e-18
    return np.where(E >= Eth, s, 0.0)


def sig_HI(E):
    return sigma_vfky96(E, ETH_HI, 4.298e-1, 5.475e4, 3.288e1, 2.963, 0, 0, 0)


def sig_HeI(E):
    return sigma_vfky96(E, ETH_HeI, 1.361e1, 9.492e2, 1.469, 3.188, 2.039,
                        4.434e-1, 2.136)


def band_bins(nnu, emin, emax, nsub=32):
    edge = np.exp(np.linspace(np.log(emin), np.log(emax), nnu + 1))
    ec = np.sqrt(edge[:-1] * edge[1:])
    lum = np.zeros(nnu)
    for i in range(nnu):
        es = edge[i] + (edge[i + 1] - edge[i]) * (np.arange(nsub) + 0.5) / nsub
        x = es * EV2ERG / (KB * TSTAR)
        lum[i] = np.sum(es**3 / np.expm1(x)) * (edge[i + 1] - edge[i]) / nsub
    lum *= LBAND / lum.sum()
    return ec, lum


def gamma_ana(r_cm, ec, lum):
    kap = NH * (XHI * sig_HI(ec) + YHE * XHEI * sig_HeI(ec))
    tau = np.outer(np.minimum(r_cm, RSPH_CM), kap)
    return np.sum(lum * sig_HI(ec) / (ec * EV2ERG) * np.exp(-tau), axis=1) \
        / (4 * np.pi * r_cm**2)


def gamma_ana_cellavg(xyz_pc, half_pc, ec, lum, nsub=4):
    off = (np.arange(nsub) + 0.5) / nsub - 0.5
    ox, oy, oz = np.meshgrid(off, off, off, indexing="ij")
    pts = np.stack([o.ravel() for o in (ox, oy, oz)])
    acc = np.zeros(xyz_pc.shape[1])
    for k in range(pts.shape[1]):
        p = xyz_pc + 2.0 * half_pc * pts[:, k:k + 1]
        r = np.sqrt((p**2).sum(axis=0)) * PC2CM
        acc += gamma_ana(r, ec, lum)
    return acc / pts.shape[1]


def walltime(logf):
    for ln in open(logf):
        m = re.search(r"Total Execution Time\s*:\s*([0-9.]+)", ln)
        if m:
            return float(m.group(1)) * 60.0
    return np.nan


# reference on the run geometry (all runs share LeafXYZ)
f0 = h5py.File(sorted(glob.glob("tests/qmc/v3/*_rates.h5"))[0], "r")
xyz = f0["LeafXYZ"]["data"][:]
e_bin = f0["E_bin"]["data"][:]
r_pc = np.sqrt((xyz**2).sum(axis=0))
ec, lum = band_bins(NNU, EMIN, EMAX)
assert np.allclose(ec, e_bin, rtol=1e-9), "bins disagree (is ion_align_edges off?)"
ana = gamma_ana_cellavg(xyz, HALF_PC, ec, lum)
sel = (r_pc > 0.15) & (r_pc < 0.85)

# radial-profile bins (a smooth global quantity)
rbins = np.linspace(0.15, 0.85, 15)
ib = np.digitize(r_pc[sel], rbins)
ana_prof = np.array([ana[sel][ib == k].mean() for k in range(1, len(rbins))])


def metrics(fn):
    g = h5py.File(fn, "r")["Gamma_HI"]["data"][:]
    rel = g[sel] / ana[sel] - 1.0
    rms_cell = np.sqrt(np.mean(rel**2))
    bias = rel.mean()
    prof = np.array([g[sel][ib == k].mean() for k in range(1, len(rbins))])
    rms_prof = np.sqrt(np.mean((prof / ana_prof - 1.0)**2))
    return rms_cell, bias, rms_prof


print("=== V3 variance study: Gamma_HI error, MC vs RQMC (fixed-state G0) ===")
for pt, nphot in (("p17", 2**17), ("p20", 2**20)):
    print(f"\n--- nphotons = {nphot} (2^{int(np.log2(nphot))}) ---")
    res = {}
    for mode in ("mc", "qmc"):
        fns = sorted(glob.glob(f"tests/qmc/v3/{mode}_{pt}_*_rates.h5"))
        logs = [fn.replace("_rates.h5", ".log") for fn in fns]
        rc = np.array([metrics(fn) for fn in fns])
        wt = np.array([walltime(lg) for lg in logs])
        res[mode] = rc
        label = "MC  " if mode == "mc" else "RQMC"
        print(f"  {label} replicates (rms_cell, bias, rms_profile):")
        for i, (a, b, c) in enumerate(rc):
            print(f"     r{i+1}: rms_cell={a*100:6.3f}%  bias={b*100:+6.3f}%  "
                  f"rms_prof={c*100:6.3f}%")
        print(f"     mean rms_cell={rc[:,0].mean()*100:6.3f}%  "
              f"mean |bias|={np.abs(rc[:,1]).mean()*100:6.3f}%  "
              f"mean rms_prof={rc[:,2].mean()*100:6.3f}%  "
              f"<wall>={np.nanmean(wt):.2f}s")
    mc, qmc = res["mc"], res["qmc"]
    for j, nm in ((0, "cell-wise"), (2, "radial-profile")):
        rmc = np.sqrt(np.mean(mc[:, j]**2))
        rqm = np.sqrt(np.mean(qmc[:, j]**2))
        gain = (rmc / rqm)**2
        print(f"  >> {nm:14s} RMS: MC={rmc*100:6.3f}%  RQMC={rqm*100:6.3f}%  "
              f"variance ratio (eff. photon gain) = {gain:5.2f}x")
