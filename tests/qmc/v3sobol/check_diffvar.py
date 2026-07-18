#!/usr/bin/env python3
"""V3(e): diffuse-field replicate scatter, ordinary MC vs RQMC.

For each replicate compute V_ion (H-ionized volume) and R_eff over the gas
sphere, then the replicate std for the MC set (vary iseed) and the RQMC set
(vary qmc_seed).  The RQMC scatter should be <= the MC scatter for these smooth
global quantities (docs predicts a smaller / less stable gain for diffuse).
"""
import os
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
PC = 3.0856776e18
RMAX = 1.9


def vion_reff(fn):
    f = h5py.File(fn, "r")
    d = {k: f[k]["data"][:] for k in f
         if isinstance(f[k], h5py.Group) and "data" in f[k]}
    xyz = d["LeafXYZ"]
    vol = (d["LeafSize"] * PC) ** 3
    r = np.sqrt((xyz ** 2).sum(axis=0))
    m = r <= RMAX
    Vion = np.sum(vol[m] * (1.0 - d["x_HI"][m])) / PC**3
    Reff = (3.0 * Vion / (4.0 * np.pi)) ** (1.0 / 3.0)
    return Vion, Reff


mc = [vion_reff(os.path.join(HERE, f"dv_mc_{s}_rates.h5")) for s in (11, 22, 33)]
qm = [vion_reff(os.path.join(HERE, f"dv_qmc_{s}_rates.h5")) for s in (101, 202, 303)]
mc = np.array(mc)
qm = np.array(qm)

print("=== V3(e) diffuse-field replicate scatter: MC (iseed) vs RQMC (qmc_seed) ===")
print("    diffuse Stromgren, case A + diffuse, aligned grid, 5e5 photons, 40 iter, front at R_eff~0.72 pc")
for j, name in ((0, "V_ion [pc^3]"), (1, "R_eff [pc]")):
    print(f"\n  {name}:")
    print(f"    MC   replicates: {mc[:,j]}")
    print(f"    RQMC replicates: {qm[:,j]}")
    smc, sqm = mc[:, j].std(ddof=1), qm[:, j].std(ddof=1)
    rmc = smc / mc[:, j].mean()
    rqm = sqm / qm[:, j].mean()
    print(f"    MC   std = {smc:.4e}  (rel {rmc*100:.3f}%)")
    print(f"    RQMC std = {sqm:.4e}  (rel {rqm*100:.3f}%)")
    ratio = (smc / sqm) if sqm > 0 else float("inf")
    verdict = "RQMC <= MC (variance reduced)" if sqm <= smc else \
              "RQMC > MC (no gain here - expected weaker/less stable for diffuse)"
    print(f"    std ratio MC/RQMC = {ratio:.2f}x  -> {verdict}")
