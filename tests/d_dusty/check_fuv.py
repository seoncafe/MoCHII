"""Check the add_fuv option: radial ionization structure of H and of the
low-threshold metals (C, Mg) across the I-front.  Without FUV the metals
beyond the front recombine toward neutral; with add_fuv the FUV photons
(dust extinction only) keep them singly ionized there (PDR-side C II and
Mg II).

Run:  python3 check_fuv.py d_fuv_opt_rates.h5
"""

import sys
import numpy as np
import h5py

fname = sys.argv[1] if len(sys.argv) > 1 else "d_fuv_opt_rates.h5"

with h5py.File(fname, "r") as f:
    xyz = f["LeafXYZ/data"][:]           # (3, nleaf), pc
    xHI = f["x_HI/data"][:]
    xc = f["x_c_stages/data"][:]         # (nleaf, nstage)
    xmg = f["x_mg_stages/data"][:]
r = np.sqrt((xyz**2).sum(axis=0))

edges = np.linspace(0.0, 4.0, 41)
mid = 0.5 * (edges[1:] + edges[:-1])
print(f"{'r[pc]':>6} {'x_HII':>8} {'C I':>8} {'C II':>8} {'C III':>8}"
      f" {'Mg I':>8} {'Mg II':>8} {'Mg III':>8}")
for i in range(len(mid)):
    m = (r >= edges[i]) & (r < edges[i + 1])
    if not m.any():
        continue
    print(f"{mid[i]:6.2f} {1-xHI[m].mean():8.4f}"
          f" {xc[m, 0].mean():8.4f} {xc[m, 1].mean():8.4f}"
          f" {xc[m, 2].mean():8.4f}"
          f" {xmg[m, 0].mean():8.4f} {xmg[m, 1].mean():8.4f}"
          f" {xmg[m, 2].mean():8.4f}")

# headline numbers: front radius and the PDR-side C II / Mg II
xHII = 1.0 - xHI
vion = (xHII > 0.5)
# R_eff from the ionized volume (uniform leaves at level 6: equal volumes)
R_eff = (3.0 * xHII.sum() * (8.0 / xHI.size * 8.0**2) / (4.0 * np.pi))**(1/3)
pdr = (r > 3.3) & (r < 3.9)
print(f"\nR_eff(V_ion, equal-volume approx) = {R_eff:.3f} pc")
print(f"PDR shell 3.3-3.9 pc:  <x_HII> = {xHII[pdr].mean():.4f}")
print(f"  C  I/II/III = {xc[pdr,0].mean():.3f} {xc[pdr,1].mean():.3f} "
      f"{xc[pdr,2].mean():.3f}")
print(f"  Mg I/II/III = {xmg[pdr,0].mean():.3f} {xmg[pdr,1].mean():.3f} "
      f"{xmg[pdr,2].mean():.3f}")
