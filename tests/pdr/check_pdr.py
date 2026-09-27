"""PDR-physics gate: with metal_ne + metal_heat + grain_pe on, the zone
beyond the ionization front must carry (1) n_e from the metal cascade
(n_e ~ A_C n_H, carbon being the dominant electron donor), (2) a T_e
responding to the photoelectric heating instead of pinning at te_min,
while (3) the H II interior stays unchanged against the reference run
without the switches (tests/d_dusty/d_fuv_vol).

Verdict.  Of the three, only (2) is stated sharply enough to test without
choosing a tolerance: "pinning at te_min" is T_e at the lower bound of the
thermal bisection, so (2) passes when no leaf of the PDR shell (2.5-3.9 pc)
has T_e within 1e-6 of par%te_min, read from pdr_L5.in.  (1) "n_e ~ A_C n_H"
and (3) "unchanged" carry no tolerance here or in tests/README.md; they are
printed for reading and do not enter the exit status.  Exit status 0 when (2)
passes, 1 when it fails or an input is missing.

Run:  python3 check_pdr.py   (in tests/pdr)
"""

import os
import re
import sys
import numpy as np
import h5py

REF = "../d_dusty/d_fuv_vol_rates.h5"
RUN = "pdr_L5_rates.h5"
RUN_INPUT = "pdr_L5.in"

for fn in (RUN, REF, RUN_INPUT):
    if not os.path.exists(fn):
        print(f"FAIL: missing input {fn} (run in tests/pdr after pdr_L5.in)")
        sys.exit(1)
te_min = None
for ln in open(RUN_INPUT):
    m = re.match(r"\s*par%te_min\s*=\s*([0-9.eEdD+-]+)", ln, re.I)
    if m:
        te_min = float(m.group(1).lower().replace("d", "e"))
if te_min is None:
    print(f"FAIL: {RUN_INPUT} sets no par%te_min to test criterion (2) against")
    sys.exit(1)


def rd(f, name):
    return f[name + "/data"][:]


with h5py.File(RUN, "r") as f:
    xyz = rd(f, "LeafXYZ")
    ne = rd(f, "n_e")
    te = rd(f, "T_e")
    xhi = rd(f, "x_HI")
    # x_HII as stored (rates files since 2026-09-26), else 1 - x_HI
    xhii = rd(f, "x_HII") if "x_HII" in f else 1.0 - xhi
    xc = rd(f, "x_c_stages")
    if xc.shape[0] != xyz.shape[1]:
        xc = xc.T
with h5py.File(REF, "r") as f:
    xyz0 = rd(f, "LeafXYZ")
    ne0 = rd(f, "n_e")
    te0 = rd(f, "T_e")
    xhi0 = rd(f, "x_HI")

r = np.sqrt((xyz**2).sum(axis=0))
r0 = np.sqrt((xyz0**2).sum(axis=0))

print(f"{'r[pc]':>6} {'x_HII':>8} {'n_e':>10} {'n_e(ref)':>10} "
      f"{'T_e':>8} {'T_e(ref)':>9} {'C II':>7}")
edges = np.linspace(0.0, 4.0, 21)
for i in range(20):
    m = (r >= edges[i]) & (r < edges[i + 1])
    m0 = (r0 >= edges[i]) & (r0 < edges[i + 1])
    if not m.any():
        continue
    mid = 0.5 * (edges[i] + edges[i + 1])
    print(f"{mid:6.1f} {np.median(xhii[m]):8.4f} {np.median(ne[m]):10.3e} "
          f"{np.median(ne0[m0]):10.3e} {np.median(te[m]):8.0f} "
          f"{np.median(te0[m0]):9.0f} {np.median(xc[m, 1]):7.3f}")

# headline checks (run and ref may be on different grid levels)
pdr = (r > 2.5) & (r < 3.9)
pdr0 = (r0 > 2.5) & (r0 < 3.9)
hii = (r < 1.8)
hii0 = (r0 < 1.8)
print(f"\nPDR shell (2.5-3.9 pc):")
print(f"  <n_e>  = {ne[pdr].mean():.3e} cm^-3  "
      f"(A_C n_H = {2.2e-4*100:.3e}; ref run: {ne0[pdr0].mean():.3e})")
print(f"  <T_e>  = {te[pdr].mean():8.0f} K  (ref run: {te0[pdr0].mean():8.0f} K)")
print(f"H II interior (r < 1.8 pc):")
print(f"  |<x_HI>run - <x_HI>ref| = "
      f"{abs(xhi[hii].mean()-xhi0[hii0].mean()):.2e}")
print(f"  |<T_e>run/<T_e>ref - 1| = "
      f"{abs(te[hii].mean()/te0[hii0].mean()-1):.2e}")

pinned = int(np.count_nonzero(np.abs(te[pdr] / te_min - 1.0) < 1e-6))
ok = pinned == 0
print(f"\n(2) PDR-shell leaves at te_min = {te_min:g} K: {pinned} of "
      f"{int(pdr.sum())} (T_e {te[pdr].min():.0f}-{te[pdr].max():.0f} K)")
print("(1), (3): no tolerance; reported above, not judged")
print(f"PDR GATE (criterion 2): {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
