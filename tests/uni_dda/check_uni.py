"""Uniform-DDA gate: grid_type='uniform' (raster + integer-step DDA
traversal) must reproduce grid_type='amr' (single-level octree) on the
same Stromgren problem with the same seed, leaf by leaf after matching
positions; and the timing comparison measures what the DDA path buys.

Run:  python3 check_uni.py
"""

import re
import numpy as np
import h5py


def load(base):
    with h5py.File(f"{base}_rates.h5", "r") as f:
        xyz = f["LeafXYZ/data"][:]
        xhi = f["x_HI/data"][:]
        ne = f["n_e/data"][:]
    return xyz, xhi, ne


xyz_a, xhi_a, ne_a = load("uni_amr")
xyz_d, xhi_d, ne_d = load("uni_dda")

# match leaves by position (tree order vs raster order)
def keys(xyz):
    # cell centers sit at half-integer multiples of dcell: floor-index
    # them (np.round's half-even rule collides neighboring cells).
    q = np.floor((xyz + 4.0) / 0.125).astype(int)   # L6: dcell = 8/64
    return q[0] + 1000 * (q[1] + 1000 * q[2])


ka = np.argsort(keys(xyz_a))
kd = np.argsort(keys(xyz_d))
dx = np.abs(xhi_a[ka] - xhi_d[kd])
print(f"leaf-matched |d x_HI|: max = {dx.max():.3e}, "
      f"median = {np.median(dx):.3e}")

vion_a = (1.0 - xhi_a).sum()
vion_d = (1.0 - xhi_d).sum()
print(f"V_ion ratio (dda/amr) = {vion_d/vion_a:.6f}")
R_a = (3.0 * vion_a * (8.0**3 / xhi_a.size) / (4 * np.pi))**(1 / 3)
R_d = (3.0 * vion_d * (8.0**3 / xhi_d.size) / (4 * np.pi))**(1 / 3)
print(f"R_eff: amr {R_a:.4f} pc, dda {R_d:.4f} pc "
      f"(G1 gate reference 3.0534 pc at this setup)")

# timing: per-iteration transport wall time from the progress stamps
def times(log):
    t = [float(m.group(1)) for m in
         re.finditer(r"transport\.\.\.\s+@\s+([0-9.]+) mins", open(log).read())]
    tot = re.search(r"Total Execution Time\s*:\s*([0-9.]+)", open(log).read())
    return np.diff(t), float(tot.group(1)) if tot else np.nan


it_a, tot_a = times("uni_amr.log")
it_d, tot_d = times("uni_dda.log")
print(f"per-iteration wall time [min]: amr median {np.median(it_a):.3f}, "
      f"dda median {np.median(it_d):.3f}  -> speedup x{np.median(it_a)/np.median(it_d):.2f}")
print(f"total run [min]: amr {tot_a:.2f}, dda {tot_d:.2f}")
ok = dx.max() < 5e-3 and abs(vion_d/vion_a - 1) < 1e-3
print("GATE:", "PASS" if ok else "FAIL")
