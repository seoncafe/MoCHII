"""Car-grid namelist + density-model gate.

Two checks on the same R=4 pc, nH=100 Stromgren sphere:
  1. car_namelist (grid built from par%nx/xmax + par%nH_const + par%rmax,
     NO file) must reproduce uni_dda (the file-built 'car' grid) leaf by
     leaf -- same raster order, same seed -> bit-identical.
  2. car_amr_ovr (the octree file grid with nH set via par%nH_const and the
     sphere via par%rmax, an override) must reproduce uni_amr (the plain
     octree run whose file already carries the same nH sphere).

Run:  python3 check_car.py
"""

import numpy as np
import h5py


def load(base):
    with h5py.File(f"{base}_rates.h5", "r") as f:
        return f["LeafXYZ/data"][:], f["x_HI/data"][:], f["n_e/data"][:]


def keys(xyz):
    q = np.floor((xyz + 4.0) / 0.125).astype(int)   # L6: dcell = 8/64
    return q[0] + 1000 * (q[1] + 1000 * q[2])


def compare(name, a, b, tol):
    xyz_a, xhi_a, ne_a = load(a)
    xyz_b, xhi_b, ne_b = load(b)
    ka = np.argsort(keys(xyz_a))
    kb = np.argsort(keys(xyz_b))
    dx = np.abs(xhi_a[ka] - xhi_b[kb])
    va, vb = (1.0 - xhi_a).sum(), (1.0 - xhi_b).sum()
    print(f"{name}: leaf-matched |d x_HI| max = {dx.max():.3e}, "
          f"median = {np.median(dx):.3e}; V_ion ratio = {vb/va:.8f}")
    ok = dx.max() <= tol and abs(vb / va - 1) < 1e-6
    print(f"  GATE: {'PASS' if ok else 'FAIL'}  (tol {tol:.0e})")
    return ok


# 1. namelist car vs file car -> bit-identical
ok1 = compare("namelist-car vs file-car", "car_namelist", "uni_dda", 0.0)
# 2. amr density override vs plain octree -> bit-identical
ok2 = compare("amr-override vs octree", "car_amr_ovr", "uni_amr", 0.0)

print("\nOVERALL:", "PASS" if (ok1 and ok2) else "FAIL")
