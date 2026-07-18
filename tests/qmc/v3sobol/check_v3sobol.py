#!/usr/bin/env python3
"""V3: the recorded physics gates re-run under launch_sequence='sobol'.

Reuses the EXISTING gate checkers (tests/multi_src/check_multi.report_superpos
and tests/ext_field/check_ext.report) on the sobol rates files, so the pass
criteria are identical to the 'random' gates.  Also (d) MPI-count independence
of the sobol launch: the same run at np=4 and np=8 must agree to the ALLREDUCE
floor.
"""
import os, sys
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tests", "multi_src"))
sys.path.insert(0, os.path.join(ROOT, "tests", "ext_field"))
import check_multi
import check_ext

R = lambda t: os.path.join(HERE, f"{t}_rates.h5")

print("=== V3(a) two-point superposition, launch='sobol' ===")
print("    (criterion: median<3% and |mean-1|<0.02; recorded 'random': "
      "median 0.76%, mean 0.9989)")
a = check_multi.report_superpos("two-point (sobol)", R("v3s_twop"),
                                srcs=[(-1.0, 0.0, 0.0), (1.0, 0.5, 0.0)],
                                lums=[3.0e37, 1.0e37])

print("\n=== V3(b) external gates, launch='sobol' ===")
print("    (criterion: |<J>/J - 1| < 0.05; recorded 'random': rec 0.9981, sph 0.9991)")
def rec_mask(xyz, p):
    return np.ones(xyz.shape[1], bool)
def sph_mask(xyz, p):
    r = np.sqrt((xyz ** 2).sum(axis=0))
    return r < 0.6 * 1.2
b1 = check_ext.report("external_rec (sobol)", R("v3s_extrec"),
                      os.path.join(HERE, "v3s_extrec.in"), rec_mask)
b2 = check_ext.report("external_sph (sobol)", R("v3s_extsph"),
                      os.path.join(HERE, "v3s_extsph.in"), sph_mask)

print("\n=== V3(c) mixed point + external, launch='sobol' ===")
print("    (criterion: median<3% and |mean-1|<0.02; recorded 'random': median 0.75%)")
c = check_multi.report_superpos("mixed (sobol)", R("v3s_mixed"),
                                srcs=[(0.0, 0.0, 0.0)], lums=[1.0e37], Jext=1.0e-5)

print("\n=== V3(d) MPI-count independence of the sobol launch (np=4 vs np=8) ===")
def worst(tag, t4, t8):
    def ds(fn):
        f = h5py.File(fn, "r")
        return {k: f[k]["data"][:] for k in f
                if isinstance(f[k], h5py.Group) and "data" in f[k]}
    a, b = ds(R(t4)), ds(R(t8))
    w, wk, wr = 0.0, "", 0.0
    for k in sorted(set(a) & set(b)):
        x, y = np.asarray(a[k], float), np.asarray(b[k], float)
        if x.shape != y.shape:
            continue
        d = np.abs(x - y)
        m = float(np.nanmax(d)) if d.size else 0.0
        rel = float(np.nanmax(d / np.where(np.abs(x) > 0, np.abs(x), 1.0))) if d.size else 0.0
        if m > w:
            w, wk = m, k
        wr = max(wr, rel)
    ok = wr < 1e-12
    print(f"  {tag}: worst |diff| np4-vs-np8 = {w:.3e} (max rel {wr:.1e} in {wk}) "
          f"-> {'PASS (ALLREDUCE floor)' if ok else 'FAIL'}")
    return ok
d1 = worst("two-point", "v3s_twop_np4", "v3s_twop_np8")
d2 = worst("external",  "v3s_extrec_np4", "v3s_extrec_np8")

allok = a and b1 and b2 and c and d1 and d2
print(f"\nV3 GATE: {'PASS' if allok else 'FAIL'}")
