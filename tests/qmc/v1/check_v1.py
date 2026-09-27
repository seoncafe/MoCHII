#!/usr/bin/env python3
"""V1: HEAD (pre-QMC) vs NEW binary, both -fp-model precise, launch='random'.

(a),(b): worst |diff| over every rates dataset must be 0 (the random path is
arithmetically unchanged).  (c): on the log grid the Task-0 binary-search bin
mapping is unchanged too, so any difference must stay at the rounding floor,
max relative difference < 1e-13 (the 8-rank diffuse run differs there by
-ipo code generation, not by a bin flip).  (d): the aligned-grid diffuse run
SHOULD differ - that is the bug fix; quantify V_ion and the He I/He II
fractions HEAD vs NEW.

Verdict: PASS when (a), (b) and (c) hold; (d) is a measurement and does not
enter it.  Exit status 0 on PASS, 1 on FAIL or when
a rates file is missing.  run_v1.sh names two binaries under /tmp that no
longer exist, so the rates files are not in the tree and this check cannot
be run as it stands.  Run from the repository root.
"""
import glob, os, sys
import numpy as np
import h5py

PC = 3.0856776e18


def datasets(fn):
    out = {}
    with h5py.File(fn, "r") as f:
        for k in f:
            if isinstance(f[k], h5py.Group) and "data" in f[k]:
                out[k] = f[k]["data"][:]
    return out


def worst_diff(tag, base):
    a = datasets(f"tests/qmc/v1/{base}_head_rates.h5")
    b = datasets(f"tests/qmc/v1/{base}_new_rates.h5")
    keys = sorted(set(a) & set(b))
    worst, wk, wrel = 0.0, "", 0.0
    for k in keys:
        x, y = np.asarray(a[k], float), np.asarray(b[k], float)
        if x.shape != y.shape:
            print(f"  {tag}: shape mismatch {k}")
            continue
        d = np.abs(x - y)
        m = float(np.nanmax(d)) if d.size else 0.0
        rel = d / np.where(np.abs(x) > 0, np.abs(x), 1.0)
        mr = float(np.nanmax(rel)) if rel.size else 0.0
        if m > worst:
            worst, wk = m, k
        wrel = max(wrel, mr)
    if worst == 0.0:
        status = "BIT-IDENTICAL"
    elif wrel < 1e-13:
        status = (f"rounding floor (worst |diff|={worst:.3e} in {wk}, max rel={wrel:.1e} "
                  f"= benign -ipo codegen; bins verified identical by check_task0.py)")
    else:
        status = f"DIFFERS (worst |diff|={worst:.3e} in {wk}, max rel={wrel:.1e})"
    print(f"  {tag} [{base}]: {len(keys)} datasets -> {status}")
    return worst, wrel


def vion_he(fn, rmax=1.9):
    d = datasets(fn)
    xyz = d["LeafXYZ"]
    sz = d["LeafSize"]
    vol = (sz * PC) ** 3
    # x_HII as stored (rates files since 2026-09-26), else 1 - x_HI
    xHII = d["x_HII"] if "x_HII" in d else 1.0 - d["x_HI"]
    xHeI = d["x_HeI"]
    xHeII = d["x_HeII"]
    r = np.sqrt((xyz ** 2).sum(axis=0))            # leaf-center radius [pc]
    inbox = r <= rmax                              # inside the gas sphere
    Vion = np.sum(vol[inbox] * xHII[inbox])
    Reff = (3.0 * Vion / (4.0 * np.pi)) ** (1.0 / 3.0) / PC
    # volume-weighted mean He fractions over the gas sphere
    w = vol[inbox]
    meanHeI = np.sum(w * xHeI[inbox]) / np.sum(w)
    meanHeII = np.sum(w * xHeII[inbox]) / np.sum(w)
    # x_HeIII as stored (rates files since 2026-09-26), else the complement
    xHeIII = (d["x_HeIII"] if "x_HeIII" in d
              else np.clip(1 - xHeI - xHeII, 0, None))
    meanHeIII = np.sum(w * xHeIII[inbox]) / np.sum(w)
    return Vion / PC**3, Reff, meanHeI, meanHeII, meanHeIII


need = [f"tests/qmc/v1/{b}_{w}_rates.h5"
        for b in ("v1a_multi", "v1b_ext", "v1c_diff_logedge", "v1d_diff_aligned")
        for w in ("head", "new")]
missing = [fn for fn in need if not os.path.exists(fn)]
if missing:
    print("missing rates files (run_v1.sh writes them):")
    for fn in missing:
        print("  " + fn)
    print("\nV1 GATE: FAIL (missing input)")
    sys.exit(1)

print("=== V1 regression: HEAD vs NEW, launch='random', -fp-model precise ===")
wa, _ = worst_diff("(a) multi-source np=1", "v1a_multi")
wb, _ = worst_diff("(b) external    np=1", "v1b_ext")
wc, wcr = worst_diff("(c) diffuse LOG np=8", "v1c_diff_logedge")
print(f"  >> (a),(b) worst |diff| = {max(wa,wb):.3e} (must be 0): "
      f"{'PASS' if max(wa,wb)==0.0 else 'FAIL'}")
print(f"  >> (c) worst |diff| = {wc:.3e}, max rel = {wcr:.1e}: "
      f"{'PASS (rounding floor - bin mapping identical, not a bin flip)' if wcr<1e-13 else 'FAIL'}")
ok = max(wa, wb) == 0.0 and wcr < 1e-13

print("\n=== V1(d) diffuse ALIGNED grid: physical impact of the bin fix ===")
vh_head = vion_he("tests/qmc/v1/v1d_diff_aligned_head_rates.h5")
vh_new = vion_he("tests/qmc/v1/v1d_diff_aligned_new_rates.h5")
lab = ["V_ion[pc^3]", "R_eff[pc]", "<x_HeI>", "<x_HeII>", "<x_HeIII>"]
print(f"  {'quantity':14s} {'HEAD(buggy)':>14s} {'NEW(fixed)':>14s} {'rel.change':>12s}")
for name, h, n in zip(lab, vh_head, vh_new):
    rc = (n / h - 1.0) if h != 0 else float("nan")
    print(f"  {name:14s} {h:14.6g} {n:14.6g} {rc*100:+11.3f}%")
dw, dwr = worst_diff("(d) diffuse ALIGNED", "v1d_diff_aligned")
print(f"  (the (d) datasets DIFFER by construction - the fix; worst |diff| above)")

print(f"\nV1 GATE ((a)-(c)): {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
