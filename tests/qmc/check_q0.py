#!/usr/bin/env python3
"""Q0 gate for qmc_mod: Sobol generator correctness + Owen-scramble balance.

(a) Unscrambled codes vs scipy.stats.qmc.Sobol(d=3, scramble=False): scipy uses
    the same Joe-Kuo direction numbers.  For N=2^m the first N Sobol points are
    the same SET regardless of the (direct vs Gray-code) enumeration order, so
    we compare the sorted set of 3-D integer codes exactly.

(b) Owen-scramble balance.  A nested-uniform (Owen) scramble maps a digital net
    to another net with IDENTICAL (t,m,s) equidistribution parameters, so:
      - every 1-D dyadic interval of size 2^-j (j<=m) holds exactly N/2^j points
        (each dimension is a t=0 net in 1-D);
      - the (dim1,dim2) = (van der Corput, first Sobol) projection is a
        (0,m,2)-net: exactly one point per elementary box of area 2^-m;
      - for EVERY 2-D projection and every dyadic box shape, the scrambled net
        has the SAME box-occupancy histogram as the unscrambled net.
    The last item is the exact preservation property; note the (dim1,dim3) and
    (dim2,dim3) Sobol projections are NOT (0,2)-nets (t=1), so "one point per
    box" is expected to fail there for BOTH the raw and the scrambled net -- the
    requirement is that the scramble does not change the occupancy.

(c) Repeat (b) for several seeds; the point sets must differ.
"""
import numpy as np
from scipy.stats import qmc

N = 1024
M = 10
assert 2**M == N

# ---- (a) unscrambled codes vs scipy -----------------------------------
raw = np.loadtxt("tests/qmc/q0_raw.txt", dtype=np.int64)
codes = raw[:, 1:4]                                   # (N,3) 32-bit codes
u_raw = (codes + 0.5) / 2.0**32
sob = qmc.Sobol(d=3, scramble=False).random(N)        # scipy, first N points
scipy_codes = np.rint(sob * 2.0**32).astype(np.int64)

def sort_rows(a):
    return a[np.lexsort(a.T[::-1])]

match = np.array_equal(sort_rows(codes), sort_rows(scipy_codes))
maxdiff = int(np.abs(sort_rows(codes) - sort_rows(scipy_codes)).max())
print("=== Q0(a) unscrambled Sobol vs scipy (Joe-Kuo direction numbers) ===")
print(f"  set-equal (sorted 3-D integer codes)  : {match}")
print(f"  max |code diff| over sorted set        : {maxdiff}")
a_ok = match and maxdiff == 0


def box_hist(u, a1, a2, j1, j2):
    bx = np.floor(u[:, a1] * (2**j1)).astype(np.int64)
    by = np.floor(u[:, a2] * (2**j2)).astype(np.int64)
    return np.bincount(bx * (2**j2) + by, minlength=2**(j1 + j2))


def check_balance(u, tag):
    ok = True
    interior = bool((u > 0.0).all() and (u < 1.0).all())
    ok = ok and interior
    # 1-D exact dyadic counts
    d1_ok = True
    for d in range(3):
        for j in range(1, M + 1):
            cnt = np.bincount(np.floor(u[:, d] * (2**j)).astype(int),
                              minlength=2**j)
            if not np.all(cnt == N // (2**j)):
                d1_ok = False
    ok = ok and d1_ok
    # (dim1,dim2) 2-D elementary boxes: one per box at every (j1,j2), j1+j2=M
    net01_ok = True
    for j1 in range(0, M + 1):
        if not np.all(box_hist(u, 0, 1, j1, M - j1) == 1):
            net01_ok = False
    ok = ok and net01_ok
    # preservation vs the unscrambled net for all three projections
    preserve_ok = True
    for (a1, a2) in [(0, 1), (0, 2), (1, 2)]:
        for j1 in range(0, M + 1):
            j2 = M - j1
            h_raw = np.sort(box_hist(u_raw, a1, a2, j1, j2))
            h_scr = np.sort(box_hist(u,     a1, a2, j1, j2))
            if not np.array_equal(h_raw, h_scr):
                preserve_ok = False
    ok = ok and preserve_ok
    print(f"  [{tag}] interior(0,1)={interior}  1-D exact={d1_ok}  "
          f"(d1,d2) (0,m,2)-net={net01_ok}  occupancy-preserved={preserve_ok}"
          f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


print("\n=== Q0(b,c) Owen-scramble dyadic balance / net preservation ===")
seeds = [12345, 777, 2024]
usets = {}
b_ok = True
for s in seeds:
    u = np.loadtxt(f"tests/qmc/q0_scr_{s}.txt")[:, 1:4]
    usets[s] = u
    b_ok = b_ok and check_balance(u, f"seed {s}")

diff01 = float(np.abs(usets[seeds[0]] - usets[seeds[1]]).max())
diff02 = float(np.abs(usets[seeds[0]] - usets[seeds[2]]).max())
print(f"  replicate max |du|: seed{seeds[0]} vs {seeds[1]} = {diff01:.3f}, "
      f"vs {seeds[2]} = {diff02:.3f} (should be O(1))")
c_ok = diff01 > 0.1 and diff02 > 0.1

# report which Sobol 2-D projections are NOT (0,2)-nets (informational)
print("\n  informational: raw Sobol 2-D (0,m,2)-net status")
for (a1, a2) in [(0, 1), (0, 2), (1, 2)]:
    onep = all(np.all(box_hist(u_raw, a1, a2, j1, M - j1) == 1)
               for j1 in range(0, M + 1))
    print(f"    dims({a1},{a2}): one-point-per-elementary-box = {onep}")

print("\nQ0 GATE:", "PASS" if (a_ok and b_ok and c_ok) else "FAIL")
