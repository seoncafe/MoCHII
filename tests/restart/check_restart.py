#!/usr/bin/env python3
"""Check the par%restart_file runs written by run_check.sh in directory argv[1]."""
import sys
import h5py
import numpy as np

W = sys.argv[1]
KEYS = ('x_HI', 'x_HeI', 'x_HeII', 'n_e', 'T_e')
src = h5py.File(f'{W}/fresh_rates.h5', 'r')
S = {k: src[k + '/data'][:] for k in KEYS + ('LeafXYZ', 'LeafSize')}
ok = True

same = h5py.File(f'{W}/same_rates.h5', 'r')
eq = {k: bool(np.array_equal(S[k], same[k + '/data'][:])) for k in KEYS}
print(' (a) same grid, gas_niter = 0, state returned bit for bit:', eq)
ok &= all(eq.values())

co = h5py.File(f'{W}/coarse_rates.h5', 'r')
xyz = co['LeafXYZ/data'][:]
oxyz, osz = S['LeafXYZ'], S['LeafSize']
s = osz[0]
org = (oxyz - osz / 2).min(axis=1)
n = int(round((oxyz.max(axis=1) - org + s / 2)[0] / s))
lut = np.full((n, n, n), -1)
oi = np.floor((oxyz.T - org) / s).astype(int)
lut[oi[:, 0], oi[:, 1], oi[:, 2]] = np.arange(oxyz.shape[1])
idx = np.floor((xyz.T - org) / s).astype(int)
inside = np.all((idx >= 0) & (idx < n), axis=1)
j = np.full(xyz.shape[1], -1)
j[inside] = lut[idx[inside, 0], idx[inside, 1], idx[inside, 2]]
use = (j >= 0) & (S['n_e'][np.maximum(j, 0)] > 0) & (co['n_e/data'][:] > 0)
eq = {k: bool(np.array_equal(co[k + '/data'][:][use], S[k][j[use]]))
      for k in ('x_HI', 'x_HeI', 'x_HeII', 'T_e')}
print(f' (b) 32^3 grid: {use.sum()} gas leaves take the containing previous leaf:', eq)
ok &= all(eq.values()) and use.sum() > 0

print(' GATE (restart_file):', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
