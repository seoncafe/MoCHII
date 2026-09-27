#!/usr/bin/env python3
"""Check the par%restart_file runs written by run_check.sh in directory argv[1]."""
import sys
import h5py
import numpy as np

W = sys.argv[1]
KEYS = ('x_HI', 'x_HII', 'x_HeI', 'x_HeII', 'x_HeIII', 'n_e', 'T_e')
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
      for k in ('x_HI', 'x_HII', 'x_HeI', 'x_HeII', 'x_HeIII', 'T_e')}
print(f' (b) 32^3 grid: {use.sum()} gas leaves take the containing previous leaf:', eq)
ok &= all(eq.values()) and use.sum() > 0

#--- (c) files written before 2026-09-26, without x_HII, without x_HeIII,
#--- and without both: the missing fraction from the complement, everything
#--- else as in (a), and the run says so
X3 = np.maximum(0.0, 1.0 - (S['x_HeI'] + S['x_HeII']))  # as previous_state_mod forms it
for c, missing in (('noxhii', ('x_HII',)), ('noxheiii', ('x_HeIII',)),
                   ('noboth', ('x_HII', 'x_HeIII'))):
    lg = h5py.File(f'{W}/{c}_rates.h5', 'r')
    log = open(f'{W}/{c}.log').read()
    eq = {k: bool(np.array_equal(S[k], lg[k + '/data'][:]))
          for k in KEYS if k not in missing}
    if 'x_HII' in missing:
        eq['x_HII = 1 - x_HI'] = bool(np.array_equal(lg['x_HII/data'][:],
                                                     1.0 - S['x_HI']))
        eq['x_HII note'] = 'carries no x_HII' in log
    if 'x_HeIII' in missing:
        eq['x_HeIII = 1 - x_HeI - x_HeII'] = bool(
            np.array_equal(lg['x_HeIII/data'][:], X3))
        eq['x_HeIII note'] = 'carries no x_HeIII' in log
    print(f' (c) restart from a file without {" and ".join(missing)}:', eq)
    ok &= all(eq.values())

print(' GATE (restart_file):', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
