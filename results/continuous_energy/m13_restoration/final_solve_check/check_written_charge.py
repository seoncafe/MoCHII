#!/usr/bin/env python3
"""Charge neutrality of a written state: n_e written against the charge the
written H, He and metal stages carry, n_H [x_HII + y (x_HeII + 2 x_HeIII)
+ sum_el A_el sum_j j x_j], with the uniform n_H = 100 cm^-3 of the HII40
grid (tests/g2_hii/hii40_L6_rout5p2.fits).

usage: check_written_charge.py RATES.h5 [...]
"""
import sys
import numpy as np
import h5py

ABUND = dict(c=2.2e-4, n=4.0e-5, o=3.3e-4, ne=5.0e-5, s=9.0e-6)
Y_HE = 0.1
N_H = 100.0

for path in sys.argv[1:]:
    with h5py.File(path, 'r') as h:
        ne = h['n_e/data'][:]
        q = h['x_HII/data'][:] + Y_HE*(h['x_HeII/data'][:] + 2*h['x_HeIII/data'][:])
        for el, a in ABUND.items():
            x = h[f'x_{el}_stages/data'][:]
            q += a*(x*np.arange(x.shape[1])).sum(axis=1)
        V = h['LeafSize/data'][:]**3
    g = ne > 0
    nH = ne[g]/q[g]
    r = nH/N_H - 1
    print(f'{path}: gas leaves {g.sum()}, '
          f'max |n_e/charge - 1| {np.abs(r).max():.3e}, '
          f'volume-mean {np.average(np.abs(r), weights=V[g]):.3e}, '
          f'leaves > 1e-6: {(np.abs(r) > 1e-6).sum()}')
