#!/usr/bin/env python3
"""Local change between the last written iterate (<base>_rates_iterate.h5,
one iteration before the end) and the written state (<base>_rates.h5) of
PN150 runs, for base names on the command line: counts of leaves whose T_e
moved by more than 10% / 1%, their share of R_B(H) (on the alpha_B of
tools/python/mochii_rates.py), the x_HII and T_e range of the >10% leaves,
and the leaves above the 30000 K edge of the SH95 line tables."""
import sys, numpy as np, h5py
sys.path.insert(0, '../../../../tools/python'); import mochii_rates as mr
nH = 3000.0
def g(h, k): return h[k]['data'][()].ravel()
lgT = np.linspace(3.0, 5.0, 201); lga = np.array([mr.alphaB_HII(10**t) for t in lgT])
for base in sys.argv[1:]:
    with h5py.File(base + '_rates.h5', 'r') as a, h5py.File(base + '_rates_iterate.h5', 'r') as b:
        d = a['Gamma_HI'].attrs['DIST_CM'][0]; V = (g(a, 'LeafSize')*d)**3
        T1, T0, x1, x0, ne = g(a, 'T_e'), g(b, 'T_e'), g(a, 'x_HII'), g(b, 'x_HII'), g(a, 'n_e')
    m = ne > 0
    w = np.where(m, ne*nH*x1*np.interp(np.log10(np.maximum(T1, 1.0)), lgT, lga)*V, 0.0); W = w.sum()
    dT = np.abs(T1 - T0)/np.maximum(T0, 1.0)
    big, mid, hot = m & (dT > 0.1), m & (dT > 0.01), m & (T1 > 3.0e4)
    print(f'{base}: leaves {m.sum()}')
    print(f'  |dT|/T > 10%: {big.sum()} leaves, R_B(H) share {w[big].sum()/W:.2e}', end='')
    if big.any():
        print(f', x_HII {x0[big].min():.1e}-{x0[big].max():.1e}, T_e {T0[big].min():.0f}-{T0[big].max():.0f} K')
    else:
        print()
    print(f'  |dT|/T > 1%: {mid.sum()} leaves, R_B(H) share {w[mid].sum()/W:.2e}')
    print(f'  T_e > 30000 K: {hot.sum()} leaves, R_B(H) share {w[hot].sum()/W:.2e}, max T_e {T1[m].max():.0f} K')
