#!/usr/bin/env python3
"""Absorption ledger by radius (2026-10-17): the H0 absorptions of each
source (Ledger_H0 of a rates file written with par%absorption_ledger),
summed in radius bins, each divided by the H recombination demand of the
bin, n_e n(H+) alpha_A(T) V (tools/python/mochii_rates.py).  The column
'case B' is everything but channel 1 (the H ground captures), the column
'demand B' n_e n(H+) alpha_B V.  Checks that the leaf sums reproduce the
LEDGER: totals of the log.  Usage: ledger_radial.py <rates.h5> ..."""
import sys, numpy as np, h5py
sys.path.insert(0, '../../../../tools/python'); import mochii_rates as mr
Q, nH = 5.416e47, 3000.0
SRC = ['star', 'dif1', 'dif2', 'dif3', 'dif4', 'dif5', 'dif6', 'dif7', 'dif8', 'FUV']
lgT = np.linspace(3.0, 5.0, 401)
aA = np.array([mr.alphaA_HII(10**x) for x in lgT]); aB = np.array([mr.alphaB_HII(10**x) for x in lgT])
edges = np.array([0.03, 0.045, 0.06, 0.07, 0.08, 0.085, 0.09, 0.10, 0.11, 0.12, 0.13, 0.16])
for f in sys.argv[1:]:
    with h5py.File(f, 'r') as h:
        g = lambda k: h[k]['data'][()]
        d = h['Gamma_HI'].attrs['DIST_CM'][0]
        L = g('Ledger_H0')
        if L.shape[0] != len(SRC): L = L.T                     # (source, leaf)
        V, ne = (g('LeafSize')*d)**3, g('n_e')
        r = np.sqrt((g('LeafXYZ')**2).sum(axis=0))*d/3.0857e18
        T, xp = g('T_e'), g('x_HII')
    m = ne > 0
    lt = np.log10(np.maximum(T, 1.0))
    dA = np.where(m, ne*nH*xp*np.interp(lt, lgT, aA)*V, 0.0)
    dB = np.where(m, ne*nH*xp*np.interp(lt, lgT, aB)*V, 0.0)
    print(f, ' totals per Q: ' + ' '.join(f'{s} {L[i].sum()/Q:.4f}' for i, s in enumerate(SRC) if L[i].sum() > 0))
    print(f"{'r [pc]':13s}{'dem A/Q':>9s}{'dem B/Q':>9s}" + ''.join(f'{s:>8s}' for s in SRC[:9]) + f"{'case B':>8s}")
    b = np.digitize(r, edges)
    for i in range(1, len(edges)):
        k = (b == i) & m
        DA, DB = dA[k].sum(), dB[k].sum()
        if DA <= 0: continue
        fr = [L[j][k].sum()/DB for j in range(9)]
        print(f'{edges[i-1]:.3f}-{edges[i]:.3f}{DA/Q:9.4f}{DB/Q:9.4f}' + ''.join(f'{x:8.4f}' for x in fr)
              + f'{(sum(fr) - fr[1]):8.4f}')
    print('  (source columns: H0 absorptions of the bin over its case-B demand)')
