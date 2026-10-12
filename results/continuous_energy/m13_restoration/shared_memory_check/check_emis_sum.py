#!/usr/bin/env python3
"""Leaf emissivities of the metal lines (emis_<el>_<stage>, summed onto rank 0
in chunks since 2026-10-12) against the line luminosities of _lines.txt:
sum(emis V) must equal L for every written line.  Run: python3 check_emis_sum.py emis"""
import sys, h5py, numpy as np
base = sys.argv[1]
L = {}
for row in open(base + '_lines.txt'):
    if row.startswith('#'): continue
    f = row.split()
    L[(f[0], int(f[1]), round(float(f[2]), 2))] = float(f[3])
with h5py.File(base + '_emis.h5', 'r') as h:
    get = lambda k: h[k]['data'][()] if isinstance(h[k], h5py.Group) else h[k][()]
    with h5py.File(base + '_rates.h5', 'r') as r:
        dist = r['Gamma_HI'].attrs['DIST_CM'][0]
        size = (r['LeafSize']['data'][()] if isinstance(r['LeafSize'], h5py.Group) else r['LeafSize'][()])
    vol = (size*dist)**3
    worst, nlines, nsets = 0.0, 0, 0
    for k in h:
        if not k.startswith('emis_') or k.count('_') != 2: continue
        _, el, st = k.split('_')
        em = get(k); wl = get('wl_' + el + '_' + st)
        em = em.reshape(len(wl), -1) if em.shape[0] == len(wl) else em.T
        s = em @ vol
        for w, v in zip(wl, s):
            ref = L[(el, int(st), round(float(w), 2))]
            worst = max(worst, abs(v/ref - 1)); nlines += 1
        nsets += 1
print(f'{nsets} metal stages, {nlines} lines: max |sum(emis V)/L - 1| = {worst:.3e}')
