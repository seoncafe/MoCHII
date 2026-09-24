#!/usr/bin/env python3
"""Secondary-ionization comparison: off (current binary), on (current binary,
SvS85 x DYL99), on (pre-change binary, SvS85 alone from 40 eV)."""
import sys, re, numpy as np, h5py
d = sys.argv[1] if len(sys.argv) > 1 else '.'
def load(tag):
    h = h5py.File(f'{d}/sec_{tag}_rates.h5', 'r')
    g = lambda k: h[k + '/data'][:]
    r = dict(T=g('T_e'), ne=g('n_e'), xHI=g('x_HI'), xHeI=g('x_HeI'),
             xHeII=g('x_HeII'), V=g('LeafSize')**3 if g('LeafSize').ndim == 1
             else np.prod(g('LeafSize'), axis=-1), xyz=g('LeafXYZ'))
    if r['xyz'].shape[0] == 3 and r['xyz'].ndim == 2: r['xyz'] = r['xyz'].T
    r['rad'] = np.sqrt((r['xyz']**2).sum(-1))
    lines = {}
    for ln in open(f'{d}/sec_{tag}_lines.txt'):
        m = re.search(r'L\(Hbeta\)\s*=\s*([0-9.eE+-]+)', ln)
        if m: lines['Hb'] = float(m.group(1))
    r['lines'] = lines
    return r
C = {t: load(t) for t in ('off_new', 'on_new', 'on_old')}
nH = 1.0
for t, r in C.items():
    V = r['V']; np_ = (1 - r['xHI']); w = V*np_*r['ne']
    TNpNe = np.sum(r['T']*w)/np.sum(w)
    xHI_v = np.sum(V*r['xHI'])/V.sum()
    xHeI_v = np.sum(V*r['xHeI'])/V.sum()
    xHeII_v = np.sum(V*r['xHeII'])/V.sum()
    ne_v = np.sum(V*r['ne'])/V.sum()
    neutral = r['xHI'] > 0.5
    print(f'{t:8s} <T[NpNe]>={TNpNe:9.2f} K  <xHI>_V={xHI_v:.5e}  <xHeI>_V={xHeI_v:.5e}'
          f'  <xHeII>_V={xHeII_v:.5e}  <ne>_V={ne_v:.5e}  L(Hb)={r["lines"].get("Hb",np.nan):.5e}'
          f'  V(xHI>0.5)={V[neutral].sum()/V.sum():.5f}')
