#!/usr/bin/env python3
"""Radial profiles and front radii of runs in one set.  fronts.py <set> <run> [<run> ...]"""
import sys, h5py, numpy as np
W = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib/runs/' + sys.argv[1]
def st(x):
    r = h5py.File(f'{W}/{x}/bb1e5_rates.h5', 'r'); e = h5py.File(f'{W}/{x}/bb1e5_emis.h5', 'r')
    xyz = r['LeafXYZ/data'][:]; V = r['LeafSize/data'][:]**3
    rad = np.sqrt((xyz**2).sum(0))
    d = dict(rad=rad, V=V, nH=e['n_H/data'][:], xHII=r['x_HII/data'][:], xHeI=r['x_HeI/data'][:],
             xHeIII=r['x_HeIII/data'][:], xS2=e['x_s_stages/data'][:, 2], xS1=e['x_s_stages/data'][:, 1],
             e9532=e['emis_s_3/data'][:, 5], T=r['T_e/data'][:])
    return d
def front(d, key, thr, rising):
    # volume-weighted shell means on a fine radial grid, first crossing
    edges = np.linspace(0, d['rad'].max(), 400)
    idx = np.digitize(d['rad'], edges)
    prof = np.array([np.sum(d['V'][idx == i]*d[key][idx == i])/max(np.sum(d['V'][idx == i]), 1e-300) for i in range(1, len(edges))])
    rc = 0.5*(edges[1:] + edges[:-1])
    ok = np.sum(d['V'][:, None]*(idx[:, None] == np.arange(1, len(edges))[None, :]), 0) > 0
    rc, prof = rc[ok], prof[ok]
    c = np.where((prof < thr) if not rising else (prof > thr))[0]
    if len(c) == 0: return np.nan
    i = c[0]
    if i == 0: return rc[0]
    return rc[i-1] + (rc[i]-rc[i-1])*(thr-prof[i-1])/(prof[i]-prof[i-1])
for x in sys.argv[2:]:
    d = st(x)
    gas = d['nH'] > 0
    print(f'{x:8s} R(x_HII=0.5) {front(d,"xHII",0.5,False):.5f}  R(x_HeI=0.5) {front(d,"xHeI",0.5,True):.5f}'
          f'  R(x_HeIII=0.5) {front(d,"xHeIII",0.5,False):.5f}  sum V nH xS2 {np.sum(d["V"]*d["nH"]*d["xS2"]):.6e}'
          f'  sum V e9532 {np.sum(d["V"]*d["e9532"]):.6e}')
