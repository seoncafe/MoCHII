#!/usr/bin/env python3
"""Summaries and radial profiles of the on-the-spot He III runs.

usage: analyze.py BASE [BASE ...]
For each <BASE>_rates.h5 / <BASE>_lines.txt: <T[NpNe]> (n_p n_e weighted),
L(Hbeta), He II 4686 (absolute and relative to Hbeta), and the radii where
the spherically averaged x(H0) and x(He2+) cross 0.5 (volume-weighted shells).
"""
import sys
import numpy as np
import h5py

PC = 3.0856775814913673e18
AHE = 0.1


def load(base):
    f = h5py.File(base + '_rates.h5', 'r')
    d = {k: f[k]['data'][()] for k in ('x_HI', 'x_HeI', 'x_HeII', 'n_e', 'T_e',
                               'LeafXYZ', 'LeafSize')}
    f.close()
    xyz = d['LeafXYZ']
    if xyz.shape[1] != 3:
        xyz = xyz.T
    d['r'] = np.sqrt((xyz**2).sum(axis=1))
    size = np.ravel(d['LeafSize'])
    d['vol'] = size**3
    d['xHeIII'] = np.clip(1.0 - d['x_HeI'] - d['x_HeII'], 0.0, 1.0)
    return d


def lines(base):
    out = {}
    for ln in open(base + '_lines.txt'):
        if ln.startswith('# L(Hbeta)'):
            out['Hbeta'] = float(ln.split('=')[1].split()[0])
        s = ln.split()
        if len(s) >= 6 and s[-1] == 'HeII4687':
            out['4686'] = float(s[3])
    return out


def profile(d, key, edges):
    r, v, x = d['r'], d['vol'], d[key]
    i = np.digitize(r, edges) - 1
    ok = (i >= 0) & (i < len(edges) - 1) & (d['n_e'] > 0)
    num = np.bincount(i[ok], weights=(v*x)[ok], minlength=len(edges) - 1)
    den = np.bincount(i[ok], weights=v[ok], minlength=len(edges) - 1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return num/den


def crossing(rc, y, level=0.5):
    """first radius where the profile rises through level (linear)"""
    for k in range(len(y) - 1):
        if np.isfinite(y[k]) and np.isfinite(y[k+1]) and y[k] < level <= y[k+1]:
            return rc[k] + (level - y[k])*(rc[k+1] - rc[k])/(y[k+1] - y[k])
    return np.nan


def summary(base, nH=100.0):
    d = load(base)
    L = lines(base)
    # n_p n_e weight; nH uniform in these grids is not stored, it cancels
    w = d['vol']*(1.0 - d['x_HI'])*d['n_e']
    tnpne = (w*d['T_e']).sum()/w.sum()
    rmax = d['r'].max()
    edges = np.linspace(0.0, rmax, 161)
    rc = 0.5*(edges[1:] + edges[:-1])
    xh0 = profile(d, 'x_HI', edges)
    x3 = profile(d, 'xHeIII', edges)
    rH = crossing(rc, xh0)
    r3 = crossing(rc, 1.0 - x3)
    return dict(T=tnpne, Hb=L.get('Hbeta'), He4686=L.get('4686'), rH=rH, rHe3=r3)


if __name__ == '__main__':
    for b in sys.argv[1:]:
        s = summary(b)
        print(f"{b:40s} <T[NpNe]> = {s['T']:9.2f} K  L(Hb) = {s['Hb']:.5e}  "
              f"L(4686) = {s['He4686']:.5e} ({s['He4686']/s['Hb']:.5f} Hb)  "
              f"r(xH0=0.5) = {s['rH']:.4f}  r(xHe2+=0.5) = {s['rHe3']:.4f}")
