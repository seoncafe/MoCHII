#!/usr/bin/env python3
"""Cloudy c25.00's converged PN150 structure as a MoCHII restart file, for the
fixed-structure test (2026-10-17): every gas leaf of a MoCHII L6 rates file
takes the Cloudy zone values at its radius, averaged over 4^3 points inside
the leaf (T_e, n_e, x_HI, x_HII, x_HeI, x_HeII, x_HeIII from gam.ovr; the
metal photoionization rates Gamma_<el>_stages from gam.gam, summed over the
shells of each ion, behind the stage fractions MoCHII then forms).  Zones
are spherical shells from 1e17 cm, values linear in radius between zone
centres; beyond the last zone (0.135 pc, Cloudy's 4000 K stop) the last
zone's values are held.  Phi_HeI_23S and the other datasets keep the
template's values.  Usage: make_cloudy_state.py <template rates.h5> <out.h5>"""
import shutil, sys, numpy as np, h5py
EL = dict(c=6, n=7, o=8, ne=10, mg=12, si=14, s=16)
o = np.loadtxt('gam.ovr')
rz = 1e17 + o[:, 0]
cols = dict(T_e=1, n_e=4, x_HI=6, x_HII=7, x_HeI=8, x_HeII=9, x_HeIII=10)
# gam.gam: zone blocks, each starting at Z=1, ion=1, shell=1
blocks, cur = [], None
for line in open('gam.gam'):
    if line.startswith('#'):
        continue
    w = line.split()
    Z, ion = int(w[0]), int(w[1])
    if (Z, ion, int(w[2])) == (1, 1, 1):
        cur = {}; blocks.append(cur)
    cur[(Z, ion)] = cur.get((Z, ion), 0.0) + float(w[3])
nz = min(len(blocks), len(rz))
rz = rz[:nz]

def profile(vals, r):
    return np.interp(r, rz, vals[:nz])

shutil.copyfile(sys.argv[1], sys.argv[2])
with h5py.File(sys.argv[2], 'r+') as h:
    g = lambda k: h[k]['data']
    dcm = h['Gamma_HI'].attrs['DIST_CM'][0]
    xyz, size, ne = g('LeafXYZ')[()], g('LeafSize')[()], g('n_e')[()]
    gas = np.where(ne > 0)[0]
    off = (np.arange(4) + 0.5)/4.0 - 0.5
    ox, oy, oz = np.meshgrid(off, off, off, indexing='ij')
    sub = np.stack([ox.ravel(), oy.ravel(), oz.ravel()])         # (3, 64)
    pts = xyz[:, gas, None] + size[None, gas, None]*sub[:, None, :]
    r = np.sqrt((pts**2).sum(axis=0))*dcm                         # (ngas, 64) cm
    for k, c in cols.items():
        d = g(k)[()]
        d[gas] = profile(o[:, c], r).mean(axis=1)
        g(k)[...] = d
    for el, Z in EL.items():
        G = g(f'Gamma_{el}_stages')[()]
        for it in range(1, G.shape[1] + 1):
            vals = np.array([b.get((Z, it), 0.0) for b in blocks[:nz]])
            G[gas, it - 1] = profile(vals, r).mean(axis=1)
        g(f'Gamma_{el}_stages')[...] = G
    print(f'{len(gas)} gas leaves from {nz} Cloudy zones; r(leaf centres) '
          f'{np.sqrt((xyz[:, gas]**2).sum(axis=0)).min()*dcm:.4e}-'
          f'{np.sqrt((xyz[:, gas]**2).sum(axis=0)).max()*dcm:.4e} cm, '
          f'zones {rz[0]:.4e}-{rz[-1]:.4e} cm')
