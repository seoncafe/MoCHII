#!/usr/bin/env python3
"""Photoionization rates (n_H V weighted sums, state independent) of one
element stage in radial shells, run B over run A - 1.
gamma_shells.py <set> <A> <B> <element> <stage index> [shell edges]"""
import sys, h5py, numpy as np
W = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib/runs/' + sys.argv[1]
A, B, el, j = sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
sh = [float(v) for v in sys.argv[6:]] or [0, 1.8, 2.2, 2.9, 3.6, 4.1]
ra = h5py.File(f'{W}/{A}/bb1e5_rates.h5', 'r'); rb = h5py.File(f'{W}/{B}/bb1e5_rates.h5', 'r')
xyz = ra['LeafXYZ/data'][:]; rad = np.sqrt((xyz**2).sum(0)); V = ra['LeafSize/data'][:]**3
nH = h5py.File(f'{W}/{A}/bb1e5_emis.h5', 'r')['n_H/data'][:]
key = 'Gamma_HeII/data' if el == 'he' else f'Gamma_{el}_stages/data'
ga = ra[key][:] if el == 'he' else ra[key][:, j]; gb = rb[key][:] if el == 'he' else rb[key][:, j]
out = []
for i in range(len(sh)-1):
    m = (rad >= sh[i]) & (rad < sh[i+1]) & (nH > 0)
    out.append(f'{sh[i]}-{sh[i+1]}: {np.sum(V[m]*nH[m]*gb[m])/np.sum(V[m]*nH[m]*ga[m])-1:+.2e}')
print(f'{el} {j} ({B}/{A}-1):', '  '.join(out))
