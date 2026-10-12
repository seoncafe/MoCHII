#!/usr/bin/env python3
"""J_nu by bin in radial shells, ratio of run B over run A (sums over leaves,
volume weighted); jshell.py <set> <A> <B> [r_lo r_hi ...]"""
import sys, h5py, numpy as np
W = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib/runs/' + sys.argv[1]
A, B = sys.argv[2], sys.argv[3]
sh = [float(v) for v in sys.argv[4:]] or [0, 1.7, 2.9, 4.0, 5.3]
ra = h5py.File(f'{W}/{A}/bb1e5_rates.h5', 'r'); rb = h5py.File(f'{W}/{B}/bb1e5_rates.h5', 'r')
xyz = ra['LeafXYZ/data'][:]; rad = np.sqrt((xyz**2).sum(0)); V = ra['LeafSize/data'][:]**3
E = ra['E_bin/data'][:]; dE = ra['dE_bin/data'][:]
JA = ra['J_nu/data'][:]; JB = rb['J_nu/data'][:]
gas = h5py.File(f'{W}/{A}/bb1e5_emis.h5', 'r')['n_H/data'][:] > 0
sel = np.where((E > 13) & (E < 60))[0]
hdr = '  E_lo-E_hi [eV] ' + ' '.join(f'  r {sh[i]:.1f}-{sh[i+1]:.1f}' for i in range(len(sh)-1))
print(hdr)
for ib in sel:
    row = f'{E[ib]-dE[ib]/2:7.2f}-{E[ib]+dE[ib]/2:7.2f} '
    for i in range(len(sh)-1):
        m = gas & (rad >= sh[i]) & (rad < sh[i+1])
        row += f' {np.sum(V[m]*JB[m, ib])/np.sum(V[m]*JA[m, ib])-1:+.3e}'
    print(row)
# 4 pi J_nu dnu in the S+ windows, absolute, for the shells
for lo, hi in [(23.34, 24.59), (24.59, 34.86)]:
    b = np.where((E - dE/2 >= lo - 1e-3) & (E + dE/2 <= hi + 1e-3))[0]
    row = f'band {lo}-{hi}: '
    for i in range(len(sh)-1):
        m = gas & (rad >= sh[i]) & (rad < sh[i+1])
        ja = np.sum(V[m, None]*JA[m][:, b]*dE[b]); jb = np.sum(V[m, None]*JB[m][:, b]*dE[b])
        row += f' {jb/ja-1:+.3e}'
    print(row)
