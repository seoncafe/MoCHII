#!/usr/bin/env python3
"""Change of the S+ -> S2+ photoionization rate between two runs, split by
J_nu bin, summed over leaves with weight V n_H (state independent, so the
ratio carries no selection by either run's noise), in radial shells.
sigma(S+) is the VFKY96 outer-shell fit of data/atomic/element_s.txt
(transition 2), evaluated at the bin center (an estimate: the code uses
the exact packet energy).  gamma_bins.py <set> <A> <B> [shell edges]"""
import sys, h5py, numpy as np
W = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib/runs/' + sys.argv[1]
A, B = sys.argv[2], sys.argv[3]
sh = [float(v) for v in sys.argv[4:]] or [0, 2.9, 4.0, 5.3]
def vfky(E, Eth, E0, s0, ya, P, yw, y0, y1):
    x = E/E0 - y0; y = np.sqrt(x*x + y1*y1)
    F = ((x - 1)**2 + yw**2)*y**(0.5*P - 5.5)*(1 + np.sqrt(y/ya))**(-P)
    return np.where(E >= Eth, s0*F, 0.0)
PH = [2.333788042842e+01, 8.787, 3.136e+02, 3.442, 1.281e+01, 7.354e-01, 2.782, 1.788e-01]
ra = h5py.File(f'{W}/{A}/bb1e5_rates.h5', 'r'); rb = h5py.File(f'{W}/{B}/bb1e5_rates.h5', 'r')
ea = h5py.File(f'{W}/{A}/bb1e5_emis.h5', 'r')
xyz = ra['LeafXYZ/data'][:]; rad = np.sqrt((xyz**2).sum(0)); V = ra['LeafSize/data'][:]**3
E = ra['E_bin/data'][:]; dE = ra['dE_bin/data'][:]
sig = vfky(E, *PH)
w = V*ea['n_H/data'][:]
JA = ra['J_nu/data'][:]; JB = rb['J_nu/data'][:]
GA = ra['Gamma_s_stages/data'][:, 1]; GB = rb['Gamma_s_stages/data'][:, 1]
for i in range(len(sh)-1):
    m = (rad >= sh[i]) & (rad < sh[i+1]) & (w > 0)
    ga = np.sum(w[m, None]*JA[m]*(sig*dE/E)[None, :], 0)
    gb = np.sum(w[m, None]*JB[m]*(sig*dE/E)[None, :], 0)
    tot = ga.sum()
    print(f'shell {sh[i]}-{sh[i+1]} pc: n_H V share {w[m].sum()/w[w>0].sum():.3f};'
          f' code Gamma(S+) sum ratio B/A-1 {np.sum(w[m]*GB[m])/np.sum(w[m]*GA[m])-1:+.3e};'
          f' J-based {gb.sum()/tot-1:+.3e}')
    o = np.argsort(-abs(gb - ga))[:6]
    for ib in o:
        print(f'    bin {E[ib]-dE[ib]/2:6.2f}-{E[ib]+dE[ib]/2:6.2f} eV: share of Gamma(A) {ga[ib]/tot:.3f}, contribution to dGamma/Gamma {(gb[ib]-ga[ib])/tot:+.3e}')
