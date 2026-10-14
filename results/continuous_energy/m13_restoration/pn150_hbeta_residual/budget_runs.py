#!/usr/bin/env python3
"""R_B(H), R_B(He+), R_B(He2+) per stellar Q(H) = 5.416e47 s^-1 and L(Hbeta)
for MoCHII PN150 runs (base names on the command line), all on the
alpha_B of tools/python/mochii_rates.py."""
import sys, numpy as np, h5py
sys.path.insert(0, '../../../../tools/python'); import mochii_rates as mr
Q = 5.416e47; nH = 3000.0
def g(h, k): return h[k]['data'][()] if isinstance(h[k], h5py.Group) else h[k][()]
def lhb(p):
    for r in open(p):
        if r.startswith('# L(Hbeta)'): return float(r.split('=')[1].split()[0])
for base in sys.argv[1:]:
    with h5py.File(base + '_rates.h5', 'r') as h:
        d = h['Gamma_HI'].attrs['DIST_CM'][0]; V = (g(h,'LeafSize')*d)**3
        ne = g(h,'n_e'); m = ne > 0
        T = g(h,'T_e')[m]; ne = ne[m]; V = V[m]
        xp = g(h,'x_HII')[m]; x2 = g(h,'x_HeII')[m]; x3 = g(h,'x_HeIII')[m]
    f = lambda fn: np.array([fn(t) for t in T])
    RH = np.sum(ne*nH*xp*f(mr.alphaB_HII)*V)
    R2 = np.sum(ne*nH*0.1*x2*f(mr.alphaB_HeII)*V)
    R3 = np.sum(ne*nH*0.1*x3*f(mr.alphaB_HeIII)*V)
    Tm = np.sum(ne*nH*xp*T*V)/np.sum(ne*nH*xp*V)
    L = lhb(base + '_lines.txt')
    print(f'{base:22s} R_B(H)/Q {RH/Q:.4f}  R_B(He+)/Q {R2/Q:.4f}  R_B(He2+)/Q {R3/Q:.4f}'
          f'  <T> {Tm:7.0f} K  L(Hb) {L:.4e}')
