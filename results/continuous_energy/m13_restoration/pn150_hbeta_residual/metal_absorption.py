#!/usr/bin/env python3
"""Photoionizations per second by every absorber in the MoCHII PN150 state:
H, He0, He+ (Gamma_* from the rates file) and each metal stage
(Gamma_<el>_stages x n(X^i)), against the stellar Q(H) = 5.416e47 s^-1."""
import numpy as np, h5py
def g(h, k): return h[k]['data'][()] if isinstance(h[k], h5py.Group) else h[k][()]
Q = 5.416e47; nH = 3000.0
with h5py.File('../pn150_p27_nodrift_rates.h5', 'r') as h:
    d = h['Gamma_HI'].attrs['DIST_CM'][0]
    V = (g(h, 'LeafSize')*d)**3
    ne = g(h, 'n_e'); gas = ne > 0
    tot = {}
    tot['H0'] = np.sum(g(h,'Gamma_HI')[gas]*nH*g(h,'x_HI')[gas]*V[gas])
    tot['He0'] = np.sum(g(h,'Gamma_HeI')[gas]*nH*0.1*g(h,'x_HeI')[gas]*V[gas])
    tot['He+'] = np.sum(g(h,'Gamma_HeII')[gas]*nH*0.1*g(h,'x_HeII')[gas]*V[gas])
    met = {}
    for k in h:
        if k.startswith('x_') and k.endswith('_stages'):
            el = k[2:-7]
            x = g(h, k); G = g(h, 'Gamma_'+el+'_stages')
            if x.shape[0] != G.shape[0] + 1: x = x.T; G = G.T
            ab = h[k].attrs['ABUND'][0] if 'ABUND' in h[k].attrs else h[k]['data'].attrs['ABUND'][0]
            per = [np.sum(G[i,gas]*ab*nH*x[i,gas]*V[gas]) for i in range(G.shape[0])]
            met[el] = per
for k, v in tot.items(): print(f'{k:6s} {v:11.4e}  /Q = {v/Q:.4f}')
mt = 0
for el, per in met.items():
    s = sum(per); mt += s
    print(f'{el:6s} {s:11.4e}  /Q = {s/Q:.4f}   by stage ' + ' '.join(f'{p/Q:.4f}' for p in per))
print(f'metals {mt:11.4e}  /Q = {mt/Q:.4f}')
