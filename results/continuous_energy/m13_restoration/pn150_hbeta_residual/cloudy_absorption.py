#!/usr/bin/env python3
"""Photoionizations per second by every absorber in Cloudy c25.00's PN150
(gam.in = pn150_c25.in + save gammas / element files): sum over zones of
PhotoRate (all shells) x n(X^i) x dV, against Q(H) = 5.416e47 s^-1, for the
comparison with metal_absorption.py (MoCHII)."""
import numpy as np
Q = 5.416e47; nH = 3000.0
Z = {'h':1,'he':2,'c':6,'n':7,'o':8,'ne':10,'mg':12,'si':14,'s':16}
ab = {'h':1.0,'he':0.1,'c':10**-3.523,'n':1e-4,'o':10**-3.222,'ne':10**-3.824,
      'mg':10**-4.523,'si':10**-4.523,'s':10**-4.824}
rows = [l.split() for l in open('gam.gam') if not l.startswith('#')]
# zone blocks start at '1 1 1'
blocks = []; cur = None
for r in rows:
    k = (int(r[0]), int(r[1]), int(r[2]))
    if k == (1,1,1):
        cur = {}; blocks.append(cur)
    cur[(k[0], k[1])] = cur.get((k[0], k[1]), 0.0) + float(r[3])
dep = np.loadtxt('gam.ovr')[:,0]
nz = min(len(blocks), len(dep))
r = 1e17 + dep[:nz]
edges = np.concatenate([[1e17], 0.5*(r[1:]+r[:-1]), [r[-1] + 0.5*(r[-1]-r[-2])]])
dV = 4/3*np.pi*(edges[1:]**3 - edges[:-1]**3)
tot = {}
for el, z in Z.items():
    f = np.loadtxt('gam.ele_'+el)[:nz, 1:z+2]
    per = []
    for ion in range(z):
        g = np.array([b.get((z, ion+1), 0.0) for b in blocks[:nz]])
        per.append(np.sum(g*ab[el]*nH*f[:,ion]*dV))
    tot[el] = per
print(f'zones {nz}')
for el, per in tot.items():
    s = sum(per)
    print(f'{el:3s} {s:11.4e} /Q = {s/Q:.4f}  by stage ' + ' '.join(f'{p/Q:.4f}' for p in per[:7]))
print(f"metals /Q = {sum(sum(tot[e]) for e in tot if e not in ('h','he'))/Q:.4f}")
