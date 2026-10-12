#!/usr/bin/env python3
"""Shell profile of one run: n_H V weighted fractions, T_e, [S III] 9532 share.
profile.py <run dir>"""
import sys, h5py, numpy as np
d = sys.argv[1]
r = h5py.File(d + '/bb1e5_rates.h5', 'r'); e = h5py.File(d + '/bb1e5_emis.h5', 'r')
xyz = r['LeafXYZ/data'][:]; rad = np.sqrt((xyz**2).sum(0)); V = r['LeafSize/data'][:]**3
nH = e['n_H/data'][:]; xs = e['x_s_stages/data'][:]; e9 = e['emis_s_3/data'][:, 5]
q = {k: r[k + '/data'][:] for k in ['x_HII', 'x_HeI', 'x_HeIII', 'T_e']}
ed = np.linspace(0, 5.3, 23); i = np.digitize(rad, ed); tot = np.sum(V*e9)
print('  r[pc]  x_HII   x_HeI   x_HeIII  S+    S2+   S3+    T_e   [S III] 9532 share')
for k in range(1, len(ed)):
    m = (i == k) & (nH > 0)
    if not m.any(): continue
    w = V[m]*nH[m]; f = lambda a: np.sum(w*a[m])/np.sum(w)
    print(f' {0.5*(ed[k-1]+ed[k]):5.2f}  {f(q["x_HII"]):.4f} {f(q["x_HeI"]):.2e} {f(q["x_HeIII"]):.3f}  '
          f'{f(xs[:,1]):.3f} {f(xs[:,2]):.3f} {f(xs[:,3]):.3f} {f(q["T_e"]):6.0f}  {np.sum(V[m]*e9[m])/tot:.3f}')
