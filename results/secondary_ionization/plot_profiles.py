#!/usr/bin/env python3
"""Radial profiles of the three secondary-ionization runs (volume-weighted
shells over the gas-bearing leaves): Te, x(H0), x(He0) and the differences
from the run without secondaries."""
import h5py, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
RUNS = [('off_new', 'off (current)', 'k', '-'),
        ('on_new', r'on, SvS85 $\times$ DYL99 (current)', 'C0', '-'),
        ('on_old', 'on, SvS85 from 40 eV (previous)', 'C3', '--')]
R = {}
for t, *_ in RUNS:
    h = h5py.File(f'sec_{t}_rates.h5', 'r')
    R[t] = {k: h[k + '/data'][:] for k in ('T_e', 'x_HI', 'x_HeI', 'LeafXYZ', 'LeafSize')}
o = R['off_new']; xyz = o['LeafXYZ']; xyz = xyz.T if xyz.shape[0] == 3 else xyz
r = np.sqrt((xyz**2).sum(-1)); V = o['LeafSize']**3
live = ~((o['x_HI'] == 1e-3) & (o['T_e'] == 1e4))      # leaves the solve never touched
edges = np.linspace(2.2, 2.9, 36); rc = 0.5*(edges[1:] + edges[:-1])
def prof(t, k):
    out = np.full(rc.size, np.nan)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        m = live & (r >= a) & (r < b)
        if m.sum(): out[i] = np.sum(R[t][k][m]*V[m])/V[m].sum()
    return out
fig, ax = plt.subplots(2, 2, figsize=(9, 6.5), sharex=True)
for t, lab, c, ls in RUNS:
    ax[0, 0].plot(rc, prof(t, 'T_e'), c=c, ls=ls, label=lab)
    ax[0, 1].semilogy(rc, prof(t, 'x_HI'), c=c, ls=ls)
    ax[0, 1].semilogy(rc, prof(t, 'x_HeI'), c=c, ls=ls, alpha=0.5)
    if t != 'off_new':
        ax[1, 0].plot(rc, prof(t, 'T_e') - prof('off_new', 'T_e'), c=c, ls=ls)
        ax[1, 1].plot(rc, prof(t, 'x_HI')/prof('off_new', 'x_HI') - 1, c=c, ls=ls)
ax[0, 0].set_ylabel(r'$T_{\rm e}$ [K]'); ax[0, 0].legend(fontsize=8)
ax[0, 1].set_ylabel(r'$x({\rm H^0})$ (dark), $x({\rm He^0})$ (light)')
ax[1, 0].set_ylabel(r'$\Delta T_{\rm e}$ from off [K]')
ax[1, 1].set_ylabel(r'$x({\rm H^0})/x({\rm H^0})_{\rm off} - 1$')
for a in ax[1]: a.set_xlabel('r [pc]'); a.axhline(0, c='0.6', lw=0.5)
fig.suptitle(r'HII20 shell grid, $10^5$ K blackbody, 80 iterations')
fig.tight_layout(); fig.savefig('sec_profiles.pdf'); fig.savefig('sec_profiles.png', dpi=80)
