#!/usr/bin/env python3
"""T_e, x(H0) and x(C0) across the H front, m5 against m8 (transported
production inputs), Cloudy c25.00 overlaid; fine radial bins.

Run from this directory:  python3 front_zoom.py
"""
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

src = open('../m6_new_physics/compare_new_physics.py').read()
exec(src.split("edges = np.arange")[0])
CLD = '../../../tests/cloudy_c25'


def prof(h5, edges):
    with h5py.File(h5, 'r') as h:
        xyz = h['LeafXYZ/data'][:]; V = h['LeafSize/data'][:]**3
        r = np.sqrt((xyz**2).sum(axis=0)); gas = h['n_e/data'][:] > 0
        f = {'T': h['T_e/data'][:], 'H0': h['x_HI/data'][:],
             'C0': h['x_c_stages/data'][:][:, 0], 'O0': h['x_o_stages/data'][:][:, 0]}
    idx = np.where(gas, np.digitize(r, edges) - 1, -1)
    out = {}
    for k, v in f.items():
        a = np.full(edges.size - 1, np.nan)
        for i in range(a.size):
            s = idx == i
            if s.any():
                a[i] = np.average(v[s], weights=V[s])
        out[k] = a
    return out


fig, axes = plt.subplots(4, 2, figsize=(10, 11), sharex='col')
for col, (tag, lo, hi) in enumerate((('hii40', 4.2, 5.0), ('hii20', 2.5, 3.0))):
    edges = np.arange(lo, hi, 0.01)
    rc = 0.5*(edges[1:] + edges[:-1])
    cl = cloudy(tag)
    c = np.loadtxt(f'{CLD}/{tag}_c25.ele_c', comments='#')
    o = np.loadtxt(f'{CLD}/{tag}_c25.ele_o', comments='#')
    rcl_c = (R0 + c[:, 0])/PC; rcl_o = (R0 + o[:, 0])/PC
    for run, sty in (('../m5_converged', dict(c='0.55', lw=1, label='m5')),
                     ('.', dict(c='C3', lw=1.5, label='m8'))):
        p = prof(f'{run}/{tag}_converged_rates.h5', edges)
        axes[0, col].plot(rc, p['T'], **sty)
        axes[1, col].plot(rc, p['H0'], **sty)
        axes[2, col].plot(rc, p['C0'], **sty)
        axes[3, col].plot(rc, p['O0'], **sty)
    axes[0, col].plot(cl['r'], cl['T'], c='k', ls=':', label='Cloudy c25.00')
    axes[1, col].plot(cl['r'], cl['H0'], c='k', ls=':')
    axes[2, col].plot(rcl_c, c[:, 1], c='k', ls=':')
    axes[3, col].plot(rcl_o, o[:, 1], c='k', ls=':')
    axes[0, col].set_title(tag.upper()); axes[0, col].legend(fontsize=8)
    for row, lab in enumerate((r'$T_e$ [K]', r'$x({\rm H}^0)$', r'$x({\rm C}^0)$',
                               r'$x({\rm O}^0)$')):
        axes[row, col].set_ylabel(lab)
        if row:
            axes[row, col].set_yscale('log'); axes[row, col].set_ylim(1e-6, 1.5)
    axes[3, col].set_xlabel(r'$r$ [pc]'); axes[3, col].set_xlim(lo, hi)
fig.tight_layout()
fig.savefig('front_zoom.pdf')
