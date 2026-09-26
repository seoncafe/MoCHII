#!/usr/bin/env python3
"""Spherically averaged T_e, x(H0), x(He+), x(He2+) of the three runs
(ots_on, ots_off, transported) -> heii_ots_profiles.pdf/.png.
Run from this directory after the three runs."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze import load

RUNS = [('transported', 'transported (reference)', 'k', '-'),
        ('ots_on', r'on the spot, \texttt{heii\_on\_the\_spot = T}', 'C0', '--'),
        ('ots_off', r'on the spot, \texttt{heii\_on\_the\_spot = F}', 'C3', ':')]
FIELDS = [('T_e', r'$T_e$ [K]', False),
          ('x_HI', r'$x({\rm H^0})$', True),
          ('x_HeII', r'$x({\rm He^+})$', False),
          ('xHeIII', r'$x({\rm He^{2+}})$', False)]


def shells(d, key, edges):
    r, v, x = d['r'], d['vol'], d[key]
    i = np.digitize(r, edges) - 1
    ok = (i >= 0) & (i < len(edges) - 1) & (d['n_e'] > 0)
    num = np.bincount(i[ok], weights=(v*x)[ok], minlength=len(edges) - 1)
    den = np.bincount(i[ok], weights=v[ok], minlength=len(edges) - 1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return num/den


fig, axs = plt.subplots(2, 2, figsize=(10, 7.5), sharex=True)
for base, label, col, ls in RUNS:
    d = load(base)
    edges = np.linspace(0.0, 5.2, 105)
    rc = 0.5*(edges[1:] + edges[:-1])
    for ax, (key, ylab, logy) in zip(axs.flat, FIELDS):
        ax.plot(rc, shells(d, key, edges), color=col, ls=ls, lw=1.4, label=label)
        ax.set_ylabel(ylab)
        if logy:
            ax.set_yscale('log')
for ax in axs[1]:
    ax.set_xlabel(r'$r$ [pc]')
axs[0, 0].legend(fontsize=8, loc='upper right')
fig.suptitle(r'$10^5$ K blackbody, HII40 geometry: He\,\textsc{ii} recombination photons')
fig.tight_layout()
fig.savefig('heii_ots_profiles.pdf')
fig.savefig('heii_ots_profiles.png', dpi=110)
