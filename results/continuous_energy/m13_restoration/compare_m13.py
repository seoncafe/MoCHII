#!/usr/bin/env python3
"""Lexington HII40 / HII20 under the atomic physics of 2026-10-11 (m13: stages
I-VII of every metal, full CHIANTI 11 atoms, phfit2 subshell absorption,
diffuse channels 6 and 7, Ne III from Nahar 2019) against the m11 states
they restarted from; Cloudy c25.00 overlaid, including the metal stages m11
did not carry (O3+, Ne3+, ...).

Run from this directory:  python3 compare_m13.py
Writes compare_m13.txt and m13_profiles.pdf.
"""
import re
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

src = open('../m6_new_physics/compare_new_physics.py').read()
exec(src.split("edges = np.arange")[0])     # t_npne, ionized_ratio, lines, profiles, front, cloudy

M11 = '../m11_hecx_solver/{t}_converged'
M13 = '{t}_converged'
ELEMENTS = ('c', 'n', 'o', 'ne', 's')
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII']
edges = np.arange(0.0, 6.0, 0.02)
rc = 0.5*(edges[1:] + edges[:-1])
report = []


def iterations(base):
    txt = open(base + '.log').read()
    return len(re.findall(r'^     iteration +\d+:', txt, re.M))


def stage_profiles(h5, el):
    """Volume-weighted spherical average of every stage of element el."""
    with h5py.File(h5, 'r') as h:
        xyz = h['LeafXYZ/data'][:]; V = h['LeafSize/data'][:]**3
        x = h[f'x_{el}_stages/data'][:]
        gas = h['n_e/data'][:] > 0
    r = np.sqrt((xyz**2).sum(axis=0))
    idx = np.where(gas, np.digitize(r, edges) - 1, -1)
    out = np.full((edges.size - 1, x.shape[1]), np.nan)
    for i in range(edges.size - 1):
        s = idx == i
        if s.any():
            out[i] = np.average(x[s], weights=V[s], axis=0)
    return out


def volume_mean_stages(h5, el):
    """Volume-weighted mean stage fractions over the ionized gas (x_HII > 0.5)."""
    with h5py.File(h5, 'r') as h:
        V = h['LeafSize/data'][:]**3; x = h[f'x_{el}_stages/data'][:]
        g = (h['n_e/data'][:] > 0) & (h['x_HII/data'][:] > 0.5)
    return np.average(x[g], weights=V[g], axis=0)


def cloudy_stages(tag, el):
    path = f'{CL}/{tag}_c25.ele_{el}'
    try:
        d = np.loadtxt(path)
    except OSError:
        return None
    return (R0 + d[:, 0])/PC, d[:, 1:8]


def state(base):
    hb, L = lines(base + '_lines.txt')
    pr = profiles(base + '_rates.h5', edges)
    return dict(T=t_npne(base + '_rates.h5'), Hb=hb, L=L, pr=pr,
                he=ionized_ratio(base + '_rates.h5'),
                fH=front(rc, 1 - pr['H0']), fHe=front(rc, pr['He1']),
                st={el: stage_profiles(base + '_rates.h5', el) for el in ELEMENTS},
                vm={el: volume_mean_stages(base + '_rates.h5', el) for el in ELEMENTS})


res = {}
for tag in ('hii40', 'hii20'):
    a = state(M11.format(t=tag)); b = state(M13.format(t=tag))
    res[tag] = (a, b)
    report.append(f'=== {tag.upper()} ===')
    report.append(f'  iterations after restart: {iterations(M13.format(t=tag))}')
    report.append(f'  <T[NpNe]> {a["T"]:8.2f} -> {b["T"]:8.2f} K ({b["T"] - a["T"]:+7.2f})'
                  f'  L(Hb) {a["Hb"]:.4e} -> {b["Hb"]:.4e} ({100*(b["Hb"]/a["Hb"] - 1):+6.2f}%)')
    report.append(f'  <He+>/<H+> {a["he"]:.4f} -> {b["he"]:.4f}'
                  f'  H front {a["fH"]:.4f} -> {b["fH"]:.4f}  He front {a["fHe"]:.4f} -> {b["fHe"]:.4f} pc')
    report.append('  volume-mean stage fractions over x_HII > 0.5 (m11 -> m13):')
    for el in ELEMENTS:
        s11 = ' '.join(f'{v:.4f}' for v in a['vm'][el])
        s13 = ' '.join(f'{v:.4f}' for v in b['vm'][el])
        report.append(f'    {el:2s} {s11}  ->  {s13}')
    moved = []
    for key, v in b['L'].items():
        v0 = a['L'].get(key)
        if v0 and v0 > 1e-3 and v > 0 and abs(v/v0 - 1) > 5e-3:
            moved.append((abs(v/v0 - 1), key, v0, v, v/v0 - 1))
    new = sorted(((v, key) for key, v in b['L'].items() if key not in a['L'] and v > 1e-3),
                 reverse=True)
    moved.sort(reverse=True)
    report.append(f'  lines (rel. Hb > 1e-3) moving > 0.5%, m11 -> m13: {len(moved)}')
    for _, key, v0, v, d in moved[:30]:
        report.append(f'    {key[0]:>2s} {key[1]} {key[2]:12.2f}  {v0:.4e} -> {v:.4e}  ({100*d:+.2f}%)')
    report.append(f'  lines in m13 only (rel. Hb > 1e-3): {len(new)}')
    for v, key in new[:20]:
        report.append(f'    {key[0]:>2s} {key[1]} {key[2]:12.2f}  {v:.4e}')

fig, axes = plt.subplots(3 + len(ELEMENTS), 2, figsize=(10, 3.0*(3 + len(ELEMENTS))), sharex='col')
for col, tag in enumerate(('hii40', 'hii20')):
    a, b = res[tag]
    cl = cloudy(tag)
    for row, (key, ylab) in enumerate((('T', r'$T_e$ [K]'), ('H0', r'$x({\rm H}^0)$'),
                                       ('He1', r'$x({\rm He}^+)$'))):
        ax = axes[row, col]
        ax.plot(rc, a['pr'][key], c='0.55', lw=1, label='m11')
        ax.plot(rc, b['pr'][key], c='C3', lw=1.5, label='m13')
        ax.plot(cl['r'], cl[key], c='k', ls=':', lw=1.2, label='Cloudy c25.00')
        ax.set_ylabel(ylab)
        if key == 'H0':
            ax.set_yscale('log'); ax.set_ylim(1e-5, 1.5)
        if row == 0:
            ax.set_title(tag.upper()); ax.legend(fontsize=7)
    for k, el in enumerate(ELEMENTS):
        ax = axes[3 + k, col]
        c = cloudy_stages(tag, el)
        for j in range(b['st'][el].shape[1]):
            ax.plot(rc, b['st'][el][:, j], c=f'C{j}', lw=1.4,
                    label=r'${\rm ' + el.capitalize() + r'}\,$' + ROMAN[j])
            if j < a['st'][el].shape[1]:
                ax.plot(rc, a['st'][el][:, j], c=f'C{j}', lw=0.8, ls='--')
            if c is not None:
                ax.plot(c[0], c[1][:, j], c=f'C{j}', lw=1.0, ls=':')
        ax.set_yscale('log'); ax.set_ylim(1e-4, 1.3)
        ax.set_ylabel(r'$x({\rm ' + el.capitalize() + r'})$')
        ax.legend(fontsize=6, ncol=4, loc='lower left')
    axes[-1, col].set_xlabel(r'$r$ [pc] (solid m13, dashed m11, dotted Cloudy)')
    axes[-1, col].set_xlim(0.8 if tag == 'hii40' else 0.9, 5.2 if tag == 'hii40' else 3.0)
fig.tight_layout()
fig.savefig('m13_profiles.pdf')
print('\n'.join(report))
open('compare_m13.txt', 'w').write('\n'.join(report) + '\n')
