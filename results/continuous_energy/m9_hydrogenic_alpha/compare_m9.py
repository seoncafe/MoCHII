#!/usr/bin/env python3
"""The m8 paper set restarted with the hydrogenic recombination of
2026-09-26 (m9): alpha_B of H II and He III from the zero-density Storey &
Hummer cascade (1.13% below the Badnell-minus-Milne value for H II at
1e4 K), the joint H0/He0 solve of the He II on-the-spot cascade, and the
written-continuum changes (2s share by case, table gaps), against the m8
states it restarted from; Cloudy c25.00 overlaid.  Scalar table for all
eight HII40/HII20 cases, line ratios for the production pair, profiles for
the production and on-the-spot pairs.

Run from this directory:  python3 compare_m9.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

src = open('../m6_new_physics/compare_new_physics.py').read()
exec(src.split("edges = np.arange")[0])     # t_npne, ionized_ratio, lines, profiles, front, cloudy

CASES = ('converged', 'ots_converged', 'lowden_converged', 'nohecx_converged')
edges = np.arange(0.0, 6.0, 0.02)
rc = 0.5*(edges[1:] + edges[:-1])
report = []


def state(base):
    hb, L = lines(base + '_lines.txt')
    pr = profiles(base + '_rates.h5', edges)
    return dict(T=t_npne(base + '_rates.h5'), Hb=hb, L=L, pr=pr,
                he=ionized_ratio(base + '_rates.h5'),
                fH=front(rc, 1 - pr['H0']), fHe=front(rc, pr['He1']))


res = {}
for tag in ('hii40', 'hii20'):
    report.append(f'=== {tag.upper()} ===')
    for case in CASES:
        a = state(f'../m8_defaults/{tag}_{case}')
        b = state(f'{tag}_{case}')
        res[(tag, case)] = (a, b)
        report.append(f'  {case:17s} <T[NpNe]> {a["T"]:8.2f} -> {b["T"]:8.2f} K ({b["T"] - a["T"]:+7.2f})'
                      f'  L(Hb) {a["Hb"]:.4e} -> {b["Hb"]:.4e} ({100*(b["Hb"]/a["Hb"] - 1):+6.2f}%)'
                      f'  <He+>/<H+> {a["he"]:.4f} -> {b["he"]:.4f}'
                      f'  H front {a["fH"]:.4f} -> {b["fH"]:.4f}  He front {a["fHe"]:.4f} -> {b["fHe"]:.4f} pc')
    a, b = res[(tag, 'converged')]
    moved = []
    for key, v in b['L'].items():
        v0 = a['L'].get(key)
        if v0 and v0 > 1e-3 and v > 0 and abs(v/v0 - 1) > 5e-3:
            moved.append((abs(v/v0 - 1), key, v0, v, v/v0 - 1))
    moved.sort(reverse=True)
    report.append('  production lines (rel. Hb > 1e-3) moving > 0.5%, m8 -> m9:')
    for _, key, v0, v, d in moved[:16]:
        report.append(f'    {key[0]:>2s} {key[1]} {key[2]:12.2f}  {v0:.4e} -> {v:.4e}  ({100*d:+.2f}%)')

fig, axes = plt.subplots(4, 2, figsize=(10, 12.5), sharex='col')
sty = {('converged', 0): dict(c='0.55', ls='-', lw=1, label='m8, transported'),
       ('converged', 1): dict(c='C3', ls='-', lw=1.5, label='m9, transported'),
       ('ots_converged', 0): dict(c='C0', ls=':', lw=1, label='m8, on the spot'),
       ('ots_converged', 1): dict(c='C0', ls='--', lw=1.3, label='m9, on the spot')}
for col, tag in enumerate(('hii40', 'hii20')):
    cl = cloudy(tag)
    for row, (key, ylab) in enumerate((('T', r'$T_e$ [K]'), ('H0', r'$x({\rm H}^0)$'),
                                       ('He1', r'$x({\rm He}^+)$'))):
        ax = axes[row, col]
        for (case, i), s in sty.items():
            ax.plot(rc, res[(tag, case)][i]['pr'][key], **s)
        ax.plot(cl['r'], cl[key], c='k', ls=':', lw=1.2, label='Cloudy c25.00')
        ax.set_ylabel(ylab)
        if key == 'H0':
            ax.set_yscale('log'); ax.set_ylim(1e-5, 1.5)
        if row == 0:
            ax.set_title(tag.upper()); ax.legend(fontsize=7)
    ax = axes[3, col]
    for case, c in (('converged', 'C3'), ('ots_converged', 'C0')):
        a, b = res[(tag, case)]
        ax.plot(rc, b['pr']['T'] - a['pr']['T'], c=c, lw=1.3,
                label='m9 $-$ m8, ' + ('transported' if case == 'converged' else 'on the spot'))
    ax.axhline(0, c='k', lw=0.5)
    ax.set_ylabel(r'$\Delta T_e$ [K]'); ax.legend(fontsize=7)
    ax.set_xlabel(r'$r$ [pc]')
    ax.set_xlim(0.8 if tag == 'hii40' else 0.9, 5.2 if tag == 'hii40' else 3.0)
fig.tight_layout()
fig.savefig('m9_profiles.pdf')
print('\n'.join(report))
open('compare_m9.txt', 'w').write('\n'.join(report) + '\n')
