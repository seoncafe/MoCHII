#!/usr/bin/env python3
"""Case-B on-the-spot He I recombination photons (par%hei_on_the_spot):
the Lexington HII40/HII20 on-the-spot inputs (m5 *_ots_converged) restarted
with the treatment on and off, against the transported-continuum runs of
m6_new_physics (*_new, the physical reference the on-the-spot approximation
stands in for) and Cloudy c25.00.

Run from this directory:  python3 compare_hei_on_the_spot.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

src = open('../m6_new_physics/compare_new_physics.py').read()
exec(src.split("edges = np.arange")[0])     # t_npne, lines, profiles, front, cloudy

RUNS = {'off': '{t}_ots_off', 'on': '{t}_ots_on',
        'transported': '../m6_new_physics/{t}_new'}
edges = np.arange(0.0, 6.0, 0.02)
rc = 0.5*(edges[1:] + edges[:-1])
fig, axes = plt.subplots(4, 2, figsize=(10, 12.5), sharex='col')
report = []
for col, tag in enumerate(('hii40', 'hii20')):
    report.append(f'=== {tag.upper()} ===')
    res = {}
    for k, pat in RUNS.items():
        base = pat.format(t=tag)
        hb, L = lines(base + '_lines.txt')
        pr = profiles(base + '_rates.h5', edges)
        res[k] = dict(T=t_npne(base + '_rates.h5'), Hb=hb, L=L, pr=pr,
                      fH=front(rc, 1 - pr['H0']), fHe=front(rc, pr['He1']))
    for k in RUNS:
        r = res[k]
        report.append(f'  {k:12s} <T[NpNe]> {r["T"]:8.2f} K   L(Hb) {r["Hb"]:.5e}   '
                      f'H front {r["fH"]:.4f} pc   He front {r["fHe"]:.4f} pc')
    on, off, tr = res['on'], res['off'], res['transported']
    report.append(f'  on - off: dT {on["T"] - off["T"]:+.2f} K, L(Hb) '
                  f'{100*(on["Hb"]/off["Hb"] - 1):+.3f}%, H front '
                  f'{on["fH"] - off["fH"]:+.4f} pc, He front {on["fHe"] - off["fHe"]:+.4f} pc')
    report.append(f'  distance to transported: off dT {off["T"] - tr["T"]:+.2f} K, '
                  f'He front {off["fHe"] - tr["fHe"]:+.4f} pc, H front {off["fH"] - tr["fH"]:+.4f} pc; '
                  f'on dT {on["T"] - tr["T"]:+.2f} K, He front {on["fHe"] - tr["fHe"]:+.4f} pc, '
                  f'H front {on["fH"] - tr["fH"]:+.4f} pc')
    moved = []
    for key, v in on['L'].items():
        v0 = off['L'].get(key)
        if v0 and v0 > 1e-3 and v > 0 and abs(v/v0 - 1) > 5e-3:
            moved.append((abs(v/v0 - 1), key, v0, v, v/v0 - 1, tr['L'].get(key, np.nan)))
    moved.sort(reverse=True)
    report.append('  lines (rel. Hb > 1e-3) moving > 0.5%, on vs off (transported shown):')
    for _, key, v0, v, d, vt in moved[:14]:
        report.append(f'    {key[0]:>2s} {key[1]} {key[2]:12.2f}  {v0:.4e} -> {v:.4e} '
                      f'({100*d:+.2f}%)   transported {vt:.4e}')
    cl = cloudy(tag)
    sty = {'off': dict(c='C0', ls='--', lw=1.2), 'on': dict(c='C3', ls='-', lw=1.5),
           'transported': dict(c='0.4', ls='-', lw=1)}
    lab = {'off': 'on the spot, He I photons to He', 'on': 'on the spot, He I photons to H',
           'transported': 'transported (m6)'}
    for row, (key, ylab, clk) in enumerate((('T', r'$T_e$ [K]', 'T'),
                                           ('H0', r'$x({\rm H}^0)$', 'H0'),
                                           ('He1', r'$x({\rm He}^+)$', 'He1'))):
        ax = axes[row, col]
        for k in RUNS:
            ax.plot(rc, res[k]['pr'][key], label=lab[k], **sty[k])
        ax.plot(cl['r'], cl[clk], c='k', ls=':', lw=1.2, label='Cloudy c25.00')
        ax.set_ylabel(ylab)
        if key == 'H0':
            ax.set_yscale('log'); ax.set_ylim(1e-5, 1.5)
        if row == 0:
            ax.set_title(tag.upper()); ax.legend(fontsize=7)
    ax = axes[3, col]
    ax.plot(rc, on['pr']['T'] - off['pr']['T'], c='C3', lw=1.5, label='on $-$ off')
    ax.plot(rc, off['pr']['T'] - tr['pr']['T'], c='C0', ls='--', lw=1, label='off $-$ transported')
    ax.plot(rc, on['pr']['T'] - tr['pr']['T'], c='0.3', lw=1, label='on $-$ transported')
    ax.axhline(0, c='k', lw=0.5)
    ax.set_ylabel(r'$\Delta T_e$ [K]'); ax.legend(fontsize=7)
    ax.set_xlabel(r'$r$ [pc]')
    ax.set_xlim(0.8 if tag == 'hii40' else 0.9, 5.2 if tag == 'hii40' else 3.0)
fig.tight_layout()
fig.savefig('hei_on_the_spot_profiles.pdf')
print('\n'.join(report))
open('compare_hei_on_the_spot.txt', 'w').write('\n'.join(report) + '\n')
