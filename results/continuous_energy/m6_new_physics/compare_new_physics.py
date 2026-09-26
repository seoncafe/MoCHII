#!/usr/bin/env python3
"""Lexington HII40/HII20 with the 2026-09-25 physics (charge-exchange heating,
He I 2^3S coupling, collider-resolved two-photon quenching) against the same
binary with the three switches off (the control build, restarted from the
same m5 state) and against m5_converged itself; Cloudy c25.00 overlaid.

Run from this directory:  python3 compare_new_physics.py
"""
import re
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PC = 3.0856775814913673e18
R0 = 10**18.477121                      # Cloudy inner radius [cm]
CL = '../../../tests/cloudy_c25'
RUNS = {'m5': '../m5_converged/{t}_converged', 'control': '{t}_ctl',
        'new': '{t}_new'}


def t_npne(h5):
    with h5py.File(h5, 'r') as h:
        T = h['T_e/data'][:]; ne = h['n_e/data'][:]; xHI = h['x_HI/data'][:]
    g = (T > 0) & (ne > 0); w = (1.0 - xHI)*ne
    return float(np.sum(T[g]*w[g])/np.sum(w[g]))


def ionized_ratio(h5):
    """<He+>/<H+>, volume weighted with n_H constant, over the gas leaves."""
    with h5py.File(h5, 'r') as h:
        V = h['LeafSize/data'][:]**3; xHI = h['x_HI/data'][:]
        xHeII = h['x_HeII/data'][:]; ne = h['n_e/data'][:]
    g = ne > 0                          # empty leaves carry n_e = 0 (T_e > 0)
    return float(np.sum(V[g]*xHeII[g])/np.sum(V[g]*(1 - xHI[g])))


def lines(path):
    out, hb = {}, None
    for ln in open(path):
        m = re.search(r'L\(Hbeta\)\s*=\s*([0-9.eE+-]+)', ln)
        if m:
            hb = float(m.group(1)); continue
        if ln.startswith('#'):
            continue
        p = ln.split()
        if len(p) >= 5:
            out[(p[0], int(p[1]), float(p[2]))] = float(p[4])
    return hb, out


def profiles(h5, edges):
    with h5py.File(h5, 'r') as h:
        xyz = h['LeafXYZ/data'][:]; V = h['LeafSize/data'][:]**3
        r = np.sqrt((xyz**2).sum(axis=0))
        f = {'T': h['T_e/data'][:], 'H0': h['x_HI/data'][:],
             'He1': h['x_HeII/data'][:], 'He0': h['x_HeI/data'][:]}
        gas = h['n_e/data'][:] > 0          # empty leaves carry n_e = 0
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


def front(rc, x, level=0.5):
    """First radius where x falls through `level` (linear interpolation)."""
    for i in range(len(x) - 1):
        if np.isfinite(x[i]) and np.isfinite(x[i+1]) and x[i] >= level > x[i+1]:
            return rc[i] + (x[i] - level)/(x[i] - x[i+1])*(rc[i+1] - rc[i])
    return np.nan


def cloudy(tag):
    o = np.genfromtxt(f'{CL}/{tag}_c25.ovr', names=True, comments=None,
                      deletechars='')
    names = o.dtype.names
    r = (R0 + o[names[0]])/PC
    return dict(r=r, T=o['Te'], H0=o['HI'], He1=o['HeII'])


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
                      he=ionized_ratio(base + '_rates.h5'),
                      fH=front(rc, 1 - pr['H0']), fHe=front(rc, pr['He1']))
    for k in RUNS:
        r = res[k]
        report.append(f'  {k:8s} <T[NpNe]> {r["T"]:8.2f} K   L(Hb) {r["Hb"]:.5e}   '
                      f'<He+>/<H+> {r["he"]:.4f}   H front {r["fH"]:.4f} pc   '
                      f'He front {r["fHe"]:.4f} pc')
    n, c = res['new'], res['control']
    report.append(f'  new - control: dT {n["T"] - c["T"]:+.2f} K, L(Hb) '
                  f'{100*(n["Hb"]/c["Hb"] - 1):+.3f}%, He front '
                  f'{n["fHe"] - c["fHe"]:+.4f} pc, H front {n["fH"] - c["fH"]:+.4f} pc')
    report.append(f'  control - m5: dT {c["T"] - res["m5"]["T"]:+.2f} K, L(Hb) '
                  f'{100*(c["Hb"]/res["m5"]["Hb"] - 1):+.3f}% (restart band)')
    # lines relative to Hbeta that move by more than 0.5% new vs control
    moved = []
    for key, v in n['L'].items():
        v0 = c['L'].get(key)
        if v0 and v0 > 1e-4 and v > 0:
            d = v/v0 - 1
            if abs(d) > 5e-3:
                moved.append((abs(d), key, v0, v, d))
    moved.sort(reverse=True)
    report.append('  lines (rel. Hb > 1e-4) moving > 0.5%, new vs control:')
    for _, key, v0, v, d in moved[:15]:
        report.append(f'    {key[0]:>2s} {key[1]} {key[2]:12.2f}  {v0:.4e} -> {v:.4e}  ({100*d:+.2f}%)')
    cl = cloudy(tag)
    sty = {'m5': dict(c='0.6', ls='-', lw=1), 'control': dict(c='C0', ls='--', lw=1.2),
           'new': dict(c='C3', ls='-', lw=1.5)}
    for row, (key, lab, clk) in enumerate((('T', r'$T_e$ [K]', 'T'),
                                          ('H0', r'$x({\rm H}^0)$', 'H0'),
                                          ('He1', r'$x({\rm He}^+)$', 'He1'))):
        ax = axes[row, col]
        for k in RUNS:
            ax.plot(rc, res[k]['pr'][key], label=k, **sty[k])
        ax.plot(cl['r'], cl[clk], c='k', ls=':', lw=1.2, label='Cloudy c25.00')
        ax.set_ylabel(lab)
        if key == 'H0':
            ax.set_yscale('log'); ax.set_ylim(1e-5, 1.5)
        if row == 0:
            ax.set_title(tag.upper())
            ax.legend(fontsize=8)
    ax = axes[3, col]
    ax.plot(rc, res['new']['pr']['T'] - res['control']['pr']['T'], c='C3',
            lw=1.5, label='new $-$ control')
    ax.plot(rc, res['control']['pr']['T'] - res['m5']['pr']['T'], c='0.5',
            lw=1, label='control $-$ m5 (restart band)')
    ax.axhline(0, c='k', lw=0.5)
    ax.set_ylabel(r'$\Delta T_e$ [K]')
    ax.legend(fontsize=8)
    axes[3, col].set_xlabel(r'$r$ [pc]')
    axes[3, col].set_xlim(0.8 if tag == 'hii40' else 0.9, 5.2 if tag == 'hii40' else 3.0)
# HII40, each switch alone against the control build: the volume-weighted
# <T[NpNe]>, L(Hbeta), and the shell profile difference averaged over the
# He+ zone (1.2-4 pc) and at its minimum across the H front (4.4-4.9 pc).
report.append('=== HII40 single-switch attribution (each switch alone vs control) ===')
Tc = t_npne('hii40_ctl_rates.h5'); hbc, _ = lines('hii40_ctl_lines.txt')
prc = profiles('hii40_ctl_rates.h5', edges)['T']
in_he = (rc >= 1.2) & (rc <= 4.0); at_front = (rc >= 4.4) & (rc <= 4.9)
for k, base in (('only_cxh', 'hii40_only_cxh'), ('only_meta', 'hii40_only_meta'),
                ('only_q2p', 'hii40_only_q2p'), ('new', 'hii40_new')):
    hb, _ = lines(base + '_lines.txt')
    d = profiles(base + '_rates.h5', edges)['T'] - prc
    report.append(f'{k:10s} d<T[NpNe]> {t_npne(base + "_rates.h5") - Tc:+7.2f} K  '
                  f'L(Hb) {100*(hb/hbc - 1):+.3f}%  mean dT(1.2-4 pc) '
                  f'{np.nanmean(d[in_he]):+6.2f} K  min dT(4.4-4.9 pc) '
                  f'{np.nanmin(d[at_front]):+7.2f} K')
fig.tight_layout()
fig.savefig('new_physics_profiles.pdf')
print('\n'.join(report))
open('compare_new_physics.txt', 'w').write('\n'.join(report) + '\n')
