#!/usr/bin/env python3
"""Markdown table of item effects against seed noise for one run set.
summarize.py <set>.  Columns are ratios B/A - 1 for the pairs below that
exist in the set: (A, B) = (a, d) channel 6, (b, d) channel 7, (c, d) Ne III
table, (d, noise) seed noise, (a_noise, noise) channel 6 at qmc_seed 777."""
import sys, os, h5py, numpy as np
W = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib/runs/' + sys.argv[1]
PAIRS = [('a', 'd', 'ch. 6 (d/a)'), ('b', 'd', 'ch. 7 (d/b)'), ('c', 'd', 'Ne III table (d/c)'),
         ('d', 'noise', 'seed noise (noise/d)'), ('a_noise', 'noise', 'ch. 6 at seed 777')]
PAIRS = [p for p in PAIRS if os.path.exists(f'{W}/{p[0]}/bb1e5_emis.h5') and os.path.exists(f'{W}/{p[1]}/bb1e5_emis.h5')]
LINES = [('s', '3', '9532.25', '[S III] 9532'), ('s', '3', '335008.38', '[S III] 33.5 um'),
         ('s', '3', '187055.74', '[S III] 18.7 um'), ('s', '2', '6718.29', '[S II] 6716'),
         ('o', '3', '5008.24', '[O III] 5007'), ('o', '2', '7322.20', '[O II] 7320'),
         ('n', '2', '6585.27', '[N II] 6584'), ('n', '1', '5201.71', '[N I] 5200'),
         ('ne', '3', '3869.85', '[Ne III] 3869'), ('ne', '4', '2422.39', '[Ne IV] 2422'),
         ('ne', '5', '3426.86', '[Ne V] 3426'), ('c', '2', '2325.41', 'C II] 2325'),
         ('he', '2', '4687.00', 'He II 4686'), ('he', '1', '5877.27', 'He I 5876')]
STAGES = [('c', 2, 'C2+'), ('n', 2, 'N2+'), ('o', 2, 'O2+'), ('ne', 3, 'Ne3+'), ('s', 1, 'S+'), ('s', 2, 'S2+')]
cache = {}
def lines(fn):
    L = {}; hb = None
    for l in open(fn):
        if l.startswith('# L(Hbeta) ='): hb = float(l.split('=')[1].split()[0])
        if l.startswith('#') or not l.strip(): continue
        p = l.split()
        try: L[(p[0], p[1], p[2])] = float(p[3])
        except Exception: pass
    return L, hb
def load(x):
    if x in cache: return cache[x]
    r = h5py.File(f'{W}/{x}/bb1e5_rates.h5', 'r'); e = h5py.File(f'{W}/{x}/bb1e5_emis.h5', 'r')
    d = {k: e[k + '/data'][:] for k in ['T_e', 'n_e', 'n_H', 'x_HeIII', 'LeafSize']}
    xyz = r['LeafXYZ/data'][:]; d['rad'] = np.sqrt((xyz**2).sum(0))
    for el in ['c', 'n', 'o', 'ne', 's']: d['x_' + el] = e[f'x_{el}_stages/data'][:]
    d['GS'] = r['Gamma_s_stages/data'][:, 1]
    d['L'], d['Hb'] = lines(f'{W}/{x}/bb1e5_lines.txt')
    cache[x] = d; return d
rows = {}
order = []
def put(name, vals):
    if name not in rows: order.append(name); rows[name] = []
    rows[name].append(vals)
for A, B, lab in PAIRS:
    a, b = load(A), load(B)
    V = a['LeafSize']**3; gas = a['n_H'] > 0
    put('L(Hbeta)', b['Hb']/a['Hb'] - 1)
    w = V*a['n_e']**2
    put('<T_e> (n_e^2 V)', np.sum(w*b['T_e'])/np.sum(w*a['T_e']) - 1)
    put('<x(He2+)> (n_e^2 V)', np.sum(w*b['x_HeIII'])/np.sum(w*a['x_HeIII']) - 1)
    for el, j, nm in STAGES:
        put(f'N({nm})', np.sum(V*a['n_H']*b['x_' + el][:, j])/np.sum(V*a['n_H']*a['x_' + el][:, j]) - 1)
    for lo, hi in ([(2.9, 4.0), (4.0, 5.3)] if sys.argv[1] == 'it2' else [(2.9, 3.6), (3.6, 4.1)]):
        m = gas & (a['rad'] >= lo) & (a['rad'] < hi)
        put(f'Gamma(S+), n_H V weighted, {lo}-{hi} pc', np.sum(V[m]*a['n_H'][m]*b['GS'][m])/np.sum(V[m]*a['n_H'][m]*a['GS'][m]) - 1)
    for k1, k2, k3, nm in LINES:
        k = (k1, k2, k3)
        put(nm + ' (L/Hbeta)', (b['L'][k]/b['Hb'])/(a['L'][k]/a['Hb']) - 1)
    big = [(b['L'][k]/b['Hb'])/(v/a['Hb']) - 1 for k, v in a['L'].items() if k in b['L'] and v/a['Hb'] > 1e-3]
    put('median abs, lines > 1e-3 Hbeta', np.median(np.abs(big)))
print('| quantity | ' + ' | '.join(p[2] for p in PAIRS) + ' |')
print('|---|' + '---|'*len(PAIRS))
for n in order:
    print(f'| {n} | ' + ' | '.join(f'{v:+.2e}' if not n.startswith('median') else f'{v:.1e}' for v in rows[n]) + ' |')
