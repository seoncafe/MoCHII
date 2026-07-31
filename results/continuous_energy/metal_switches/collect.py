#!/usr/bin/env python3
"""Fill the G-M2 switch matrix.

The front positions come from tests/continuous_energy/analyze_hii40_residuals.py
itself -- volume-weighted 0.01-pc shells and its own crossing rule -- so these
numbers are comparable with the recorded residual measurements instead of being
a second implementation of the same diagnostic.
"""
import os, sys, re
import numpy as np, h5py
ROOT = '/nfs/mocafe/kiseon/RT_Codes/MoCHII'
sys.path.insert(0, os.path.join(ROOT, 'tests/continuous_energy'))
import analyze_hii40_residuals as ana

def t_npne(h5):
    with h5py.File(h5, 'r') as h:
        T = h['T_e/data'][:]; ne = h['n_e/data'][:]; xHI = h['x_HI/data'][:]
    g = (T > 0) & (ne > 0)
    w = 100.0 * (1.0 - xHI) * ne
    return float(np.sum(T[g] * w[g]) / np.sum(w[g]))

def l_hbeta(path):
    if not os.path.exists(path): return float('nan')
    for ln in open(path):
        m = re.search(r'L\(Hbeta\)\s*=\s*([0-9.eE+-]+)', ln)
        if m: return float(m.group(1))
    return float('nan')

CASES = [('on', 'on', 'production/{b}_continuous'),
         ('on', 'off', 'metal_switches/{b}_heat_off'),
         ('off', 'on', 'metal_switches/{b}_ne_off'),
         ('off', 'off', 'metal_switches/{b}_both_off')]

for tag, b in (('HII40', 'hii40'), ('HII20', 'hii20')):
    print(f'\n=== {tag}: G-M2 switch matrix ===')
    print(f'{"ne":4s} {"heat":5s} {"<T[NpNe]> K":>12s} {"L(Hbeta)":>12s} '
          f'{"r(He0=0.5)":>11s} {"r(C2=0.5)":>10s} {"r(H0=0.1)":>10s}')
    base = None
    for ne_s, ht_s, stem in CASES:
        stem = stem.format(b=b)
        h5 = os.path.join(ROOT, 'results/continuous_energy', stem + '_rates.h5')
        if not os.path.exists(h5):
            print(f'{ne_s:4s} {ht_s:5s} {"(pending)":>12s}'); continue
        prof = ana.read_mochii(h5)
        row = [t_npne(h5), l_hbeta(os.path.join(ROOT,'results/continuous_energy',stem+'_lines.txt')),
               ana.outward_crossing(prof['He0'], 0.5),
               ana.outward_crossing(prof['C2'], 0.5),
               ana.outward_crossing(prof['H0'], 0.1)]
        if base is None: base = row
        d = f'   [dT = {row[0]-base[0]:+7.1f} K, dL/L = {100*(row[1]/base[1]-1):+6.3f}%, dr(He0) = {row[2]-base[2]:+7.4f} pc]'
        print(f'{ne_s:4s} {ht_s:5s} {row[0]:12.1f} {row[1]:12.4e} {row[2]:11.4f} {row[3]:10.4f} {row[4]:10.4f}{"" if base is row else d}')
