#!/usr/bin/env python3
"""Phase 4 attribution: compare the run of build d (all items on) with the
control builds and the seed-noise run.  Usage: compare.py <set> (it2 | rs)."""
import sys, h5py, numpy as np
W = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib/runs'
SET = sys.argv[1]
OTHERS = [x for x in ['a', 'b', 'c', 'noise'] if len(sys.argv) < 3 or x in sys.argv[2:]]
ABUND = {'c': 2.2e-4, 'n': 4.0e-5, 'o': 3.3e-4, 'ne': 5.0e-5, 's': 9.0e-6}

def load(x):
    r = h5py.File(f'{W}/{SET}/{x}/bb1e5_rates.h5', 'r')
    e = h5py.File(f'{W}/{SET}/{x}/bb1e5_emis.h5', 'r')
    d = {}
    for k in ['T_e', 'n_e', 'n_H', 'x_HI', 'x_HII', 'x_HeI', 'x_HeII', 'x_HeIII', 'LeafSize']:
        d[k] = e[k + '/data'][:]
    for el in ABUND:
        d['x_' + el] = e[f'x_{el}_stages/data'][:]
        d['G_' + el] = r[f'Gamma_{el}_stages/data'][:]
    d['s3'] = e['emis_s_3/data'][:]
    d['s2'] = e['emis_s_2/data'][:]
    d['wl_s3'] = e['wl_s_3/data'][:]
    d['J'] = r['J_nu/data'][:]
    d['E'] = r['E_bin/data'][:]
    d['dE'] = r['dE_bin/data'][:]
    d['lines'], d['Hb'] = lines(f'{W}/{SET}/{x}/bb1e5_lines.txt')
    d['log'] = open(f'{W}/{SET}/{x}/run.log').read()
    return d

def lines(fn):
    L = {}; hb = None
    for l in open(fn):
        if l.startswith('# L(Hbeta) ='):
            hb = float(l.split('=')[1].split()[0])
        if l.startswith('#') or not l.strip():
            continue
        p = l.split()
        key = (p[0], p[1], p[2]) if len(p) >= 4 else None
        try:
            L[key] = float(p[3].replace('E+0', 'E+').replace('E-0', 'E-'))
        except Exception:
            pass
    return L, hb

D = load('d')
V = D['LeafSize']**3
gas = D['n_H'] > 0

def chan(log):
    rows = [l for l in log.splitlines() if 'diffuse L by channel' in l]
    return np.array([float(v) for v in rows[-1].split(':')[1].split()])

print(f'== set {SET}: effect of an item = d/X - 1 (X = build without it); noise = noise/d - 1')
print('d final-pass channel 1-7 L [erg/s]:', ' '.join(f'{v:.4e}' for v in chan(D['log'])))
for X in OTHERS:
    A = load(X)
    sgn = (X == 'noise')
    B1, B0 = (A, D) if sgn else (D, A)       # ratio B1/B0 - 1
    print(f'\n######## {X}  ({"noise/d" if sgn else "d/" + X} - 1)')
    print(f'{X} final-pass channel L:', ' '.join(f'{v:.4e}' for v in chan(A['log'])))
    print(f'L(Hbeta) {B1["Hb"]/B0["Hb"]-1:+.3e}')
    w = (V*D['n_e']**2)[gas]
    for k in ['T_e', 'n_e', 'x_HII', 'x_HeII', 'x_HeIII']:
        print(f'  <{k}>_ne2V {np.sum(w*B1[k][gas])/np.sum(w*B0[k][gas])-1:+.3e}')
    # stage populations: sum V n_H x
    print('  stage populations sum V n_H x (share in d > 1e-3):')
    for el in ABUND:
        s = []
        tot = np.sum(V*D['n_H'])
        for j in range(D['x_' + el].shape[1]):
            p0 = np.sum(V*D['n_H']*B0['x_' + el][:, j]); p1 = np.sum(V*D['n_H']*B1['x_' + el][:, j])
            if np.sum(V*D['n_H']*D['x_' + el][:, j])/tot > 1e-3:
                s.append(f'{el}{j}+ {p1/p0-1:+.2e}')
        print('   ', '  '.join(s))
    # lines
    big = []
    for k, v in B0['lines'].items():
        if k in B1['lines'] and v/B0['Hb'] > 1e-3:
            big.append((k, (B1['lines'][k]/B1['Hb'])/(v/B0['Hb']) - 1, v/B0['Hb']))
    big.sort(key=lambda t: -abs(t[1]))
    print(f'  lines > 1e-3 Hbeta: {len(big)}, median |rel L/Hb| {np.median([abs(t[2-1]) for t in big]):.2e}')
    for k, dr, r in big[:12]:
        print(f'    {k[0]:>3s} {k[1]:>2s} {float(k[2]):12.2f}  {dr:+.3e}  L/Hb {r:.3e}')
    for k in [('s', '3', '9532.25'), ('s', '3', '335008.38'), ('s', '3', '187055.74'), ('s', '2', '6718.29'), ('s', '2', '6732.67')]:
        print(f'    S line {k}: {(B1["lines"][k]/B1["Hb"])/(B0["lines"][k]/B0["Hb"])-1:+.3e} (L/Hb)  {B1["lines"][k]/B0["lines"][k]-1:+.3e} (L)')

    # [S III] 9532 leaf decomposition (ln-epsilon split)
    i9 = int(np.argmin(abs(D['wl_s3'] - 9532.25)))
    e0 = B0['s3'][:, i9]; e1 = B1['s3'][:, i9]
    m = (e0 > 0) & (e1 > 0)
    dL = np.sum(V[m]*(e1[m] - e0[m])); L0 = np.sum(V*e0)
    em = 0.5*(e0 + e1)
    dlx = np.log(B1['x_s'][m, 2]/B0['x_s'][m, 2])
    dln = np.log(B1['n_e'][m]/B0['n_e'][m])
    dle = np.log(e1[m]/e0[m])
    dlf = dle - dlx - dln
    print(f'  [S III] 9532 sum V eps: dL/L {dL/L0:+.3e}; split (eps-weighted dln): x(S2+) {np.sum(V[m]*em[m]*dlx)/L0:+.3e}'
          f'  n_e {np.sum(V[m]*em[m]*dln)/L0:+.3e}  rest (T_e, n_e in the level solve) {np.sum(V[m]*em[m]*dlf)/L0:+.3e}')
    # regions by x_HII (d)
    xh = D['x_HII']; xh3 = D['x_HeIII']
    regs = [('He2+ zone x_HeIII>0.5', xh3 > 0.5), ('He+ zone x_HeIII<0.5, x_HII>0.9', (xh3 <= 0.5) & (xh > 0.9)),
            ('H front 0.1<x_HII<0.9', (xh > 0.1) & (xh <= 0.9)), ('neutral x_HII<0.1', xh <= 0.1)]
    for nm, r in regs:
        rr = r & m
        print(f'    {nm:34s} share of L {np.sum(V[r]*e0[r])/L0:.3f}  dL/L {np.sum(V[rr]*(e1[rr]-e0[rr]))/L0:+.3e}'
              f'  dln x(S2+) {np.sum(V[rr]*em[rr]*np.log(B1["x_s"][rr,2]/B0["x_s"][rr,2]))/L0:+.3e}')
    # S+ -> S2+ photoionization and field in the S+ window, weighted by the 9532 emissivity
    ww = (V*em)[m]
    for j, nm in [(1, 'Gamma(S+)'), (2, 'Gamma(S2+)')]:
        g0 = B0['G_s'][m, j]; g1 = B1['G_s'][m, j]
        print(f'    eps-weighted mean dln {nm}: {np.sum(ww*np.log(g1/g0))/np.sum(ww):+.3e}')
    for k in ['T_e', 'n_e', 'x_HI']:
        print(f'    eps-weighted mean dln {k}: {np.sum(ww*np.log(B1[k][m]/B0[k][m]))/np.sum(ww):+.3e}')
    # ratio x(S2+)/x(S+) vs Gamma(S+)
    r0 = B0['x_s'][m, 2]/B0['x_s'][m, 1]; r1 = B1['x_s'][m, 2]/B1['x_s'][m, 1]
    print(f'    eps-weighted mean dln[x(S2+)/x(S+)]: {np.sum(ww*np.log(r1/r0))/np.sum(ww):+.3e}')
    # J_nu in bins, weighted by n(S+) V and by eps9532 V
    wsp = V*D['n_H']*D['x_s'][:, 1]
    E = D['E']
    sel = (E > 20) & (E < 60)
    print('    J_nu bin ratio (B1/B0 - 1), weights: n(S+)V | eps9532 V | all-gas V')
    wE = np.zeros_like(V); wE[m] = V[m]*em[m]
    for ib in np.where(sel)[0]:
        j0 = B0['J'][:, ib]; j1 = B1['J'][:, ib]
        a = np.sum(wsp*j1)/np.sum(wsp*j0) - 1
        b = np.sum(wE*j1)/np.sum(wE*j0) - 1
        c = np.sum(V*j1)/np.sum(V*j0) - 1
        print(f'      E {E[ib]:7.2f} dE {D["dE"][ib]:6.2f}  {a:+.3e}  {b:+.3e}  {c:+.3e}')
