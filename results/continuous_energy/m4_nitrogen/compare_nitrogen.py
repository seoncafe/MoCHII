#!/usr/bin/env python3
"""Nitrogen charge exchange (Huang et al. 2023 rows -> Kingdon & Ferland 1996):
ionization profiles and lines of the Lexington runs before (m2_atomic) and
after (m4_nitrogen), with Cloudy c25.00 (tests/cloudy_c25) overlaid."""
import re, sys, numpy as np, h5py, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '../metal_switches')
PC = 3.0856775814913673e18
R0 = 10**18.477121
CL = '../../../tests/cloudy_c25'
def t_npne(h5):
    with h5py.File(h5, 'r') as h:
        T = h['T_e/data'][:]; ne = h['n_e/data'][:]; xHI = h['x_HI/data'][:]
    g = (T > 0) & (ne > 0); w = 100.0*(1.0 - xHI)*ne
    return float(np.sum(T[g]*w[g])/np.sum(w[g]))
def lines(path):
    out = {}
    for ln in open(path):
        m = re.search(r'L\(Hbeta\)\s*=\s*([0-9.eE+-]+)', ln)
        if m: out['Hb'] = float(m.group(1)); continue
        p = ln.split()
        if len(p) >= 5 and p[0] in ('n', 'o', 'c', 's') and not ln.startswith('#'):
            out[(p[0], int(p[1]), float(p[2]))] = float(p[4])
    return out
def pick(L, el, st, lam, tol=3.0):
    v = [val for k, val in L.items() if k != 'Hb' and k[0] == el and k[1] == st
         and abs(k[2] - lam) < tol*max(1.0, lam/5000.)]
    return sum(v) if v else np.nan
def profiles(h5, edges):
    with h5py.File(h5, 'r') as h:
        xyz = h['LeafXYZ/data'][:]; V = h['LeafSize/data'][:]**3
        r = np.sqrt((xyz**2).sum(axis=0))
        f = {'H0': h['x_HI/data'][:], 'N0': h['x_n_stages/data'][:, 0],
             'N1': h['x_n_stages/data'][:, 1], 'N2': h['x_n_stages/data'][:, 2],
             'T': h['T_e/data'][:]}
    idx = np.digitize(r, edges) - 1; out = {}
    for k, v in f.items():
        a = np.full(edges.size - 1, np.nan)
        for i in range(a.size):
            s = idx == i
            if s.any(): a[i] = np.average(v[s], weights=V[s])
        out[k] = a
    return out
def cloudy(tag):
    n = np.genfromtxt(f'{CL}/{tag}_c25.ele_n', names=True)
    h = np.genfromtxt(f'{CL}/{tag}_c25.ele_h', names=True)
    rn = (R0 + n['depth'])/PC; rh = (R0 + h['depth'])/PC
    return dict(rN=rn, N0=n['N'], N1=n['N_1'] if 'N_1' in n.dtype.names else n[n.dtype.names[2]],
                N2=n[n.dtype.names[3]], rH=rh, H0=h['H'])
LAB = [('n', 2, 6585.27, 6583.45, '[N II] 6584'), ('n', 2, 6549.86, 6548.05, '[N II] 6548'),
       ('n', 2, 5756.19, 5754.59, '[N II] 5755'), ('n', 2, 1218026.8, None, '[N II] 122um'),
       ('n', 1, 5199.35, None, '[N I] 5198+5200'), ('n', 3, 573238., None, '[N III] 57um')]
fig, axes = plt.subplots(2, 2, figsize=(10, 7))
for col, tag in enumerate(('hii40', 'hii20')):
    old = f'../m2_atomic/{tag}_atomic'; new = f'{tag}_nitrogen'
    print(f'=== {tag.upper()} ===')
    To, Tn = t_npne(old + '_rates.h5'), t_npne(new + '_rates.h5')
    Lo, Ln = lines(old + '_lines.txt'), lines(new + '_lines.txt')
    print(f'  <T[NpNe]> {To:.2f} -> {Tn:.2f} K ({Tn-To:+.2f})   L(Hb) {Lo["Hb"]:.5e} -> {Ln["Hb"]:.5e} ({100*(Ln["Hb"]/Lo["Hb"]-1):+.3f}%)')
    rel = {}
    for ln in open(f'{CL}/{tag}_c25_rel.lin'):
        p = ln.rstrip('\n').split('\t')
        if len(p) >= 2 and p[0].strip():
            rel[p[0].strip()] = float(p[-1])
    for el, st, lam, clam, name in LAB:
        o = pick(Lo, el, st, lam); n = pick(Ln, el, st, lam)
        if name == '[N I] 5198+5200': o += pick(Lo, 'n', 1, 5201.71); n += pick(Ln, 'n', 1, 5201.71)
        c = np.nan
        if clam:
            key = [k for k in rel if k.lower().startswith('n  2') and abs(float(k.split()[2][:-1]) - clam) < 0.5]
            if key: c = rel[key[0]]
        print(f'  {name:16s} {o:.4e} -> {n:.4e} rel Hb ({100*(n/o-1):+7.2f}%)   Cloudy {c:.4e}')
    edges = np.arange(0.0, 6.0, 0.02); rc = 0.5*(edges[1:] + edges[:-1])
    po, pn = profiles(old + '_rates.h5', edges), profiles(new + '_rates.h5', edges)
    c = cloudy(tag)
    ax = axes[0, col]
    for k, cc in (('N0', 'C0'), ('N1', 'C1'), ('N2', 'C2')):
        ax.semilogy(rc, po[k], c=cc, ls='--', lw=1)
        ax.semilogy(rc, pn[k], c=cc, ls='-', lw=1.5, label={'N0': r'N$^0$', 'N1': r'N$^+$', 'N2': r'N$^{2+}$'}[k])
        ax.semilogy(c['rN'], c[k], c=cc, ls=':', lw=1.5)
    ax.semilogy(rc, pn['H0'], c='k', lw=0.8, label=r'H$^0$')
    ax.semilogy(c['rH'], c['H0'], c='k', ls=':', lw=0.8)
    ax.set_ylim(1e-4, 1.5); ax.set_title(tag.upper() + r': solid new, dashed old, dotted Cloudy c25')
    ax.set_xlabel('r [pc]'); ax.set_ylabel('ion fraction'); ax.legend(fontsize=8, loc='lower left')
    fr = pn['H0'] > 0.01
    lo_ = rc[np.nanargmax(np.where(pn['H0'] > 0.02, 1, 0))] - 0.4 if fr.any() else 0
    ax.set_xlim(max(0.8*np.nanmin(c['rN']), rc[np.where(pn['H0'] > 1e-3)[0][0]] - 0.6), min(6, c['rN'].max()*1.02))
    ax2 = axes[1, col]
    ax2.semilogy(rc, po['N0']/po['H0'], c='C0', ls='--', label='old (Huang rows)')
    ax2.semilogy(rc, pn['N0']/pn['H0'], c='C0', ls='-', label='new (KF96)')
    hc = np.interp(c['rN'], c['rH'], c['H0'])
    ax2.semilogy(c['rN'], c['N0']/hc, c='C0', ls=':', label='Cloudy c25')
    ax2.set_xlim(ax.get_xlim()); ax2.set_ylim(1e-2, 1e2)
    ax2.set_xlabel('r [pc]'); ax2.set_ylabel(r'$x({\rm N^0})/x({\rm H^0})$')
    ax2.legend(fontsize=8)
fig.tight_layout(); fig.savefig('nitrogen_profiles.pdf'); fig.savefig('nitrogen_profiles.png', dpi=80)
