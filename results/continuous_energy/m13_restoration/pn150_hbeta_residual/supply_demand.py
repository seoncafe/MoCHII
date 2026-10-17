#!/usr/bin/env python3
"""Fixed-structure test (2026-10-17): on a written state, the ionizations the
transported field supplies against the recombinations the state demands,
per stellar Q(H) = 5.416e47 s^-1, for H0 -> H+, He0 -> He+, He+ -> He2+:
supply  = sum n(X) [Gamma_X + n_e C_X(T)] V   (photo + collisional, Voronov 1997)
demand  = sum n_e n(X+) alpha_A(T) V          (tools/python/mochii_rates.py)
Charge exchange is not included (Cloudy without it moves R_B(H)/Q by 0.7%).
On a converged MoCHII state supply/demand is the calibration of the test;
on Cloudy's structure (make_cloudy_state.py, one MoCHII pass with the state
held) it says whether MoCHII's field sustains Cloudy's ionization.  Also by
radius bins.  Usage: supply_demand.py <rates.h5> ..."""
import sys, numpy as np, h5py
sys.path.insert(0, '../../../../tools/python'); import mochii_rates as mr
Q, nH, yHe = 5.416e47, 3000.0, 0.1
VOR = dict(H=(13.6, 0, 0.291e-7, 0.232, 0.39), He=(24.6, 0, 0.175e-7, 0.180, 0.35),
           Hep=(54.4, 1, 0.205e-8, 0.265, 0.25))
lgT = np.linspace(3.0, 5.0, 401)
def tab(fn): t = np.array([fn(10**x) for x in lgT]); return lambda T: np.interp(np.log10(T), lgT, t)
aH, aHe, aHep = tab(mr.alphaA_HII), tab(mr.alphaA_HeII), tab(mr.alphaA_HeIII)
cH, cHe, cHep = (tab(lambda T, v=v: mr.voronov(T, *v)) for v in (VOR['H'], VOR['He'], VOR['Hep']))
for f in sys.argv[1:]:
    with h5py.File(f, 'r') as h:
        g = lambda k: h[k]['data'][()]
        d = h['Gamma_HI'].attrs['DIST_CM'][0]
        V, ne = (g('LeafSize')*d)**3, g('n_e'); m = ne > 0
        r = np.sqrt((g('LeafXYZ')**2).sum(axis=0))[m]*d/3.0857e18
        T, ne, V = g('T_e')[m], ne[m], V[m]
        x = {k: g(k)[m] for k in ('x_HI', 'x_HII', 'x_HeI', 'x_HeII', 'x_HeIII')}
        G = {k: g(k)[m] for k in ('Gamma_HI', 'Gamma_HeI', 'Gamma_HeII')}
    rows = [('H0 -> H+', nH*x['x_HI'], G['Gamma_HI'], cH(T), nH*x['x_HII'], aH(T)),
            ('He0 -> He+', nH*yHe*x['x_HeI'], G['Gamma_HeI'], cHe(T), nH*yHe*x['x_HeII'], aHe(T)),
            ('He+ -> He2+', nH*yHe*x['x_HeII'], G['Gamma_HeII'], cHep(T), nH*yHe*x['x_HeIII'], aHep(T))]
    print(f)
    edges = np.array([0.03, 0.06, 0.08, 0.09, 0.10, 0.11, 0.12, 0.125, 0.13, 0.135, 0.16])
    for name, n0, Gm, C, n1, a in rows:
        sp = n0*Gm*V; sc = n0*ne*C*V; dm = ne*n1*a*V
        print(f'  {name:12s} photo {sp.sum()/Q:.4f}  coll {sc.sum()/Q:.4f}  demand {dm.sum()/Q:.4f}  '
              f'supply/demand {(sp.sum() + sc.sum())/dm.sum():.4f}')
        b = np.digitize(r, edges)
        print('     r [pc]   ' + ' '.join(f'{edges[i-1]:.3f}-{edges[i]:.3f}' for i in range(1, len(edges))))
        print('     s/d      ' + ' '.join(f'{(sp[b == i].sum() + sc[b == i].sum())/max(dm[b == i].sum(), 1e-300):11.4f}' for i in range(1, len(edges))))
        print('     d/Q      ' + ' '.join(f'{dm[b == i].sum()/Q:11.4f}' for i in range(1, len(edges))))
