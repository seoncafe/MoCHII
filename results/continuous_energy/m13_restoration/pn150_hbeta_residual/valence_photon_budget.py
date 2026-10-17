#!/usr/bin/env python3
"""Photons of Cloudy's heavy-ion recombination approximation (one at 0.75 and
one at 0.25 of the ionization potential of the recombined ion, for every
RR + DR recombination whose product keeps >= 3 electrons) per stellar Q(H),
evaluated the same way on a MoCHII state and on Cloudy's own PN150 structure:
sum over cells of n_e n(X^(i+1)) alpha_RR+DR(T) V, with alpha from
data/atomic/element_<el>.txt (the Badnell totals both codes ship), counting
the photons above 13.598 eV, above 24.587 eV and above 54.418 eV, and the
total recombinations.  Cloudy: gam.ovr (T_e, n_e, depth) and gam.ele_<el>
(stage fractions) of gam.in, zones as spherical shells from 1e17 cm.
Usage: valence_photon_budget.py [MoCHII rates.h5 ...]"""
import sys, numpy as np, h5py
sys.path.insert(0, 'vf_check')
from channel8_check import read_element, alpha, ABUND
Q, nH, EV = 5.416e47, 3000.0, 1.602176634e-12
EDGES = (13.598, 24.587, 54.418)

def tally(get_x, T, ne, V):
    """get_x(el) -> array (ncell, nstage) of stage fractions."""
    out = dict(rec=0.0, n=np.zeros(3), by_el={})
    for el, ab in ABUND.items():
        Z, tr = read_element(el)
        x = get_x(el)
        nel = np.zeros(3)
        for it, t in enumerate(tr, start=1):
            if Z - it + 1 < 3 or it >= x.shape[1]:
                continue
            rate = np.sum(ne*nH*ab*x[:, it]*alpha(t, T)*V)
            out['rec'] += rate
            for f in (0.75, 0.25):
                e = f*t['eth']
                nel += rate*np.array([e >= E for E in EDGES])
        out['n'] += nel
        out['by_el'][el] = nel[0]
    return out

def report(name, r):
    print(f"{name:14s} recombinations {r['rec']/Q:.4f} Q; photons > 13.598 / 24.587 / "
          f"54.418 eV: {r['n'][0]/Q:.4f} {r['n'][1]/Q:.4f} {r['n'][2]/Q:.4f} Q;  "
          + ' '.join(f"{el} {v/Q:.4f}" for el, v in r['by_el'].items()))

o = np.loadtxt('gam.ovr'); dep, T, ne = o[:, 0], o[:, 1], o[:, 4]
r = 1e17 + dep
edge = np.concatenate([[1e17], 0.5*(r[1:] + r[:-1]), [r[-1] + 0.5*(r[-1] - r[-2])]])
V = 4.0/3.0*np.pi*(edge[1:]**3 - edge[:-1]**3)
def cloudy_x(el):
    a = np.loadtxt(f'gam.ele_{el}')
    Z = read_element(el)[0]
    return a[:, 1:Z + 2]          # stage fractions X^0 .. X^Z
report('Cloudy c25', tally(cloudy_x, T, ne, V))

for f in sys.argv[1:]:
    with h5py.File(f, 'r') as h:
        g = lambda k: h[k]['data'][()]
        d = h['Gamma_HI'].attrs['DIST_CM'][0]
        Vm, Tm, nem = (g('LeafSize')*d)**3, g('T_e'), g('n_e')
        m = nem > 0
        xs = {el: g(f'x_{el}_stages')[m] for el in ABUND}
    report(f.replace('_rates.h5', ''), tally(lambda el: xs[el], Tm[m], nem[m], Vm[m]))
