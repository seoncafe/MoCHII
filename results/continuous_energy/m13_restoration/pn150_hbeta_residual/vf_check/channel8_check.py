#!/usr/bin/env python3
"""Independent evaluation of diffuse channel 8 (par%metal_recombination_photons
= 'valence_fractions') on a saved PN150 state: for every metal transition
whose recombined ion keeps >= 3 electrons, n_e n(X^(i+1)) alpha_RR+DR(T) V
times 0.75 and 0.25 of the ionization potential, split into in band
[eion_min, eion_max) and outside.  The rates are read from
data/atomic/element_<el>.txt (RR, DR rows) and evaluated as species_mod
alpha_rec does.  Usage: channel8_check.py <rates.h5> [eion_min eion_max]"""
import os, sys, numpy as np, h5py
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    '../../../../../data/atomic/')
ABUND = dict(c=2.999163e-4, n=1.0e-4, o=5.997911e-4, ne=1.499685e-4,
             mg=2.999163e-5, si=2.999163e-5, s=1.499685e-5)
EV = 1.602176634e-12

def read_element(el):
    Z, tr, cur = None, [], None
    for line in open(ROOT + f'element_{el}.txt'):
        w = line.split()
        if not w or w[0].startswith('#'):
            continue
        if w[0] == 'ELEMENT':
            Z = int(w[2])
        elif w[0] == 'TRANSITION':
            cur = dict(rr=None, rr2=None, dr=[], dr2=None); tr.append(cur)
        elif w[0] == 'ETH':
            cur['eth'] = float(w[1])
        elif w[0] == 'RR':
            cur['rr'] = [float(x) for x in w[1:7]]
        elif w[0] == 'RR2':
            cur['rr2'] = [float(x) for x in w[1:3]]
        elif w[0] == 'DR':
            n = int(w[1]); v = [float(x) for x in w[2:2 + 2*n]]
            cur['dr'] = list(zip(v[:n], v[n:]))
        elif w[0] == 'DR2':
            cur['dr2'] = [float(x) for x in w[1:5]]
    return Z, tr

def alpha(t, T):
    if t['rr2'] is not None:
        a = t['rr2'][0]*(T/1e4)**(-t['rr2'][1])
    else:
        A, B, T0, T1, C, T2 = t['rr']
        s0, s1 = np.sqrt(T/T0), np.sqrt(T/T1); b = B + C*np.exp(-T2/T)
        a = A/(s0*(1 + s0)**(1 - b)*(1 + s1)**(1 + b))
    for c, e in t['dr']:
        a = a + T**-1.5*c*np.exp(-e/T)
    if t['dr2'] is not None:
        A, B, T0, T1 = t['dr2']
        a = a + A*T**-1.5*np.exp(-T0/T)*(1 + B*np.exp(-T1/T))
    return a

def main():
    emin, emax = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (13.598, 200.0)
    with h5py.File(sys.argv[1], 'r') as h:
        g = lambda k: h[k]['data'][()]
        d = h['Gamma_HI'].attrs['DIST_CM'][0]
        V, T, ne, nH = (g('LeafSize')*d)**3, g('T_e'), g('n_e'), 3000.0
        m = ne > 0
        inb = out_lo = out_hi = 0.0
        for el, ab in ABUND.items():
            Z, tr = read_element(el)
            x = g(f'x_{el}_stages')
            for it, t in enumerate(tr, start=1):
                if Z - it + 1 < 3:
                    continue
                rate = ne[m]*nH*ab*x[m, it]*alpha(t, T[m])*V[m]
                for f in (0.75, 0.25):
                    e = f*t['eth']; L = np.sum(rate)*e*EV
                    if e < emin:
                        out_lo += L
                    elif e >= emax:
                        out_hi += L
                    else:
                        inb += L
    print(f'channel 8 in band {inb:.6e} erg/s, below eion_min {out_lo:.6e}, at/above eion_max {out_hi:.6e}')

if __name__ == '__main__':
    main()
