"""Cloudy c25.00 heavy-element recombination photons (ion_recomb.cpp:91-95,
cont_createpointers.cpp:1062-1066): every RR+DR recombination emits one photon
at 0.75 IP and one at 0.25 IP. Ionizing (>= 13.598 eV) part on the same leaves."""
import json, numpy as np, leaf_luminosities as LL
rc = json.load(open('recomb_coefficients.json')); Tg = np.array(rc['T'])
for case in LL.CASES:
    leaf, log = LL.load(case); on = leaf['nH'] > 0
    T, ne, nH, V = leaf['T'][on], leaf['ne'][on], leaf['nH'][on], leaf['V'][on]
    tot = 0.0
    for el, ab in LL.ABUND.items():
        x = leaf['x_' + el][on]
        for st in range(2, 8):
            r = rc.get(f'{el}{st-1}')
            if r is None: continue
            a = np.interp(np.log10(T), np.log10(Tg), np.array(r['alpha_rr']) + np.array(r['alpha_dr']))
            ip = r['ip']; e = (0.75 * ip if 0.75 * ip >= 13.598 else 0) + (0.25 * ip if 0.25 * ip >= 13.598 else 0)
            tot += np.sum(V * ne * ab * nH * x[:, st - 1] * a) * e * LL.EV
    print(case, 'Cloudy-style 0.75/0.25 IP pair, ionizing part [erg/s]: %.3e' % tot)
