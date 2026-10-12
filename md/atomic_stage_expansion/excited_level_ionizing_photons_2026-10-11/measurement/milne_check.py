"""Ground-term RR rate: Milne on the VFKY96 outer-shell fit (channel 6 data)
against the Badnell partial into the ground term (recomb_coefficients.json)."""
import json, numpy as np
SNAP = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase4_20261011/data_snapshot/atomic'
g = {}
for L in open(f'{SNAP}/ground_term_weights.txt'):
    t = L.split()
    if t and t[0] == 'GROUND': g[(t[1], int(t[3]))] = int(t[4])
def photo(el):
    out = {}; cur = None
    for L in open(f'{SNAP}/element_{el}.txt'):
        t = L.split()
        if not t: continue
        if t[0] == 'TRANSITION': cur = int(t[1])
        if t[0] == 'PHOTO': out[cur] = [float(x) for x in t[1:]]
        if t[0] == 'PHINNER': out[(cur, 'einn')] = float(t[1])
    return out
def vfky(E, p):
    Eth, E0, s0, ya, P, yw, y0, y1 = p
    x = E / E0 - y0; y = np.sqrt(x * x + y1 * y1)
    F = ((x - 1) ** 2 + yw ** 2) * y ** (0.5 * P - 5.5) * (1 + np.sqrt(y / ya)) ** (-P)
    return np.where(E >= Eth, s0 * F * 1e-18, 0.0)
me, c, k, eV = 9.1093837e-28, 2.99792458e10, 1.380649e-16, 1.602176634e-12
rc = json.load(open('recomb_coefficients.json')); Tg = np.array(rc['T'])
for T in (1.0e4, 1.6e4):
    print('T =', T)
    for el in ['c', 'n', 'o', 'ne', 's']:
        P = photo(el)
        for q in range(1, 7):
            p = P[q]; Eth = p[0]; kT = k * T / eV
            E = Eth + np.linspace(0, 40 * kT, 20001)
            sig = vfky(E, p)
            gl, gu = g[(el, q)], g[(el, q + 1)]
            a1 = np.sqrt(2 / np.pi) / (c * c * (me * k * T) ** 1.5) * gl / gu * \
                 np.trapz((E * eV) ** 2 * sig * np.exp(-(E - Eth) / kT), E * eV)
            r = rc[f'{el}{q}']; i = np.argmin(abs(Tg - T))
            # Badnell ground-term partial rate: rr_gt holds rate x mean photon energy; divide by mean energy ~ IP + kT
            agt = r["rr_gt_rate"][i] / r["alpha_rr_partial"][i] * r["alpha_rr"][i]
            print(f'  {el}{q}  Milne(VFKY96) {a1:.3e}  Badnell gt {agt:.3e}  ratio {a1/agt if agt>0 else np.nan:.3f}')
