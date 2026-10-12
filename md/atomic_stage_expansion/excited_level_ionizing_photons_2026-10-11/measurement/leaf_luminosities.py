"""Leaf-integrated luminosities above 13.598 eV for the two Phase 4 runs:
(3) collisionally excited metal lines (MoCHII nlevel solve, CHIANTI atoms),
(1)/(2) metal recombination photons (recomb_coefficients.json).
"""
import json, os, sys
import numpy as np, h5py
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase4_20261011/runs/work'
EUV = os.environ.get('EUV_DIR', os.path.join(HERE, 'euv'))
PC = 3.0856775814913673e18
EV = 1.602176634e-12
ABUND = dict(c=2.2e-4, n=4.0e-5, o=3.3e-4, ne=5.0e-5, s=9.0e-6)
CASES = {
    'hii40': ('hii40/hii40_p4_rates.h5', '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/g2_hii/hii40_L6.fits',
              dict(total=5.8916e38, ch=[5.668e38, 9.430e36, 1.060e35, 1.214e37, 2.056e35, 4.501e35, 2.669e29])),
    'bb1e5': ('bb1e5/bb1e5_p4_rates.h5', '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/g2_hii/hii40_L6_rout5p2.fits',
              dict(total=1.1781e39, ch=[7.567e38, 1.171e38, 5.394e37, 1.476e38, 1.008e38, 2.003e36, 2.444e31])),
    'hii40_m11': ('/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m11_hecx_solver/hii40_converged_rates.h5',
                  '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/g2_hii/hii40_L6_rout5p2.fits',
                  dict(total=7.2451e38, ch=[0, 0, 0, 0, 0, 0, 0])),
}


def load(case):
    h5, ft, logL = CASES[case]
    f = h5py.File(os.path.join(RUNS, h5), 'r')
    xyz = f['LeafXYZ/data'][:].T
    s = f['LeafSize/data'][:]
    t = fits.open(ft)[1].data
    key = lambda a: np.round(a / s[0] * 2).astype(np.int64)
    idx = {(a, b, c): i for i, (a, b, c) in enumerate(zip(key(t['x']), key(t['y']), key(t['z'])))}
    k = key(xyz)
    order = np.array([idx[(a, b, c)] for a, b, c in k])
    nH = np.asarray(t['nH'], float)[order]
    leaf = dict(T=f['T_e/data'][:], ne=f['n_e/data'][:], nH=nH, V=(s * PC) ** 3)
    for el in ABUND:
        xs = f[f'x_{el}_stages/data'][:]
        pad = np.zeros((xs.shape[0], 7)); pad[:, :xs.shape[1]] = xs
        leaf['x_' + el] = pad
    return leaf, logL


def euv_table(el, stage):
    p = os.path.join(EUV, f'{el}_{stage}.dat')
    a = np.loadtxt(p)
    lt = np.unique(a[:, 0]); ln = np.unique(a[:, 1])
    eps = a[:, 2].reshape(len(lt), len(ln)); he = (a[:, 5] + a[:, 6]).reshape(len(lt), len(ln))
    return lt, ln, eps, he


def interp2(lt, ln, tab, T, ne):
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((lt, ln), np.log10(np.maximum(tab, 1e-300)), bounds_error=False, fill_value=None)
    x = np.column_stack([np.clip(np.log10(T), lt[0], lt[-1]), np.clip(np.log10(np.maximum(ne, 1.0)), ln[0], ln[-1])])
    return 10 ** f(x)


def main():
    rc = json.load(open(os.path.join(HERE, os.environ.get('RECOMB_JSON', 'recomb_coefficients.json'))))
    Tg = np.array(rc['T'])
    report = {}
    for case in CASES:
        leaf, logL = load(case)
        on = leaf['nH'] > 0
        T, ne, nH, V = leaf['T'][on], leaf['ne'][on], leaf['nH'][on], leaf['V'][on]
        rows = []
        for el, ab in ABUND.items():
            x = leaf['x_' + el][on]
            for st in range(1, 8):
                nion = ab * nH * x[:, st - 1]
                # (3) collisional lines of stage st
                L3 = L3he = 0.0
                if os.path.exists(os.path.join(EUV, f'{el}_{st}.dat')):
                    lt, ln, eps, he = euv_table(el, st)
                    L3 = np.sum(V * ne * nion * interp2(lt, ln, eps, T, ne))
                    L3he = np.sum(V * ne * nion * interp2(lt, ln, he, T, ne))
                # (1)/(2): recombination of stage st (charge q = st - 1) into stage st - 1
                q = st - 1
                r = rc.get(f'{el}{q}') if q >= 1 else None
                rec = {}
                if r is not None:
                    def co(name, kind):
                        part = np.array(r[f'alpha_{kind}_partial'])
                        fit = np.array(r[f'alpha_{kind}'])
                        per = np.array(r[name]) / np.maximum(part, 1e-300) * fit   # eV cm^3 s^-1, scaled to the MoCHII rate
                        return np.interp(np.log10(T), np.log10(Tg), per)
                    w = V * ne * nion * EV
                    for nm, kind in [('rr_gt', 'rr'), ('rr_fb_exc', 'rr'), ('rr_fb_exc_he', 'rr'), ('rr_cas', 'rr'),
                                     ('rr_cas_he', 'rr'), ('rr_tot_energy', 'rr'), ('dr_stab', 'dr'), ('dr_stab_he', 'dr'),
                                     ('dr_cas', 'dr'), ('dr_cas_he', 'dr'), ('dr_tot_energy', 'dr')]:
                        rec[nm] = float(np.sum(w * co(nm, kind)))
                    a_rr = np.interp(np.log10(T), np.log10(Tg), r['alpha_rr'])
                    a_dr = np.interp(np.log10(T), np.log10(Tg), r['alpha_dr'])
                    kT = 8.617333e-5 * T
                    rec['bound_ip_kT'] = float(np.sum(w * (a_rr + a_dr) * (r['ip'] + kT)))
                rows.append(dict(ion=f'{el}_{st}', L3=float(L3), L3he=float(L3he), **rec))
        report[case] = dict(rows=rows, log=logL)
    json.dump(report, open(os.path.join(HERE, os.environ.get('LEAF_JSON', 'leaf_luminosities.json')), 'w'), indent=1)
    for case, R in report.items():
        tot, ch6 = R['log']['total'], R['log']['ch'][5]
        print(f'\n== {case}: diffuse total {tot:.4e}, channel 6 {ch6:.4e} erg/s')
        print(f'{"ion":7s} {"(3) coll":>10s} {"gt(ch6)":>10s} {"(1) fbexc":>10s} {"(2) casRR":>10s} {"(2) stab":>10s} {"(2) casDR":>10s} {"recEtot":>10s} {"bound":>10s}')
        S = {}
        for r in R['rows']:
            for k, v in r.items():
                if k != 'ion':
                    S[k] = S.get(k, 0) + v
            if r.get('rr_tot_energy', 0) > 0 or r['L3'] > 0:
                print(f"{r['ion']:7s} {r['L3']:10.3e} {r.get('rr_gt',0):10.3e} {r.get('rr_fb_exc',0):10.3e} {r.get('rr_cas',0):10.3e} "
                      f"{r.get('dr_stab',0):10.3e} {r.get('dr_cas',0):10.3e} {r.get('rr_tot_energy',0)+r.get('dr_tot_energy',0):10.3e} {r.get('bound_ip_kT',0):10.3e}")
        print('SUM    ', ' '.join(f'{k}={v:.3e}' for k, v in S.items()))


if __name__ == '__main__':
    main()
