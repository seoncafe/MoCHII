#!/usr/bin/env python3
"""Share of free-free absorption in the ionizing-band opacity of a MoCHII run.

kappa_ff(nu) = n_e sum_Z n_Z FF_SPECTRAL Z^2 g_ff(nu,T,Z) T^-1/2 c^2 (1 - e^-u)
               / (8 pi h nu^3),   u = h nu/kT
(Kirchhoff's law on the code's own emission coefficient: the absorption
coefficient net of stimulated emission), with g_ff from the code
(freefree_gaunt_table.txt, written by tabulate_freefree_gaunt.x).
Against the bound-free opacity of H I, He I, He II alone (metals and dust
left out, so the share is an upper bound).  n_H of a leaf is taken from n_e
as n_e/(x_HII + y_He (x_HeII + 2 x_HeIII)) -- the metal electrons left out,
which overstates n_H and so the share again.  Z = 1: H+ and He+, Z = 2:
He2+ (metal ions left out of the free-free: abundance x Z^2 <= 1e-3 of it).

usage: measure_freefree_absorption.py rates.h5 eion_max [abund_He]
"""
import sys
import numpy as np
import h5py

FF_SPECTRAL = 6.841593e-38
H = 6.62607015e-27; C = 2.99792458e10; K = 1.380649e-16; EV = 1.602176634e-12

def sigma_hydrogenic(E, Z):
    E1 = 13.605693 * Z * Z * (0.999466 if Z == 1 else 0.999864)
    s = np.zeros_like(E)
    m = E >= E1
    x = E[m] / E1
    eps = np.sqrt(np.maximum(x - 1.0, 1e-30))
    s[m] = 6.304e-18 / Z**2 * x**-4 * np.exp(4.0 - 4.0 * np.arctan(eps) / eps) \
        / (1.0 - np.exp(-2.0 * np.pi / eps))
    return s

def sigma_hei(E):
    # VFKY96 He I (phfit2 PH2): E0 13.61, s0 949.2 Mb, ya 1.469, P 3.188, yw 2.039, y0 .4434, y1 2.136
    E0, s0, ya, P, yw, y0, y1 = 13.61, 949.2, 1.469, 3.188, 2.039, 0.4434, 2.136
    x = E / E0 - y0
    y = np.sqrt(x * x + y1 * y1)
    f = ((x - 1) ** 2 + yw * yw) * y ** (0.5 * P - 5.5) * (1 + np.sqrt(y / ya)) ** (-P)
    return np.where(E >= 24.587, s0 * f * 1e-18, 0.0)

def main():
    fn = sys.argv[1]; emax = float(sys.argv[2])
    yhe = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
    here = __file__.rsplit('/', 1)[0]
    tab = np.loadtxt(here + '/freefree_gaunt_table.txt')
    lt_nodes = np.unique(tab[:, 0]); e_nodes = np.unique(tab[:, 1])
    g = tab[:, 2:].reshape(len(lt_nodes), len(e_nodes), 6)
    f = h5py.File(fn, 'r')
    T = f['T_e/data'][:]; ne = f['n_e/data'][:]
    xHI = f['x_HI/data'][:]; xHII = f['x_HII/data'][:]
    xHeI = f['x_HeI/data'][:]; xHeII = f['x_HeII/data'][:]; xHeIII = f['x_HeIII/data'][:]
    ok = ne > 0
    T, ne, xHI, xHII, xHeI, xHeII, xHeIII = (a[ok] for a in (T, ne, xHI, xHII, xHeI, xHeII, xHeIII))
    nH = ne / (xHII + yhe * (xHeII + 2 * xHeIII))
    E = np.geomspace(13.60, emax, 400)  # above the H edge (13.5984 eV)
    lt = np.clip(np.log10(T), lt_nodes[0], lt_nodes[-1])
    it = np.clip(np.searchsorted(lt_nodes, lt) - 1, 0, len(lt_nodes) - 2)
    w = (lt - lt_nodes[it]) / (lt_nodes[1] - lt_nodes[0])
    worst = 0.0; where = None; med = []
    for e in E:
        ie = np.clip(np.searchsorted(e_nodes, e) - 1, 0, len(e_nodes) - 2)
        we = (np.log(e) - np.log(e_nodes[ie])) / (np.log(e_nodes[ie + 1]) - np.log(e_nodes[ie]))
        def gz(z):
            g00 = g[it, ie, z - 1]; g10 = g[it + 1, ie, z - 1]
            g01 = g[it, ie + 1, z - 1]; g11 = g[it + 1, ie + 1, z - 1]
            return (1 - w) * (1 - we) * g00 + w * (1 - we) * g10 + (1 - w) * we * g01 + w * we * g11
        nu = e * EV / H; u = e * EV / (K * T)
        pref = FF_SPECTRAL * C * C * (-np.expm1(-u)) / (8 * np.pi * H * nu ** 3) / np.sqrt(T)
        kff = ne * pref * (nH * (xHII + yhe * xHeII) * gz(1) + 4 * nH * yhe * xHeIII * gz(2))
        kbf = nH * (xHI * sigma_hydrogenic(np.array([e]), 1)[0] + yhe * xHeI * sigma_hei(np.array([e]))
                    + yhe * xHeII * sigma_hydrogenic(np.array([e]), 2)[0])
        share = kff / (kff + kbf)
        i = np.argmax(share)
        med.append(np.median(share))
        if share[i] > worst:
            worst = share[i]; where = (e, T[i], ne[i], xHI[i], xHeI[i], xHeII[i], xHeIII[i])
    print(f'{fn}: leaves with n_e > 0: {ok.sum()}')
    print(f'  largest kappa_ff/(kappa_ff + kappa_bf(H, He)) = {worst:.3e}')
    print('  at E = %.2f eV, T = %.0f K, n_e = %.3g, x_HI = %.2e, x_HeI = %.2e, x_HeII = %.2e, x_HeIII = %.2e' % where)
    print(f'  largest leaf-median share over the band = {max(med):.3e}')

if __name__ == '__main__':
    main()
