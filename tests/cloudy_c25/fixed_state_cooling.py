#!/usr/bin/env python3
"""Fixed-state cooling comparison of MoCHII against Cloudy c25.00 (HII40).

What it computes
----------------
Cloudy's converged HII40 model (tests/cloudy_c25/hii40_cool.*, the unmodified
physics block of tsuite/auto/hii_paris.in with output commands appended) is
read zone by zone: T_e, n_e, n_H and the H/He stages from the overview
(.ovr), the metal stages from the element files (.ele_c/n/o/ne/s), the
abundances from the composition block of the .out file, and the zone radius
and thickness from the .rad file.  MoCHII's cooling terms are evaluated on
exactly that state and both sides are integrated over the sphere,
dV = 4 pi r^2 dr, channel by channel.  No thermal or ionization feedback
enters either side, so a ratio measures the atomic data and the cooling
expressions alone.

Cloudy side: the .cool file lists, for every zone, the fraction of the total
cooling Ctot carried by each coolant above 1e-4 of the total (set WeakHeatCool
-4); a channel is the sum of fraction x Ctot x dV.  Labels (sources in
~/CLOUDY/c25.00/source): 'X  n' the collisionally excited lines of ion
X^(n-1)+; 'ISclinH H' / 'ISclinHeHe' the H I / He I collisional lines
(iso_cool.cpp:567); 'ISrcolH H' / 'ISrcolHeHe' / 'ISrcolH He' the free-bound
recombination cooling onto H+ / He+ / He++, net of the photoionization heating
of the excited levels (iso_cool.cpp:317); 'FF c' net free-free
(cool_eval.cpp:466); 'IScionH H' net collisional ionization of H; 'hvFB'
heavy-element free-bound (cool_eval.cpp:514).

MoCHII side, zone by zone (n_ion from Cloudy's fractions):
  metal lines  from the compiled MoCHII modules through the driver
               fixed_state_nlevel.f90 (C I-III, N I-III, O I-III, Ne II-III,
               S II-IV):
                 fit      n_e n_ion lam_fit      Tier-1 zero-density fit
                                                 (cooling_model='low_density')
                 localne  n_e n_ion lam_localne  the (T, n_e) table of the
                                                 default cooling_model='local_ne'
                 nlev_e   n_ion P_e              n-level emissivities at the
                                                 local n_e, electron impact
                 nlev_eH  n_ion P_eH             the same with atomic-hydrogen
                                                 impact (C I, C II, O I) in one
                                                 statistical equilibrium,
                                                 n_HI = n_H x(H0)
               and two atomic-hydrogen impact terms added to the thermal loop:
                 hjoint   n_ion (P_eH - P_e)     species_mod::metal_cooling_H
                                                 since 2026-09-24
                 h2lv     two-level [O I] 63 um / [C II] 158 um H impact,
                                                 metal_cooling_H before
                                                 2026-09-24
  H I lines    n_e n_HI lam_fit(h_1) or lam_localne(h_1)
  recombination n_e [n(H+) b_HII + n(He+) b_HeII + n(He++) b_HeIII]:
                 current   b = k T [(3/2) alpha + d alpha/d ln T] on the
                           balance alpha, case A or B (cooling_mod since
                           2026-09-24; tools/python/mochii_rates.py)
                 previous  Hui & Gnedin (1997) fits for H II and He III, k T
                           alpha for He II (cooling_mod before 2026-09-24,
                           reproduced below), with today's alpha
                 jul25     'previous' with the alpha_B(He I) of 2026-07-25,
                           when the comparison was first made: Badnell alpha_A
                           minus the Mao & Kaastra (2016) alpha_1 in place of
                           the Milne-relation alpha_1 (case B only; case A
                           is unaffected)
  free-free    1.42554e-27 sqrt(T) n_e sum n_ion Z^2 gbar(T, Z) over H/He and
               the metal ions (par%metal_freefree = .true., the default)
  coll. ioniz. n_e sum n_i C_i(T) E_th,i for H I, He I, He II
Not evaluated: the metal recombination and collisional-ionization cooling
(species_mod::metal_ionization_cooling, added 2026-09-24; Cloudy's 'hvFB' is
printed for scale), Cl lines (Cloudy's Cl stages are not saved), He I
collisional lines (MoCHII carries none).

Temperature scaling: T_e of every zone is multiplied by f, with n_e and the
fractions held, until the MoCHII total of every channel above equals Cloudy's
volume-integrated heating; Te(NeNp) scales by the same f.  The paper
(mochii.tex, Section "Cooling-budget verification at a fixed state") quotes
this for the zero-density line cooling without naming the recombination case,
so every combination is printed; the recorded 7784 K is the zero-density
case-B entry scaled from Cloudy's printed 8210 K.  The totals are tabulated on
a grid of f in one driver call, the root is found on a cubic spline of
ln(total) in ln f, and the root is then checked by an exact evaluation
(printed residual).

Usage (from anywhere; `make` must have been run at the repository root):
    tests/cloudy_c25/run_fixed_state.sh            # build the driver, then run
    python3 tests/cloudy_c25/fixed_state_cooling.py [--no-scaling]
"""

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'tools', 'python'))
import mochii_rates as mr  # noqa: E402

BASE = os.path.join(HERE, 'hii40_cool')
DRIVER = os.path.join(HERE, 'fixed_state_nlevel.x')

#--- metal ions of fixed_state_nlevel.f90 (h_1 is the H I line cooling)
METAL_IONS = ['c_1', 'c_2', 'c_3', 'n_1', 'n_2', 'n_3', 'o_1', 'o_2', 'o_3',
              'ne_2', 'ne_3', 's_2', 's_3', 's_4']
HCOLL_IONS = ['c_1', 'c_2', 'o_1']          # data/atomic/hcoll_*.txt
ROMAN = {1: 'I', 2: 'II', 3: 'III', 4: 'IV'}
EL_SYM = {'c': 'C', 'n': 'N', 'o': 'O', 'ne': 'Ne', 's': 'S'}

ETH_HI, ETH_HEI, ETH_HEII = 13.598, 24.587, 54.416     # src/define.f90
FF_PREF = 1.42554e-27                                   # cooling_mod

#--- line-cooling models: (metal column, H I column, H-impact addition)
LINE_MODELS = {
    'zero_density':      ('fit', 'fit', None),
    'zero_density_h2lv': ('fit', 'fit', 'h2lv'),
    'local_ne':          ('localne', 'localne', None),
    'local_ne_h2lv':     ('localne', 'localne', 'h2lv'),     # loop before 09-24
    'local_ne_hjoint':   ('localne', 'localne', 'hjoint'),   # loop now
    'nlevel_e':          ('nlev_e', 'fit', None),
    'nlevel_eH':         ('nlev_eH', 'fit', None),
}
LINE_LABELS = {
    'zero_density':      'Tier-1 zero-density fits',
    'zero_density_h2lv': 'Tier-1 + two-level H impact (pre-09-24)',
    'local_ne':          'local_ne table, electrons only',
    'local_ne_h2lv':     'local_ne + two-level H (loop pre-09-24)',
    'local_ne_hjoint':   'local_ne + joint H impact (loop now)',
    'nlevel_e':          'n-level at local n_e, e only',
    'nlevel_eH':         'n-level at local n_e, e + H',
}
REC_MODELS = ['previous_A', 'previous_B', 'jul25_B', 'current_A', 'current_B']


# ---------------------------------------------------------------------------
# Cloudy output
# ---------------------------------------------------------------------------
def read_table(fname):
    """Numeric rows of a tab-separated Cloudy save file; header -> names."""
    names, rows = None, []
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                if names is None:
                    names = [x.strip() for x in line[1:].rstrip('\n').split('\t')]
                continue
            rows.append([float(x) for x in line.split()])
    return names, np.array(rows)


def read_cool(fname):
    """.cool: depth, Htot, Ctot and {label: fraction array}."""
    depth, htot, ctot, frac = [], [], [], []
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            depth.append(float(p[0]))
            htot.append(float(p[2]))
            ctot.append(float(p[3]))
            frac.append({p[i].strip(): float(p[i + 1])
                         for i in range(5, len(p) - 1, 2) if p[i].strip()})
    labels = sorted({k for d in frac for k in d})
    fr = {k: np.array([d.get(k, 0.0) for d in frac]) for k in labels}
    return np.array(depth), np.array(htot), np.array(ctot), fr


def read_out():
    """Q(H) [s^-1], Cloudy's volume Te(NeNp) [K], abundances n_X/n_H."""
    with open(BASE + '.out') as f:
        lines = f.readlines()
    QH, te, abund = None, None, {}
    for i, line in enumerate(lines):
        if 'Q(1.0-1.8):' in line and QH is None:
            p = line.replace(':', ' ').split()
            QH = sum(10.0 ** float(p[p.index(k) + 1]) for k in
                     ('Q(1.0-1.8)', 'Q(1.8-4.0)', 'Q(4.0-20)', 'Q(20--)'))
        if line.strip().startswith('Volume:'):
            te = float(line.split()[3])          # Te(NeNp) column
        if 'Gas Phase Chemical Composition' in line and not abund:
            p = ' '.join(lines[i + 1:i + 3]).replace(':', ' ').split()
            abund = {p[k].lower(): 10.0 ** float(p[k + 1]) for k in range(0, len(p) - 1, 2)}
    return QH, te, abund


def cloudy_state():
    names, ovr = read_table(BASE + '.ovr')
    c = {n: i for i, n in enumerate(names)}
    _, rad = read_table(BASE + '.rad')
    s = {'depth': ovr[:, c['depth']], 'T': ovr[:, c['Te']],
         'nH': ovr[:, c['hden']], 'ne': ovr[:, c['eden']],
         'xH0': ovr[:, c['HI']], 'xH1': ovr[:, c['HII']],
         'xHe0': ovr[:, c['HeI']], 'xHe1': ovr[:, c['HeII']],
         'xHe2': ovr[:, c['HeIII']]}
    if not np.allclose(rad[:, 2], s['depth'], rtol=1e-4):
        raise RuntimeError('.rad and .ovr zones do not line up')
    #--- .rad radius is the zone midpoint (r_0 + depth), dr the thickness
    s['r'], s['dr'] = rad[:, 1], rad[:, 3]
    s['dV'] = 4.0 * np.pi * s['r'] ** 2 * s['dr']
    s['QH'], s['te_nenp_cloudy'], s['abund'] = read_out()
    s['nHI'] = s['nH'] * s['xH0']
    #--- metal stage fractions: .ele_<el> columns X, X+, X+2, ...
    s['frac'] = {}
    for el in EL_SYM:
        names, ele = read_table(BASE + '.ele_' + el)
        if not np.allclose(ele[:, 0], s['depth'], rtol=1e-5):
            raise RuntimeError('.ele_' + el + ' zones do not line up')
        sym = EL_SYM[el]
        stages = [sym, sym + '+'] + [sym + '+%d' % j for j in range(2, 20)]
        for j, nm in enumerate(stages):
            if nm in names:
                s['frac'][(el, j + 1)] = ele[:, names.index(nm)]
    return s


def n_ion(s, ion):
    el, st = ion.split('_')
    return s['nH'] * s['abund'][el] * s['frac'][(el, int(st))]


def cloudy_channels(s):
    depth, htot, ctot, fr = read_cool(BASE + '.cool')
    if not np.allclose(depth, s['depth'], rtol=1e-5):
        raise RuntimeError('.cool zones do not line up')
    L = {k: np.sum(v * ctot * s['dV']) for k, v in fr.items()}
    return L, np.sum(htot * s['dV']), np.sum(ctot * s['dV'])


# ---------------------------------------------------------------------------
# MoCHII terms
# ---------------------------------------------------------------------------
def run_driver(T, ne, nHI):
    """Coefficients of every ion from fixed_state_nlevel.x -> {ion: {col: array}}."""
    if not os.path.exists(DRIVER):
        raise RuntimeError(DRIVER + ' not built; use run_fixed_state.sh')
    with tempfile.TemporaryDirectory() as tmp:
        fin, fout = os.path.join(tmp, 'state.txt'), os.path.join(tmp, 'coeff.txt')
        with open(fin, 'w') as f:
            f.write('# T[K] ne[cm-3] nHI[cm-3]\n%d\n' % len(T))
            for a, b, c in zip(T, ne, nHI):
                f.write('%.16e %.16e %.16e\n' % (a, b, c))
        subprocess.run([DRIVER, fin, fout], cwd=ROOT, check=True)
        with open(fout) as f:
            hdr = f.readline()[1:].split()
        data = np.loadtxt(fout, ndmin=2)
    out = {}
    for j, name in enumerate(hdr[1:]):
        ion, col = name.split(':')
        out.setdefault(ion, {})[col] = data[:, j + 1]
    return out


def split_coef(coef, nz, k):
    """Zones k*nz..(k+1)*nz-1 of a driver result on stacked states."""
    return {ion: {c: v[k * nz:(k + 1) * nz] for c, v in d.items()}
            for ion, d in coef.items()}


#--- recombination cooling before 2026-09-24 (cooling_mod at git HEAD then)
def hg97_A(T, TTR, scale=1.0):
    lam = 2.0 * TTR / T
    return scale * 1.778e-29 * T * lam ** 1.965 / (1.0 + (lam / 0.541) ** 0.502) ** 2.697


def hg97_B(T, TTR, scale=1.0):
    lam = 2.0 * TTR / T
    return scale * 3.435e-30 * T * lam ** 1.970 / (1.0 + (lam / 2.250) ** 0.376) ** 3.720


def alpha1_HeII_mao2016(T):
    """alpha_1(He I) of 2026-07-25: Mao & Kaastra (2016) coefficients
    (src/recomb_mod.f90 at 297fe79)."""
    return mr.rr_mao(T, 3.7790e-2, 1.1400, 3.9760e-3, 1.2030e2, 9.9590e-1,
                     3.6800, 4.1030e-1)


def recombination_betas(T, model):
    """(b_HII, b_HeII, b_HeIII) [erg cm^3 s^-1] of a REC_MODELS entry."""
    tag, case = model.split('_')
    if tag == 'current':
        A = case == 'A'
        return ((mr.betaA_HII(T) if A else mr.betaB_HII(T)),
                mr.beta_HeII(T, A), mr.beta_HeIII(T, A))
    kT = mr.KB_CGS * T
    if case == 'A':
        return (hg97_A(T, mr.T_TR_HI), kT * mr.alphaA_HeII(T),
                hg97_A(T, mr.T_TR_HEII, 2.0))
    aB = mr.alphaB_HeII(T) if tag == 'previous' else mr.alphaA_HeII(T) - alpha1_HeII_mao2016(T)
    return hg97_B(T, mr.T_TR_HI), kT * aB, hg97_B(T, mr.T_TR_HEII, 2.0)


def two_level_hcoll(T, nHI, ion):
    """H-impact cooling [erg s^-1 per ion] of species_mod::metal_cooling_H
    before 2026-09-24: two-level [O I] 63 um and [C II] 158 um (the same
    expressions tests/hcoll/check_hcoll.f90 prints for the record)."""
    kB = mr.KB_CGS
    if ion == 'o_1':
        qul = 9.2e-11 * (T / 100.0) ** 0.67
        qlu = 0.6 * qul * np.exp(-227.7 / T)
        return nHI * qlu * 227.7 * kB / (1.0 + nHI * (qul + qlu) / 8.865e-5)
    qul = 7.6e-10 * (T / 100.0) ** 0.14
    qlu = 2.0 * qul * np.exp(-91.25 / T)
    return nHI * qlu * 91.25 * kB / (1.0 + nHI * (qul + qlu) / 2.321e-6)


def mochii_terms(s, T, coef):
    """Volume cooling [erg cm^-3 s^-1] of every MoCHII channel at the zone
    temperatures T on the fixed n_e and fractions of state s."""
    ne, nH = s['ne'], s['nH']
    aHe = s['abund']['he']
    nH0, nH1 = nH * s['xH0'], nH * s['xH1']
    nHe0, nHe1, nHe2 = (nH * aHe * s['xHe0'], nH * aHe * s['xHe1'],
                        nH * aHe * s['xHe2'])
    t = {}
    for ion in METAL_IONS:
        ni, c = n_ion(s, ion), coef[ion]
        t[ion] = {'fit': ne * ni * c['lam_fit'],
                  'localne': ne * ni * c['lam_localne'],
                  'nlev_e': ni * c['P_e'], 'nlev_eH': ni * c['P_eH'],
                  'hjoint': ni * (c['P_eH'] - c['P_e']),
                  'h2lv': (ni * two_level_hcoll(T, s['nHI'], ion)
                           if ion in ('o_1', 'c_2') else np.zeros_like(T))}
    c = coef['h_1']
    t['h_1'] = {'fit': ne * nH0 * c['lam_fit'], 'localne': ne * nH0 * c['lam_localne']}
    for m in REC_MODELS:
        bH, bHe1, bHe2 = recombination_betas(T, m)
        t['rec_' + m] = {'H': ne * nH1 * bH, 'He1': ne * nHe1 * bHe1,
                         'He2': ne * nHe2 * bHe2}
    t['ff_HHe'] = FF_PREF * np.sqrt(T) * ne * ((nH1 + nHe1) * mr.gbar_ff(T, 1)
                                               + 4.0 * nHe2 * mr.gbar_ff(T, 2))
    ffm = np.zeros_like(T)
    for (el, st), x in s['frac'].items():
        Z = st - 1
        if Z >= 1:
            ffm += nH * s['abund'][el] * x * Z * Z * mr.gbar_ff(T, Z)
    t['ff_metal'] = FF_PREF * np.sqrt(T) * ne * ffm
    t['ci'] = ne * (nH0 * mr.ci_HI(T) * ETH_HI + nHe0 * mr.ci_HeI(T) * ETH_HEI
                    + nHe1 * mr.ci_HeII(T) * ETH_HEII) * mr.EV2ERG
    return t


def vint(s, x):
    return float(np.sum(x * s['dV']))


def ion_line(t, ion, model):
    mcol, _, hadd = LINE_MODELS[model]
    x = t[ion][mcol]
    return x + t[ion][hadd] if hadd else x


def line_total(t, s, model):
    hcol = LINE_MODELS[model][1]
    return sum(vint(s, ion_line(t, i, model)) for i in METAL_IONS) + vint(s, t['h_1'][hcol])


def total_cooling(t, s, model, rec):
    r = t['rec_' + rec]
    return (line_total(t, s, model) + vint(s, r['H'] + r['He1'] + r['He2'])
            + vint(s, t['ff_HHe'] + t['ff_metal'] + t['ci']))


# ---------------------------------------------------------------------------
def report_fixed_state(s, t, L, Hint, Cint):
    T0 = s['T']
    w = s['ne'] * s['nH'] * s['xH1'] * s['dV']
    s['te_nenp'] = np.sum(T0 * w) / np.sum(w)
    print('Cloudy c25.00 HII40 fixed state: %d zones, r = %.4e .. %.4e cm'
          % (len(T0), s['r'][0] - s['dr'][0] / 2, s['r'][-1] + s['dr'][-1] / 2))
    print('  integrated heating %.5e erg/s, cooling %.5e, C/H = %.5f'
          % (Hint, Cint, Cint / Hint))
    print('  Q(H) = %.5e s^-1, heating per ionizing photon %.4f eV'
          % (s['QH'], Hint / s['QH'] / mr.EV2ERG))
    print('  Te(NeNp) from the zones %.1f K (Cloudy .out, volume: %.0f K)'
          % (s['te_nenp'], s['te_nenp_cloudy']))

    rows = [('O  3 0.0', 'o_3'), ('O  2 0.0', 'o_2'), ('S  3 0.0', 's_3'),
            ('N  2 0.0', 'n_2'), ('S  2 0.0', 's_2'), ('Ne 3 0.0', 'ne_3'),
            ('N  3 0.0', 'n_3'), ('C  2 0.0', 'c_2'), ('Ne 2 0.0', 'ne_2'),
            ('S  4 0.0', 's_4'), ('C  3 0.0', 'c_3'), ('O  1 0.0', 'o_1'),
            ('N  1 0.0', 'n_1'), ('C  1 0.0', 'c_1')]
    cols = ['zero_density', 'zero_density_h2lv', 'local_ne', 'local_ne_h2lv',
            'local_ne_hjoint', 'nlevel_e', 'nlevel_eH']
    print('\nLine channels, MoCHII / Cloudy; Cloudy in erg/s and as % of its '
          'total cooling.  Columns:')
    for k, m in enumerate(cols):
        print('  (%d) %s' % (k + 1, LINE_LABELS[m]))
    print('%-10s %10s %6s |' % ('channel', 'Cloudy', '%Ctot')
          + ''.join('%7s' % ('(%d)' % (k + 1)) for k in range(len(cols))))
    cl_line = 0.0
    for lab, ion in rows:
        c = L.get(lab, 0.0)
        cl_line += c
        el, st = ion.split('_')
        print('%-10s %10.4e %6.2f |' % ('[%s %s]' % (EL_SYM[el], ROMAN[int(st)]),
                                        c, 100 * c / Cint)
              + ''.join('%7.3f' % (vint(s, ion_line(t, ion, m)) / c) for m in cols))
    c = L.get('ISclinH H  0.0', 0.0)
    cl_line += c
    print('%-10s %10.4e %6.2f |' % ('H I coll', c, 100 * c / Cint)
          + ''.join('%7.3f' % (vint(s, t['h_1'][LINE_MODELS[m][1]]) / c) for m in cols)
          + '   (H I: fit in 1, 2, 6, 7; table in 3-5)')
    c = L.get('ISclinHeHe 1.0', 0.0)
    cl_line += c
    print('%-10s %10.4e %6.2f |  MoCHII carries no He I collisional cooling'
          % ('He I coll', c, 100 * c / Cint))

    print('\nLine totals, MoCHII / Cloudy (Cloudy lines %.4e erg/s, %.2f%% of '
          'its cooling)' % (cl_line, 100 * cl_line / Cint))
    for m in cols:
        print('  %-44s %.4f' % (LINE_LABELS[m], line_total(t, s, m) / cl_line))
    print('  zero-density excess over the n-level solve (e only), percentage '
          'points of the Cloudy line total:')
    exc = sorted(((vint(s, t[i]['fit'] - t[i]['nlev_e']) / cl_line * 100, i)
                  for i in METAL_IONS), reverse=True)
    print('    ' + ', '.join('%s %+.2f' % (i, e) for e, i in exc[:6]))

    print('\nContinua')
    ffc = L.get('FF c 0.0', 0.0)
    print('  free-free: Cloudy %.4e; MoCHII/Cloudy H/He only %.4f, with metals %.4f'
          % (ffc, vint(s, t['ff_HHe']) / ffc, vint(s, t['ff_HHe'] + t['ff_metal']) / ffc))
    rH, rHe1, rHe2 = (L.get('ISrcolH H  0.0', 0.0), L.get('ISrcolHeHe 1.0', 0.0),
                      L.get('ISrcolH He 1.0', 0.0))
    print('  recombination, Cloudy (net of excited-level photoheating): onto H+ %.4e,'
          ' He+ %.4e, He++ %.4e;' % (rH, rHe1, rHe2))
    print('    H+ + He+ %.4e, all three %.4e erg/s' % (rH + rHe1, rH + rHe1 + rHe2))
    print('  %-12s %10s %10s %10s %10s   %s' % ('MoCHII', 'onto H+', 'He+', 'He++',
                                               'sum', 'ratio to Cloudy per channel'))
    for m in REC_MODELS:
        r = t['rec_' + m]
        a, b, c = vint(s, r['H']), vint(s, r['He1']), vint(s, r['He2'])
        print('  %-12s %10.4e %10.4e %10.4e %10.4e   %.3f %.3f %.3f'
              % (m, a, b, c, a + b + c, a / rH, b / rHe1, c / rHe2))
    print('  collisional ionization: Cloudy IScionH H %.4e, MoCHII H/He %.4e'
          % (L.get('IScionH H  0.0', 0.0), vint(s, t['ci'])))
    print('  not evaluated: Cloudy hvFB (metal free-bound) %.4e, Cl lines %.4e erg/s'
          % (L.get('hvFB 0.0', 0.0), sum(v for k, v in L.items() if k.startswith('Cl'))))

    print('\nMoCHII totals on the fixed state / Cloudy total cooling (%.4e erg/s)' % Cint)
    print('  %-44s' % '' + ''.join('%12s' % m for m in REC_MODELS))
    for m in ('zero_density', 'zero_density_h2lv', 'local_ne_h2lv', 'local_ne_hjoint',
              'nlevel_e', 'nlevel_eH'):
        print('  %-44s' % LINE_LABELS[m]
              + ''.join('%12.4f' % (total_cooling(t, s, m, r) / Cint) for r in REC_MODELS))


def report_scaling(s, Hint):
    T0, nz = s['T'], len(s['T'])
    fgrid = np.linspace(0.85, 1.15, 61)
    coef = run_driver(np.concatenate([T0 * f for f in fgrid]),
                      np.tile(s['ne'], fgrid.size), np.tile(s['nHI'], fgrid.size))
    terms = [mochii_terms(s, T0 * f, split_coef(coef, nz, k)) for k, f in enumerate(fgrid)]
    models = ['zero_density', 'zero_density_h2lv', 'local_ne_h2lv', 'local_ne_hjoint']
    roots = {}
    for m in models:
        for r in REC_MODELS:
            y = np.log([total_cooling(tk, s, m, r) / Hint for tk in terms])
            sp = CubicSpline(np.log(fgrid), y)
            lf = brentq(sp, np.log(fgrid[0]), np.log(fgrid[-1]))
            roots[(m, r)] = (np.exp(lf), float(sp(lf, 1)))
    #--- exact check of every root in one driver call
    keys = list(roots)
    fr = np.array([roots[k][0] for k in keys])
    coef = run_driver(np.concatenate([T0 * f for f in fr]),
                      np.tile(s['ne'], fr.size), np.tile(s['nHI'], fr.size))
    print('\nTemperature scaling: f with MoCHII total(f T) = Cloudy heating (%.4e '
          'erg/s).\n  Te1 = f x %.1f K (Te(NeNp) integrated over the zones), '
          'Te2 = f x %.0f K (Cloudy\'s printed\n  Te(NeNp), the value the paper '
          'scales); residual = exact total / heating - 1 at the root'
          % (Hint, s['te_nenp'], s['te_nenp_cloudy']))
    for k, key in enumerate(keys):
        tk = mochii_terms(s, T0 * fr[k], split_coef(coef, nz, k))
        res = total_cooling(tk, s, *key) / Hint - 1.0
        print('  %-44s %-11s f = %.5f  Te1 = %6.1f  Te2 = %6.1f K  dlnL/dlnT = %.3f'
              '  residual %.0e' % (LINE_LABELS[key[0]], key[1], fr[k], fr[k] * s['te_nenp'],
                                   fr[k] * s['te_nenp_cloudy'], roots[key][1], res))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--no-scaling', action='store_true',
                    help='skip the temperature-scaling estimate')
    args = ap.parse_args()
    s = cloudy_state()
    L, Hint, Cint = cloudy_channels(s)
    t = mochii_terms(s, s['T'], run_driver(s['T'], s['ne'], s['nHI']))
    report_fixed_state(s, t, L, Hint, Cint)
    if not args.no_scaling:
        report_scaling(s, Hint)


if __name__ == '__main__':
    main()
