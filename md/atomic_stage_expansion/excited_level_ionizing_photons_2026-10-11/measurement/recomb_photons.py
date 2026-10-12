"""Ionizing photons (>= 13.598 eV) from metal recombination, per recombination,
from Badnell final-state-resolved partial rates (adf48 RR / adf09 DR, PRTI=1)
and the radiative cascade of the MoCHII CHIANTI atom of the recombined ion.

Per final state (parent p, n, l):
  first photon  h nu1 = E_in + IP - E_f   (E_in: electron energy for RR,
                                          resonance energy for DR)
  cascade       radiative branching through the CHIANTI atom (A values only).
Writes coefficient tables [eV cm^3 s^-1] on a T grid for each target ion.
"""
import sys, re, glob, json, os
import numpy as np
from scipy.special import exp1
sys.path.insert(0, os.path.dirname(__file__))
import adf

SNAP = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase4_20261011/data_snapshot/atomic'
CHI = '/home/kiseon/RT_Codes/CHIANTI/dbase'
BAD = os.environ.get('BADNELL_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'badnell_dl'))
KB_EV = 8.617333262e-5
CM_EV = 1.0 / 8065.543937
RY_EV = 13.605693
THR = 13.598
THR_HE = 24.587
LSYM = 'spdfghiklmnoqrtuvwxyz'
LCAP = {c.upper(): i for i, c in enumerate('SPDFGHIKLMNOQRTUV')}
ZEL = dict(c=6, n=7, o=8, ne=10, s=16)
SEQ = {0: 'ur', 1: 'h', 2: 'he', 3: 'li', 4: 'be', 5: 'b', 6: 'c', 7: 'n', 8: 'o', 9: 'f',
       10: 'ne', 11: 'na', 12: 'mg', 13: 'al', 14: 'si', 15: 'p'}


ORDER = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '4f']
FULL = {'s': 2, 'p': 6, 'd': 10, 'f': 14, 'g': 18}


def parse_occ(tokens):
    occ = {}
    for t in tokens:
        if t.isdigit():
            continue
        m = re.match(r'^(\d+)([a-z])(\d*)$', t)
        if not m:
            return None
        orb = m.group(1) + m.group(2)
        occ[orb] = occ.get(orb, 0) + (int(m.group(3)) if m.group(3) else 1)
    return occ


def conf_key(occ, N):
    """Full configuration with implicit inner shells filled to N electrons,
    closed shells removed: one key for CHIANTI and Badnell spellings."""
    if occ is None:
        return None
    occ = dict(occ)
    need = N - sum(occ.values())
    for orb in ORDER:
        if need <= 0:
            break
        if orb in occ:
            continue
        k = min(FULL[orb[-1]], need)
        occ[orb] = k
        need -= k
    if need != 0:
        return None
    return frozenset((o, c) for o, c in occ.items() if c != FULL.get(o[-1], -1) and c > 0)


def chianti_occ(s):
    s = re.sub(r'\([^)]*\)', ' ', s.lower())
    s = re.sub(r'\b\d+\b', ' ', s)          # stray numeric columns
    toks = re.findall(r'(\d)([sp])([1-6]?)|(\d)([dfgh])(1[0-4]|[1-9])?', s.replace(' ', '.'))
    out = []
    for a, b, c, d, e, f in toks:
        out.append(a + b + c if a else d + e + (f or ''))
    return parse_occ(out)


def badnell_occ(code):
    return parse_occ([t.lower() for t in code.split()])


def add_orbital(pocc, n, l, N):
    if pocc is None:
        return None
    occ = dict(pocc)
    orb = f'{n}{LSYM[l]}'
    occ[orb] = occ.get(orb, 0) + 1
    return conf_key(occ, N)


def read_element(el):
    tr = {}
    cur = None
    for line in open(f'{SNAP}/element_{el}.txt'):
        t = line.split()
        if not t or t[0].startswith('#'):
            continue
        if t[0] == 'TRANSITION':
            cur = int(t[1]); tr[cur] = {}
        elif t[0] == 'ETH':
            tr[cur]['ip'] = float(t[1])
        elif t[0] == 'RR':
            tr[cur]['rr'] = [float(x) for x in t[1:]]
        elif t[0] == 'DR':
            k = int(t[1]); v = [float(x) for x in t[2:]]
            tr[cur]['dr'] = (np.array(v[:k]), np.array(v[k:2 * k]))
    return tr


def alpha_rr(p, T):
    A, B, T0, T1 = p[:4]
    if len(p) > 4 and p[4] != 0:
        B = B + p[4] * np.exp(-p[5] / T)
    return A / (np.sqrt(T / T0) * (1 + np.sqrt(T / T0)) ** (1 - B) * (1 + np.sqrt(T / T1)) ** (1 + B))


def alpha_dr(p, T):
    c, E = p
    return T ** -1.5 * np.sum(c[None, :] * np.exp(-E[None, :] / T[:, None]), axis=1)


class Atom:
    """MoCHII nlevel atom + CHIANTI configuration labels; radiative cascade."""

    def __init__(self, el, stage):
        self.N = ZEL[el] - stage + 1
        f = f'{SNAP}/excitation/chianti/nlevel_{el}_{stage}.txt'
        lines = open(f).read().splitlines()
        ids = None
        for L in lines:
            if L.startswith('# source level IDs:'):
                ids = [int(x) for x in L.split(':')[1].split()]
        i = [k for k, L in enumerate(lines) if L.startswith('NLEV ')][0]
        nlev = int(lines[i].split()[1])
        self.g = np.zeros(nlev); self.E = np.zeros(nlev)
        for k in range(nlev):
            t = lines[i + 1 + k].split()
            self.g[k] = float(t[1]); self.E[k] = float(t[2]) * CM_EV
        j = [k for k, L in enumerate(lines) if L.startswith('NRAD ')][0]
        nrad = int(lines[j].split()[1])
        self.A = np.zeros((nlev, nlev))   # A[u, l]
        for k in range(nrad):
            t = lines[j + 1 + k].split()
            a, b, A = int(t[0]) - 1, int(t[1]) - 1, float(t[2])
            u, l = (a, b) if self.E[a] >= self.E[b] else (b, a)
            if u != l:
                self.A[u, l] += A
        # labels
        lab = {}
        for L in open(f'{CHI}/{el}/{el}_{stage}/{el}_{stage}.elvlc'):
            t = L.split()
            if not t or t[0] == '-1':
                break
            m = re.match(r'^\s*(\d+)\s+(.+?)\s+(\d+)\s+([A-Za-z])\s+([\d.]+)\s', L)
            if m:
                lab[int(m.group(1))] = (conf_key(chianti_occ(m.group(2)), self.N), int(m.group(3)), LCAP.get(m.group(4).upper(), -1), float(m.group(5)))
        self.conf = [lab[s][0] if s in lab else None for s in ids]
        self.term = [lab[s][1:3] if s in lab else None for s in ids]
        self.J = [lab[s][3] if s in lab else None for s in ids]
        self.ground_term = [(self.conf[k] == self.conf[0] and self.term[k] == self.term[0]) for k in range(nlev)]
        self.cascade()
        self.byconf = {}
        for k, c in enumerate(self.conf):
            if c is not None:
                self.byconf.setdefault(c, []).append(k)

    def cascade(self):
        n = len(self.E)
        order = np.argsort(self.E)
        self.C = np.zeros(n); self.Che = np.zeros(n)
        for u in order:
            Atot = self.A[u].sum()
            if Atot <= 0:
                continue
            for l in np.nonzero(self.A[u])[0]:
                b = self.A[u, l] / Atot
                h = self.E[u] - self.E[l]
                self.C[u] += b * ((h if h >= THR else 0.0) + self.C[l])
                self.Che[u] += b * ((h if h >= THR_HE else 0.0) + self.Che[l])

    def has_parent(self, pocc):
        k = frozenset(pocc.items())
        if not hasattr(self, '_hp'):
            self._hp = {}
        if k not in self._hp:
            self._hp[k] = any(add_orbital(pocc, nn, ll, self.N) in self.byconf for nn in range(1, 10) for ll in range(nn))
        return self._hp[k]

    def group(self, key):
        ks = self.byconf.get(key)
        if not ks:
            return None
        w = self.g[ks] / self.g[ks].sum()
        gt = sum(w[i] for i, k in enumerate(ks) if self.ground_term[k])
        return dict(E=np.dot(w, self.E[ks]), C=np.dot(w, self.C[ks]), Che=np.dot(w, self.Che[ks]), gt=gt)


def final_state(atom, pocc, Ep, pocc_ground, n, l, ip, q):
    """(E_f, C, Che, ground-term share, how) for capture into nl on the parent
    with explicit occupancy pocc (energy Ep above the parent ground)."""
    N = atom.N
    g = atom.group(add_orbital(pocc, n, l, N))
    if g is not None:
        return g['E'], g['C'], g['Che'], g['gt'], 'atom'
    Ef = ip + Ep - q * q * RY_EV / n ** 2
    if pocc is not None and atom.has_parent(pocc):
        base, core, core_he = pocc, 0.0, 0.0
    else:
        base = pocc_ground
        core = Ep if Ep >= THR else 0.0
        core_he = Ep if Ep >= THR_HE else 0.0
    for nn in range(min(n - 1 if base is pocc else n, 9), 0, -1):
        for ll in range(min(l, nn - 1), -1, -1):
            g = atom.group(add_orbital(base, nn, ll, N))
            if g is not None:
                dE = Ef - (Ep if base is not pocc else 0.0) - g['E']
                if os.environ.get('DESCENT_NO_DIRECT'):
                    dE = 0.0          # sensitivity: no direct high-n -> low-n jump photons
                C = g['C'] + core + (dE if dE >= THR else 0.0)
                Che = g['Che'] + core_he + (dE if dE >= THR_HE else 0.0)
                return Ef, C, Che, 0.0, 'descend'
    return Ef, (Ef if Ef >= THR else 0.0), (Ef if Ef >= THR_HE else 0.0), 0.0, 'none'


def nl_rates(d, Tgrid):
    """{(p, n, l): rate(T)} from ILREP/INREP on Tgrid (log-log interpolation)."""
    out = {}
    lt = np.log(d['TE'])
    def interp(v):
        return np.exp(np.interp(np.log(Tgrid), lt, np.log(np.maximum(v, 1e-300))))
    for p, blk in d['prtf'].items():
        nl = blk.get('nl', {}); nn = blk.get('n', {})
        ns_nl = sorted({k[0] for k in nl})
        reps = sorted(nn)
        if not reps:
            continue
        lr = np.log(reps)
        ints = np.arange(reps[0], reps[-1] + 1)
        cell = np.argmin(abs(ints[:, None] - np.array(reps)[None, :]), axis=1)
        logv = np.array([np.log(np.maximum(interp(nn[k]), 1e-300)) for k in reps])  # (nrep, nT)
        fac = {}
        for j, k in enumerate(reps):
            mem = ints[cell == j]
            li = np.array([np.interp(np.log(mem), lr, logv[:, t]) for t in range(logv.shape[1])])  # (nT, nmem)
            fac[k] = np.exp(li).sum(axis=1) / np.exp(logv[j])
        for n, v in nn.items():
            vn = interp(v) * fac[n]
            ls = {k[1]: interp(vv) for k, vv in nl.items() if k[0] == n}
            if ls:
                s = sum(ls.values())
                for l, vl in ls.items():
                    out[(p, n, l)] = out.get((p, n, l), 0) + vl * np.where(s > 0, vn / np.maximum(s, 1e-300), 0)
            else:
                ref = max([m for m in ns_nl if m < n], default=None)
                if ref is None:
                    out[(p, n, 0)] = out.get((p, n, 0), 0) + vn
                    continue
                ls = {k[1]: interp(vv) for k, vv in nl.items() if k[0] == ref}
                s = sum(ls.values())
                for l, vl in ls.items():
                    out[(p, n, l)] = out.get((p, n, l), 0) + vn * vl / np.maximum(s, 1e-300)
    return out


def ion_coefficients(el, q, Tgrid, files_rr, files_dr):
    Z = ZEL[el]
    tr = read_element(el)[q]          # recombination X^{q+} -> X^{(q-1)+}, recombined stage q
    ip = tr['ip']
    atom = Atom(el, q)
    kT = KB_EV * Tgrid
    res = {k: np.zeros_like(Tgrid) for k in
           ['rr_gt', 'rr_fb_exc', 'rr_fb_exc_he', 'rr_cas', 'rr_cas_he', 'rr_tot_energy',
            'dr_stab', 'dr_stab_he', 'dr_cas', 'dr_cas_he', 'dr_tot_energy',
            'alpha_rr_partial', 'alpha_dr_partial', 'rr_unmatched', 'dr_unmatched', 'rr_gt_rate']}
    cdr, Edr = tr.get('dr', (np.array([0.0]), np.array([1.0])))
    if len(cdr) == 0:
        cdr, Edr = np.array([0.0]), np.array([1.0])
        tr['dr'] = (cdr, Edr)
    lw = np.log(np.maximum(cdr, 1e-300))[None, :] - Edr[None, :] / Tgrid[:, None]
    lw = np.where(cdr[None, :] > 0, lw, -np.inf)
    mx = lw.max(axis=1, keepdims=True)
    mx = np.where(np.isfinite(mx), mx, 0.0)
    w = np.exp(lw - mx)
    w = w / np.maximum(w.sum(axis=1, keepdims=True), 1e-300)
    Eres = Edr * KB_EV
    cache = {}
    for kind, files in (('rr', files_rr), ('dr', files_dr)):
        for f in files:
            d = adf.read(f)
            par = d['parents']
            pg = badnell_occ(par[1]['code']) if 1 in par else None
            rates = nl_rates(d, Tgrid)
            channels = []
            # IC-resolved final levels (adf48/adf09 INDX block), matched to CHIANTI levels
            res_nl = {}
            lt = np.log(d['TE'])
            for k, v in d['indx'].items():
                lev = d['levels'].get(k)
                if lev is None:
                    continue
                toks = lev['code'].lower().split()
                m = re.match(r'^(\d+)([a-z])', toks[-1])
                if not m or m.group(2) not in LSYM:
                    continue
                nl = (int(m.group(1)), LSYM.index(m.group(2)))
                key = conf_key(badnell_occ(lev['code']), atom.N)
                ks = atom.byconf.get(key)
                if not ks:
                    continue
                a = np.exp(np.interp(np.log(Tgrid), lt, np.log(np.maximum(v, 1e-300))))
                exact = [i for i in ks if atom.term[i] == (lev['S'], lev['L']) and atom.J[i] == lev['J']]
                if exact:
                    i = exact[0]
                    channels.append((a, atom.E[i], atom.C[i], atom.Che[i], float(atom.ground_term[i]), 'level'))
                else:
                    g = atom.group(key)
                    channels.append((a, g['E'], g['C'], g['Che'], g['gt'], 'atom'))
                res_nl[nl] = res_nl.get(nl, 0) + a
            for nl, ares in res_nl.items():
                tot = sum(v for (p, n, l), v in rates.items() if (n, l) == nl)
                if np.isscalar(tot) and tot == 0:
                    continue
                f_left = np.clip(1.0 - ares / np.maximum(tot, 1e-300), 0.0, 1.0)
                for key3 in [k3 for k3 in rates if (k3[1], k3[2]) == nl]:
                    rates[key3] = rates[key3] * f_left
            for (p, n, l), a in rates.items():
                pinfo = par.get(p)
                pconf = badnell_occ(pinfo['code']) if pinfo else None
                Ep = pinfo['E'] * CM_EV if pinfo else 0.0
                ck = (f, p, n, l)
                if ck not in cache:
                    cache[ck] = final_state(atom, pconf, Ep, pg, n, l, ip, q)
                Ef, C, Che, gt, how = cache[ck]
                channels.append((a, Ef, C, Che, gt, how))
            for a, Ef, C, Che, gt, how in channels:
                B = ip - Ef
                res[f'alpha_{kind}_partial'] += a
                if how not in ('atom', 'level'):
                    res[f'{kind}_unmatched'] += a
                if kind == 'rr':
                    if B > 0:
                        b = B / kT
                        with np.errstate(over='ignore', invalid='ignore'):
                            eb = np.where(b < 50, np.exp(np.minimum(b, 50)) * exp1(b), (1 - 1 / b + 2 / b**2 - 6 / b**3) / b)
                        mean_hnu = kT / eb
                    else:
                        mean_hnu = kT
                    ion = mean_hnu * np.exp(-np.maximum(0.0, THR - B) / kT)
                    ion_he = mean_hnu * np.exp(-np.maximum(0.0, THR_HE - B) / kT)
                    res['rr_gt'] += a * gt * mean_hnu
                    res['rr_gt_rate'] += a * gt
                    res['rr_fb_exc'] += a * (1 - gt) * ion
                    res['rr_fb_exc_he'] += a * (1 - gt) * ion_he
                    res['rr_cas'] += a * C
                    res['rr_cas_he'] += a * Che
                    res['rr_tot_energy'] += a * (mean_hnu + Ef)
                else:
                    hnu = Eres[None, :] + B
                    ion = np.sum(w * np.where(hnu >= THR, hnu, 0.0), axis=1)
                    ion_he = np.sum(w * np.where(hnu >= THR_HE, hnu, 0.0), axis=1)
                    res['dr_stab'] += a * ion
                    res['dr_stab_he'] += a * ion_he
                    res['dr_cas'] += a * C
                    res['dr_cas_he'] += a * Che
                    res['dr_tot_energy'] += a * (np.sum(w * np.maximum(hnu, 0), axis=1) + Ef)
    res['alpha_rr'] = alpha_rr(tr['rr'], Tgrid)
    res['alpha_dr'] = alpha_dr(tr['dr'], Tgrid)
    res['ip'] = ip
    res['Eres_mean'] = (w * Eres[None, :]).sum(axis=1)
    return res


def files_for(el, q):
    Z = ZEL[el]; N = Z - q; s = SEQ[N]
    rr = sorted(glob.glob(f'{BAD}/rr/*#{s}_{el}{q}ic.dat.gz'))
    dr = sorted(glob.glob(f'{BAD}/dr/*#{s}_{el}{q}ic*.dat.gz'))
    # one file per core-transition code; skip the 222/223 variants; prefer za00 over nrb00
    pick = {}
    for f in dr:
        code = re.search(r'ic[m]?(\d+)\.dat', f).group(1)
        if code in ('222', '223'):
            continue
        if code not in pick or os.path.basename(f).startswith('za00'):
            pick[code] = f
    return rr[:1], sorted(pick.values())


if __name__ == '__main__':
    Tgrid = 10 ** np.arange(3.4, 4.801, 0.025)
    out = {}
    for el in ['c', 'n', 'o', 'ne', 's']:
        for q in range(1, 7):
            rr, dr = files_for(el, q)
            r = ion_coefficients(el, q, Tgrid, rr, dr)
            out[f'{el}{q}'] = {k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in r.items()}
            out[f'{el}{q}']['files'] = [os.path.basename(x) for x in rr + dr]
            i = np.argmin(abs(Tgrid - 1e4))
            print(f'{el}{q:d} IP={r["ip"]:7.2f}  aRR part/fit {r["alpha_rr_partial"][i]/r["alpha_rr"][i]:6.3f}'
                  f'  aDR part/fit {r["alpha_dr_partial"][i]/max(r["alpha_dr"][i],1e-300):8.3f}'
                  f'  per recomb [eV] @1e4K: gt {r["rr_gt"][i]/r["alpha_rr_partial"][i]:6.2f}'
                  f' fbexc {r["rr_fb_exc"][i]/r["alpha_rr_partial"][i]:6.2f} casRR {r["rr_cas"][i]/r["alpha_rr_partial"][i]:6.2f}'
                  f' | stab {r["dr_stab"][i]/max(r["alpha_dr_partial"][i],1e-300):6.2f} casDR {r["dr_cas"][i]/max(r["alpha_dr_partial"][i],1e-300):6.2f}'
                  f' | unmatched RR {r["rr_unmatched"][i]/r["alpha_rr_partial"][i]:5.2f} DR {r["dr_unmatched"][i]/max(r["alpha_dr_partial"][i],1e-300):5.2f}'
                  f' | Etot RR/(IP) {r["rr_tot_energy"][i]/r["alpha_rr_partial"][i]/r["ip"]:5.3f}', flush=True)
    out['T'] = Tgrid.tolist()
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     os.environ.get('RECOMB_JSON', 'recomb_coefficients.json')), 'w'))
