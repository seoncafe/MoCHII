#!/usr/bin/env python3
"""H ionizations per shell on Cloudy's PN150 structure (faint.ovr) by diffuse
photons emitted at energy E, under
 (iso)   isotropic emission, straight-line transport through the sphere
         (MoCHII's transported packets);
 (cloudy) Cloudy c25.00's default treatment (sphere, 'OU2'):
   recombination continua (rt_ots.cpp:281-286, iso_radiative_recomb.cpp:432,
   rt_recom_effic.cpp:109-114): a share 1 - E2(tau_face(E_th)) absorbed in
   the emitting zone, E2(tau_face) into the outward radial beam;
   two-photon (rt_diffuse.cpp:160, rt_ots.cpp:565): 1 - exp(-tau_face(E))
   local, exp(-tau_face(E)) outward;
   escaping lines (lines_service.cpp:571-577): all into the outward beam
   (inward half via BeamInOut), destroyed share local;
   outward beam attenuated by exp(-kappa dr) (rt_continuum.cpp:228-238).
Absorbers H0, He0 (no metals, no dust).  Estimate only."""
import sys, numpy as np
from scipy.special import expn
sys.path.insert(0, '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tools/python')
import mochii_rates as mr
d = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/pn150_hbeta_residual/'
o = np.loadtxt(d + 'faint.ovr'); pc = 3.0857e18
r = 1e17 + o[:, 0]
T, nH, ne, xH0, xHII, xHe0, xHe2 = o[:, 1], o[:, 3], o[:, 4], o[:, 6], o[:, 7], o[:, 8], o[:, 10]
edges = np.concatenate([[1e17], 0.5*(r[1:]+r[:-1]), [r[-1] + 0.5*(r[-1]-r[-2])]])
dr = edges[1:] - edges[:-1]
dV = 4/3*np.pi*(edges[1:]**3 - edges[:-1]**3)
nz = len(r)
RH = ne*nH*xHII*np.array([mr.alphaB_HII(t) for t in T])*dV
S3 = ne*nH*0.1*xHe2*np.array([mr.alphaB_HeIII(t) for t in T])*dV     # He2+ case-B recombinations
S1 = ne*nH*xHII*np.array([mr.alpha1_HII(t) for t in T])*dV           # H ground captures

def kap(E):
    kH = nH*xH0*mr.sigma_HI(E)
    kHe = nH*0.1*xHe0*mr.sigma_HeI(E) if E > 24.587 else 0*nH
    return kH + kHe, kH/np.maximum(kH + kHe, 1e-300)

def tau_face(E):
    k, _ = kap(E); t = np.cumsum(k*dr); return t - 0.5*k*dr     # to zone center

def radial_beam(E):
    k, _ = kap(E); A = np.zeros((nz, nz))
    for i in range(nz):
        f = 1.0
        # rt_continuum.cpp:237-238 adds the zone's escape then attenuates by this zone
        for j in range(i, nz):
            a = f*(1 - np.exp(-k[j]*dr[j])); A[i, j] = a; f -= a
            if f < 1e-12: break
    return A

cache = {}
def isotropic(E, nmu=160, nsub=1, nstep=8000):
    if E in cache: return cache[E]
    k, _ = kap(E); A = np.zeros((nz, nz))
    mu = (np.arange(nmu) + 0.5)/nmu*2 - 1
    s = np.linspace(0, 2*edges[-1], nstep + 1); ds = s[1] - s[0]; sm = s[:-1] + ds/2
    for i in range(nz):
        for rs in np.linspace(edges[i], edges[i+1], nsub + 2)[1:-1]:
            rr = np.sqrt(rs**2 + sm[None, :]**2 + 2*rs*sm[None, :]*mu[:, None])
            idx = np.searchsorted(edges, rr) - 1
            inside = (idx >= 0) & (idx < nz)
            kk = np.where(inside, k[np.clip(idx, 0, nz - 1)], 0.0)
            tau = np.cumsum(kk*ds, axis=1); dab = np.exp(-(tau - kk*ds)) - np.exp(-tau)
            A[i] += np.bincount(idx[inside], weights=dab[inside], minlength=nz)[:nz]/nmu/nsub
    cache[E] = A; return A

def cloudy(E, kind, Eth=None, pdest=0.0, bowen=False):
    A = radial_beam(E)
    if kind == 'rrc':
        esc = expn(2, tau_face(Eth))
    elif kind == '2q':
        esc = np.exp(-tau_face(E))
    else:                      # line: destroyed share local
        esc = np.full(nz, 1.0 - pdest)
    loc = np.diag(1.0 - esc)
    return loc + esc[:, None]*A

rpc = r/pc
bins = [(0.0324, 0.06), (0.06, 0.09), (0.09, 0.136)]
def report(name, S, w, cases):
    """cases: list of (weight per source event, E, cloudy-kernel kwargs)"""
    hi_c = np.zeros(nz); hi_i = np.zeros(nz)
    for wt, E, kw in cases:
        _, fH = kap(E)
        hi_c += wt*(S[:, None]*cloudy(E, **kw)).sum(0)*fH
        hi_i += wt*(S[:, None]*isotropic(E)).sum(0)*fH
    print(f'{name}: H ionizations / R_B(H) total  cloudy {hi_c.sum()/RH.sum():.4f}  iso {hi_i.sum()/RH.sum():.4f}')
    for a, b in bins:
        m = (rpc >= a) & (rpc < b)
        print(f'   {a:.4f}-{b:.4f} pc  cloudy/R_H {hi_c[m].sum()/RH[m].sum():.4f}  iso/R_H {hi_i[m].sum()/RH[m].sum():.4f}'
              f'  diff {(hi_c[m].sum() - hi_i[m].sum())/RH[m].sum():+.4f}')
    return hi_c, hi_i

print('He2+ rec / H rec (case B, whole model):', S3.sum()/RH.sum())
# He II cascade per case-B He2+ recombination (photon counts: estimates, see the note)
heii = report('He II fb (n>=2) 0.35/rec at 13.7 eV', S3, None, [(0.35, 13.7, dict(kind='rrc', Eth=13.604))])
heii2 = report('He II 2q 0.30 decays x (0.5 at 16 + 0.5 at 24 eV) [~1 ionizing photon per decay]', S3, None,
               [(0.15, 16.0, dict(kind='2q')), (0.15, 24.0, dict(kind='2q'))])
heii3 = report('He II Ly-alpha 0.70/rec, Pdest share 0', S3, None, [(0.70, 40.8, dict(kind='line', pdest=0.0))])
heii4 = report('He II Ly-alpha 0.70/rec, Pdest share 0.10', S3, None, [(0.70, 40.8, dict(kind='line', pdest=0.10))])
h1 = report('H I ground continuum (alpha_1) at 13.7 eV', S1, None, [(1.0, 13.7, dict(kind='rrc', Eth=13.598))])

print('\nCombined, (cloudy - iso)/R_B(H) in the bins of the fixed-structure test')
for lab, parts in [('He II, Pdest 0', [heii, heii2, heii3]), ('He II, Pdest 0.10', [heii, heii2, heii4]),
                   ('H I ground', [h1]), ('H I ground + He II (Pdest 0.10)', [h1, heii, heii2, heii4])]:
    c = sum(p[0] for p in parts); i = sum(p[1] for p in parts)
    out = []
    for a, b in [(0.0324, 0.06), (0.06, 0.09), (0.09, 0.136)]:
        m = (rpc >= a) & (rpc < b)
        out.append(f'{a:.3f}-{b:.3f}: {(c[m].sum() - i[m].sum())/RH[m].sum():+.4f}')
    print(f'  {lab:34s} ' + '   '.join(out))
