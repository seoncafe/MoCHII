#!/usr/bin/env python3
"""PN150 L(Hbeta): MoCHII (pn150_p27_nodrift) against Cloudy c25.00 --
emission measure, H recombinations on one alpha_B(T), and Hbeta per
recombination, so the 10.7% is split into 'fewer recombinations' and
'less Hbeta per recombination'."""
import sys, numpy as np, h5py
sys.path.insert(0, '../../../../tools/python')
import mochii_rates as mr
PC = 3.0856775814913673e18
HB = 2.587e35            # filled below from the lines files
def lhb(path):
    for r in open(path):
        if r.startswith('# L(Hbeta)'): return float(r.split('=')[1].split()[0])
def g(h, k): return h[k]['data'][()] if isinstance(h[k], h5py.Group) else h[k][()]
# ---- MoCHII
with h5py.File('../pn150_p27_nodrift_rates.h5', 'r') as h:
    d = h['Gamma_HI'].attrs['DIST_CM'][0]
    V = (g(h, 'LeafSize')*d)**3
    ne = g(h, 'n_e'); xp = g(h, 'x_HII'); T = g(h, 'T_e')
    xhe2 = g(h, 'x_HeII'); xhe3 = g(h, 'x_HeIII')
gas = ne > 0
nH = 3000.0
aB = np.array([mr.alphaB_HII(t) for t in T[gas]])
EM = np.sum(ne[gas]*nH*xp[gas]*V[gas])
RB = np.sum(ne[gas]*nH*xp[gas]*aB*V[gas])
Lm = lhb('../pn150_p27_nodrift_lines.txt')
Vion = np.sum(xp[gas]*V[gas])
# ---- Cloudy
o = np.loadtxt('../../../../tests/cloudy_c25/pn150_c25.ovr')
dep, Tc, hden, eden, HII = o[:,0], o[:,1], o[:,3], o[:,4], o[:,7]
r = 1e17 + dep
edges = np.concatenate([[1e17], 0.5*(r[1:]+r[:-1]), [r[-1] + 0.5*(r[-1]-r[-2])]])
dV = 4/3*np.pi*(edges[1:]**3 - edges[:-1]**3)
aBc = np.array([mr.alphaB_HII(t) for t in Tc])
EMc = np.sum(eden*hden*HII*dV); RBc = np.sum(eden*hden*HII*aBc*dV)
Lc = 2.8966e35
Vionc = np.sum(HII*dV)
hnu = 6.62607015e-27*2.99792458e10/4862.68e-8
print(f"{'':34s}{'MoCHII':>12s}{'Cloudy':>12s}{'ratio':>8s}")
for lab, a, b in [('ionized volume sum x_HII V [cm^3]', Vion, Vionc),
                  ('EM = sum n_e n_p V [cm^-3]', EM, EMc),
                  ('R_B = sum n_e n_p alpha_B V [s^-1]', RB, RBc),
                  ('L(Hbeta) [erg/s]', Lm, Lc),
                  ('Hbeta photons / R_B', Lm/hnu/RB, Lc/hnu/RBc)]:
    print(f"{lab:34s}{a:12.4e}{b:12.4e}{a/b:8.4f}")
print('Q(H) (Cloudy, sum of its four bands) = 5.416e47 s^-1')
print(f"R_B/Q(H): MoCHII {RB/5.416e47:.4f}  Cloudy {RBc/5.416e47:.4f}")
# ---- helium recombinations (case B and the ground captures) on the same alpha
yHe = 0.1
aBHe2 = np.array([mr.alphaB_HeII(t) for t in T[gas]]); aAHe2 = np.array([mr.alphaA_HeII(t) for t in T[gas]])
aBHe3 = np.array([mr.alphaB_HeIII(t) for t in T[gas]]); aAHe3 = np.array([mr.alphaA_HeIII(t) for t in T[gas]])
w = ne[gas]*nH*yHe*V[gas]
m = dict(He2B=np.sum(w*xhe2[gas]*aBHe2), He2A=np.sum(w*xhe2[gas]*aAHe2),
         He3B=np.sum(w*xhe3[gas]*aBHe3), He3A=np.sum(w*xhe3[gas]*aAHe3))
HeII_c, HeIII_c = o[:,9], o[:,10]
wc = eden*hden*yHe*dV
c = dict(He2B=np.sum(wc*HeII_c*np.array([mr.alphaB_HeII(t) for t in Tc])),
         He2A=np.sum(wc*HeII_c*np.array([mr.alphaA_HeII(t) for t in Tc])),
         He3B=np.sum(wc*HeIII_c*np.array([mr.alphaB_HeIII(t) for t in Tc])),
         He3A=np.sum(wc*HeIII_c*np.array([mr.alphaA_HeIII(t) for t in Tc])))
for k in m:
    print(f"R({k}) [s^-1]{'':22s}{m[k]:12.4e}{c[k]:12.4e}{m[k]/c[k]:8.4f}")
aAH = np.array([mr.alphaA_HII(t) for t in T[gas]])
RA = np.sum(ne[gas]*nH*xp[gas]*aAH*V[gas])
aAHc = np.array([mr.alphaA_HII(t) for t in Tc]); RAc = np.sum(eden*hden*HII*aAHc*dV)
print(f"R_A(H) [s^-1]{'':21s}{RA:12.4e}{RAc:12.4e}{RA/RAc:8.4f}")
