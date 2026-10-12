"""Thermally averaged free-free Gaunt factors for tools/python/mochii_rates.py.

gbar(T, Z) = Int_0^inf g_ff(u) exp(-u) du, u = h nu / kT, with the Gaunt rule
of src/free_free_mod.f90 (free_free_gaunt): the Hummer (1988) fit inside its
stated range (1e-4 <= u <= 10^1.5, 1e-3 <= gamma^2 = Z^2 Ry/kT <= 1e3,
T >= 10^2.5 K), the van Hoof et al. (2014) table (data/gauntff_vh14.dat,
bilinear in log gamma^2 and log u) everywhere else. Both are implemented here
independently of the Fortran. The integral is 40-point Gauss-Legendre on
every interval between the kinks of the integrand (the van Hoof u nodes and
the Hummer range edges) over 1e-14 <= u <= 700; halving the order and the
range changes it by at most 1.3e-7 (measured 2026-10-11).

Writes gbar_ff_table.txt: log10 T = 2.00 (0.01) 6.00, Z = 1..6, the grid of
cooling_mod::gff_z.  Run: python3 make_gbar_ff_table.py
"""
import numpy as np, sys
# Hummer (1988, ApJ 327, 477) Chebyshev coefficients d(8,11), as printed in
# his routine (ported verbatim to src/gaunt.f90 from MOCASSIN).
D=np.array([8.986940175e+00,-4.009515855e+00,8.808871266e-01,2.640245111e-02,-4.580645915e-02,-3.568055702e-03,2.827798067e-03,3.365860195e-04,-8.006936989e-01,9.466021705e-01,9.043402532e-02,-9.608451450e-02,-1.885629865e-02,1.050313890e-02,2.800889961e-03,-1.078209202e-03,-3.781305103e-01,1.102726332e-01,-1.543619180e-02,8.310561114e-03,2.179620525e-02,4.259726289e-03,-4.181588794e-03,-1.770208330e-03,1.877213132e-02,-1.004885705e-01,-5.483366378e-02,-4.520154409e-03,8.366530426e-03,3.700273930e-03,6.889320423e-04,9.460313195e-05,7.300158392e-02,3.576785497e-03,-4.545307025e-03,-1.017965604e-02,-9.530211924e-03,-3.450186162e-03,1.040482914e-03,1.407073544e-03,-1.744671550e-03,2.864013856e-02,1.903394837e-02,7.091074494e-03,-9.668371391e-04,-2.999107465e-03,-1.820642230e-03,-3.874082085e-04,-1.707268366e-02,-4.694254776e-03,1.311691517e-03,5.316703136e-03,5.178193095e-03,2.451228935e-03,-2.277321615e-05,-8.182359057e-04,2.567331664e-04,-9.155339970e-03,-6.997479192e-03,-3.571518641e-03,-2.096101038e-04,1.553822487e-03,1.509584686e-03,6.212627837e-04,4.098322531e-03,1.635218463e-03,-5.918883504e-04,-2.333091048e-03,-2.484138313e-03,-1.359996060e-03,-5.371426147e-05,5.553549563e-04,3.837562402e-05,2.938325230e-03,2.393747064e-03,1.328839809e-03,9.135013312e-05,-7.137252303e-04,-7.656848158e-04,-3.504683798e-04,-8.491991820e-04,-3.615327726e-04,3.148015257e-04,8.909207650e-04,9.869737522e-04,6.134671184e-04,1.068883394e-04,-2.046080100e-04]).reshape(11,8)  # D[i,j] = d(j,i)
def clenshaw(coef,x):
    # Hummer's form: b_n recursion, result 0.25*(b1-b3) for txg; for txu b1-b3
    b=np.zeros(len(coef)+2)
    for k in range(len(coef)-1,-1,-1): b[k]=x*b[k+1]-b[k+2]+coef[k]
    return b
def hummer(Z,T,u):
    # inside its stated range only (Chebyshev branch); Hummer's u is from xlf=log10(u/(Ry/kT))
    xlrkt=5.1983649-np.log10(T)
    txg=0.66666667*(2*np.log10(Z)+xlrkt)
    c=np.empty(8)
    for j in range(8):
        b=clenshaw(D[:,j],txg); c[j]=0.25*(b[0]-b[2])
    con=0.72727273*xlrkt+0.90909091
    xlf=np.log10(u)-xlrkt  # log10(h nu/Ry) = log10(u) + log10(kT/Ry)
    txu=0.72727273*xlf+con
    assert np.all(np.abs(txu)<=2.0+1e-9), (Z,T)
    out=np.empty_like(u)
    for k,x in enumerate(txu):
        b=clenshaw(c,x); out[k]=b[0]-b[2]
    return out
RYD_EV=13.605693; EV=1.602176634e-12; KB=1.380649e-16
tab=None
def vh14_load(path):
    global tab
    rows=[l.split('#')[0].split() for l in open(path)]
    rows=[r for r in rows if r]
    ng,nu=int(rows[1][0]),int(rows[1][1]); lg0=float(rows[2][0]); lu0=float(rows[3][0]); st=float(rows[4][0])
    vals=np.array([float(v) for r in rows[5:] for v in r])
    tab=(vals[:ng*nu].reshape(nu,ng).T.copy(),lg0,lu0,st)
def vh14(Z,T,u):
    g,lg0,lu0,st=tab
    RYK=157807.4
    lg=np.log10(Z*Z*RYK/T)
    nu_ryd=u/(RYD_EV*EV/(KB*T))
    lu=np.log10(nu_ryd*RYK/T)
    xg=np.clip((lg-lg0)/st,0,g.shape[0]-1-1e-9); xu=np.clip((lu-lu0)/st,0,g.shape[1]-1-1e-9)
    ig=int(xg); wg=xg-ig; iu=xu.astype(int); wu=xu-iu
    return g[ig,iu]*(1-wg)*(1-wu)+g[ig+1,iu]*wg*(1-wu)+g[ig,iu+1]*(1-wg)*wu+g[ig+1,iu+1]*wg*wu
def gff(Z,T,u):
    gamma2=Z*Z*RYD_EV*EV/(KB*T)
    hum=(T>=10**2.5)&(1e-3<=gamma2<=1e3)
    out=vh14(Z,T,u)
    if hum:
        m=(u>=1e-4)&(u<=10**1.5)
        if m.any(): out[m]=hummer(Z,T,u[m])
    return out
def gbar(T,Z):
    # breakpoints: Hummer range edges and every van Hoof u node (kinks of bilinear)
    g,lg0,lu0,st=tab
    shift=np.log10((RYD_EV*EV/(KB*T))*157807.4/T)  # lu = log10(u) + shift
    nodes=lu0+st*np.arange(g.shape[1])-shift
    lo,hi=-14.0,np.log10(700.0)
    br=np.unique(np.concatenate([[lo,hi,-4.0,1.5],nodes[(nodes>lo)&(nodes<hi)]]))
    total=0.0
    xg,wg=np.polynomial.legendre.leggauss(40)
    for a,b in zip(br[:-1],br[1:]):
        x=0.5*(b-a)*xg+0.5*(a+b); u=10**x
        total+=0.5*(b-a)*np.sum(wg*gff(Z,T,u)*u*np.exp(-u)*np.log(10))
    return total
if __name__=='__main__':
    from pathlib import Path
    here=Path(__file__).resolve().parent
    vh14_load(here.parents[1]/'data/gauntff_vh14.dat')
    lines=['# Thermally averaged free-free Gaunt factor gbar(T, Z); make_gbar_ff_table.py',
           '# columns: log10 T [K], Z = 1..6']
    for i in range(401):
        lt=2.0+4.0*i/400
        lines.append(f'{lt:.2f} '+' '.join(f'{gbar(10**lt,z):.15e}' for z in range(1,7)))
    (here/'gbar_ff_table.txt').write_text('\n'.join(lines)+'\n')
