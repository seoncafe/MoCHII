import os
import numpy as np, sys
TRAPZ = getattr(np, "trapezoid", None) or np.trapz
sys.path.insert(0,'.')
me,c,k,eV=9.1093837e-28,2.99792458e10,1.380649e-16,1.602176634e-12
d=np.array([[float(x) for x in l.split()[:2]] for l in open(os.environ.get('REFDATA_DIR','refdata_dl')+'/NeIII_ALS.dat',errors='replace') if len(l.split())>=2 and l.split()[0][0].isdigit() and '.' in l.split()[0]])
E,s=d[:,0],d[:,1]
print('range',E[0],E[-1], 'mean sigma 63.45-64.45 eV', s[(E>63.45)&(E<64.45)].mean(), ' 64.45-66.45', s[(E>64.45)&(E<66.45)].mean(),' 66.5-68', s[(E>66.5)&(E<68)].mean())
Eth=63.45
def milne(E,s,T,Eth):
    kT=k*T/eV; m=E>=Eth
    return np.sqrt(2/np.pi)/(c*c*(me*k*T)**1.5)*(9/4)*TRAPZ((E[m]*eV)**2*s[m]*1e-18*np.exp(-(E[m]-Eth)/kT),E[m]*eV)
for T in (5e3,1e4,2e4):
    print(T, 'Milne(ALS, measured, term g=9/4, above 63.45 eV) %.3e'%milne(E,s,T,Eth), '  frac of Boltzmann weight covered', 1-np.exp(-(E[-1]-Eth)/(k*T/eV)))
