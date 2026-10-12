import os
import numpy as np
TRAPZ = getattr(np, "trapezoid", None) or np.trapz
RY=13.605693; me,c,k,eV=9.1093837e-28,2.99792458e10,1.380649e-16,1.602176634e-12
L=open(os.environ.get('REFDATA_DIR','refdata_dl')+'/ne3_ne3.px.gd.txt').read().split('\n')
i=[j for j,l in enumerate(L) if l.split()[:4]==['0','4','0','1']][0]
be=-float(L[i+1].split()[0]); n=int(L[i+1].split()[1])
d=np.array([[float(x) for x in l.split()] for l in L[i+3:i+3+n] if len(l.split())==4])
E=d[:,0]; s=d[:,3]
print('BE(eV)',be*RY,'points',len(E),'E0',E[0],'sig0',s[:3])
Eth=63.423  # MoCHII ETH (Ne2+ ground term)
Ee=(E-be)*RY+Eth
for T in (5e3,1e4,2e4):
    kT=k*T/eV
    a=np.sqrt(2/np.pi)/(c*c*(me*k*T)**1.5)*(5/4)*TRAPZ((Ee*eV)**2*s*1e-18*np.exp(-(Ee-Eth)/kT),Ee*eV)
    print(f'T={T:.0e}  Milne(Nahar 2019 3P2, g=5/g+=4) {a:.4e}   ground-term equivalent x9/5 {a*9/5:.4e}')
