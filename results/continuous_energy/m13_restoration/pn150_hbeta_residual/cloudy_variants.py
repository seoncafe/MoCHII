#!/usr/bin/env python3
"""Cloudy PN150 variants: R_B(H), R_B(He+), R_B(He2+) per Q(H) on MoCHII's
alpha_B, L(Hbeta), and the photons in metal lines shortward of 911.75 A
(intrinsic list), for base names on the command line (needs print line
faint -5 log and save overview <base>.ovr)."""
import sys, re, numpy as np
sys.path.insert(0,'../../../../tools/python'); import mochii_rates as mr
Q=5.416e47; h_c=6.62607015e-27*2.99792458e10
skip=('H','He','Ca','Inci','Q(H)','nFnu','nInu','TotH','Crst','CaBo','Blnd','ROSA','SX')
pat=re.compile(r'([A-Z][a-z ] {0,1}[0-9]{1,2}[a-z ]?|[A-Za-z][A-Za-z0-9 ]{3})\s+(\d+\.\d+)A\s+(-?\d+\.\d+)\s+(\d+\.\d+)')
for tag in sys.argv[1:]:
    o=np.loadtxt(tag+'.ovr'); dep,T,hden,ne,HII,HeII,HeIII=o[:,0],o[:,1],o[:,3],o[:,4],o[:,7],o[:,9],o[:,10]
    r=1e17+dep; e=np.concatenate([[1e17],0.5*(r[1:]+r[:-1]),[r[-1]+0.5*(r[-1]-r[-2])]]); dV=4/3*np.pi*(e[1:]**3-e[:-1]**3)
    f=lambda fn: np.array([fn(t) for t in T])
    RH=np.sum(ne*hden*HII*f(mr.alphaB_HII)*dV); R2=np.sum(ne*hden*0.1*HeII*f(mr.alphaB_HeII)*dV); R3=np.sum(ne*hden*0.1*HeIII*f(mr.alphaB_HeIII)*dV)
    txt=open(tag+'.out').read(); sec=txt[txt.index('Intrinsic line intensities'):txt.index('Emergent line intensities')]
    n=0; hb=None
    for m in pat.finditer(sec):
        lab=m.group(1).strip(); wl=float(m.group(2)); lL=float(m.group(3))
        if lab=='H  1' and abs(wl-4861.32)<0.01 and hb is None: hb=10**lL
        if 10<wl<911.75 and lab.split()[0] not in skip: n+=10**lL*wl*1e-8/h_c
    print(f'{tag:8s} R_B(H)/Q {RH/Q:.4f} He+ {R2/Q:.4f} He2+ {R3/Q:.4f}  L(Hb) {hb:.4e}  metal EUV lines {n/Q:.4f} Q')
