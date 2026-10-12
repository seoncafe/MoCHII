import numpy as np, sys
TRAPZ = getattr(np, "trapezoid", None) or np.trapz
sys.path.insert(0,'.')
from sig import rows, vfky
import os
TB=os.environ.get('REFDATA_DIR','refdata_dl')
RY=13.605693; me,c,k,eV=9.1093837e-28,2.99792458e10,1.380649e-16,1.602176634e-12
G={}
for L in open('/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/data/atomic/ground_term_weights.txt'):
    t=L.split()
    if t and t[0]=='GROUND': G[(t[1],int(t[3]))]=int(t[4])
def op(name):
    L=open(f'{TB}/{name}.txt').read().split('\n')
    hdr=L[11].split(); Eb=-float(hdr[5])
    d=np.array([[float(x) for x in l.split()] for l in L[12:] if len(l.split())==2 and 'E' in l])
    return Eb,d[:,0],d[:,1]
def milne(E,sig,Eth,T,gl,gu):
    kT=k*T/eV
    m=E>=Eth
    E=E[m];s=sig[m]
    return np.sqrt(2/np.pi)/(c*c*(me*k*T)**1.5)*gl/gu*TRAPZ((E*eV)**2*s*1e-18*np.exp(-(E-Eth)/kT),E*eV)
# g_n/g_+ here equals MoCHII's g_n/(2 g_+) with its own prefactor (checked: same alpha_1 as the gate)
IONS=['s1','s2','s3','ne3','ne2','o2','o3','n2','c2','s4']
def key(n):
    import re; m=re.match(r'([a-z]+)(\d)',n); return m.group(1),int(m.group(2))
if __name__=='__main__':
    for n in IONS:
        el,q=key(n); r=rows(el)[q]; Eth=r['ph'][0]; gl,gu=G[(el,q)],G[(el,q+1)]
        Eb,Ep,s=op(n)
        Ee=Ep*RY - Ep[0]*RY + Eth   # shift OP grid so its first point sits on the experimental IP
        E2=Eth+np.linspace(0,40*8.617e-5*2e4,40001)
        out=[]
        for T in (5e3,1e4,2e4):
            aop=milne(Ee,s,Eth,T,gl,gu)
            av=milne(E2,vfky(E2,r['ph']),Eth,T,gl,gu)
            out.append((T,aop,av))
        print(f"{n:4s} Eth={Eth:7.3f} OP_first={Ep[0]*RY:7.3f} OP_bind={Eb*RY:7.3f}  sigOP(first)={s[0]:.3f} VFKY(th)={vfky(Eth,r['ph']):.3f} | "+
              '  '.join(f"T={T:.0e}: OP {a:.3e} VF {b:.3e}" for T,a,b in out))
