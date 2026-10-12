import numpy as np, sys
A='/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/data/atomic'
def rows(el):
    out={};cur=None
    for L in open(f'{A}/element_{el}.txt'):
        t=L.split()
        if not t: continue
        if t[0]=='TRANSITION': cur=int(t[1]); out[cur]={'sh':[]}
        elif t[0]=='PHOTO': out[cur]['ph']=[float(x) for x in t[1:]]
        elif t[0]=='PHINNER': out[cur]['einn']=float(t[1])
        elif t[0]=='PHSHELL': out[cur]['sh'].append([float(x) for x in t[1:]])
    return out
def vfky(E,p):
    Eth,E0,s0,ya,P,yw,y0,y1=p
    x=E/E0-y0; y=np.sqrt(x*x+y1*y1)
    return s0*((x-1)**2+yw**2)*y**(0.5*P-5.5)*(1+np.sqrt(y/ya))**(-P)
def vy(E,q):  # PHSHELL: l? Eth E0 s0 ya P yw  -> guess layout
    l,Eth,E0,s0,ya,P,yw=q
    y=E/E0; Q=5.5+l-0.5*P
    return s0*((y-1)**2+yw**2)*y**(-Q)*(1+np.sqrt(y/ya))**(-P)
if __name__=='__main__':
  for el,q in [('s',1),('s',2),('s',3),('ne',3),('o',2),('o',3),('n',2),('c',2),('ne',2)]:
    r=rows(el)[q]; Eth=r['ph'][0]
    outer=r['sh'][-1]
    print(f"{el}{q} Eth={Eth:7.3f} VFKY96 sig_th={vfky(Eth,r['ph']):8.3f} Mb  sig(Eth+2eV)={vfky(Eth+2,r['ph']):7.3f}  VY95 outer({outer[0]:.0f},{outer[1]:.2f}) sig_th={vy(Eth,outer):8.3f} Mb  sig(+2)={vy(Eth+2,outer):7.3f} einn={r['einn']}")
