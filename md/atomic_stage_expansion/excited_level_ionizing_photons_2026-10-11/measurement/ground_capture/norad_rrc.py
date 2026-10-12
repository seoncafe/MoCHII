import sys,re,numpy as np
def read(path):
    txt=open(path).read().split('\n')
    i=[k for k,l in enumerate(txt) if l.strip()=='photoionization: opacity values'][0]
    T=10**np.linspace(1,9,81)
    toks=' '.join(txt[i+1:]).split()
    out={};j=0
    while j<len(toks):
        t=toks[j]
        if re.match(r'^\d{8}\.\d{4}$',t):
            try:
                be=float(toks[j+1]); v=np.array([float(x) for x in toks[j+2:j+83]])
            except Exception: break
            out[t]=(be,v); j+=83
        else: j+=1
    return T,out
if __name__=='__main__':
    T,o=read(sys.argv[1])
    items=sorted(o.items(),key=lambda x:x[1][0])[:6]
    for k,(be,v) in items:
        print(k,be,' '.join(f'T={t:.0e}:{np.exp(np.interp(np.log(t),np.log(T),np.log(np.maximum(v,1e-300)))):.4e}' for t in (5e3,1e4,2e4)))
