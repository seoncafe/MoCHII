"""Badnell adf48 (RR) / adf09 (DR) partial rates into the ground term of the recombined ion."""
import gzip, re, numpy as np, sys
import os
B=os.environ.get('BADNELL_DIR','badnell_dl')
LV=re.compile(r'^\s*(\d+)\s+(\d+)\s+([0-9A-Z ]+?)\s+\(\s*(\d+)\)\s*(\d+)\(\s*([\d.]+)\)\s+(-?[\d.]+)')
CF=re.compile(r'^\s*(\d+)\s+(\d+)\s+([0-9A-Z ]+?)\s+\(\s*([\d.]+)\)\s+(-?[\d.]+)')
PL=re.compile(r'^\s*(\d+)\s+([0-9A-Z ]+?)\s+\(\s*(\d+)\)\s*(\d+)\(\s*([\d.]+)\)\s+(-?[\d.]+)')
def fl(s): return float(s)
def read(path):
    L=gzip.open(path,'rt',errors='replace').read().splitlines()
    d=dict(lev={},cfg={},par={},prti={})
    sect=None; i=0
    while i<len(L):
        s=L[i]
        if 'PARENT LEVEL INDEXING' in s: sect='par'
        elif 'PARENT CONFG' in s: sect=None
        elif 'RESOLVED LEVEL INDEXING' in s: sect='lev'
        elif 'AVERAGED CONFG INDEXING' in s: sect='cfg'
        elif 'NL-SHELL' in s: sect='nl'; d['nl']={}
        elif 'N-SHELL' in s: sect=None
        elif re.match(r'^\s*PRTI=',s): break
        elif sect=='par':
            m=PL.match(s)
            if m: d['par'][int(m.group(1))]=(m.group(2).strip(),int(m.group(3)),int(m.group(4)),float(m.group(5)))
        elif sect=='lev':
            m=LV.match(s)
            if m: d['lev'][int(m.group(1))]=(m.group(3).strip(),int(m.group(4)),int(m.group(5)),float(m.group(6)),float(m.group(7)))
        elif sect=='cfg':
            m=CF.match(s)
            if m: d['cfg'][int(m.group(1))]=(m.group(3).strip(),float(m.group(4)),float(m.group(5)))
        elif sect=='nl':
            t=s.split()
            if len(t)>=3 and all(x.isdigit() for x in t[:3]): d['nl'][int(t[0])]=(int(t[1]),int(t[2]))
        i+=1
    cur=None; sub=None; prtf=None
    while i<len(L):
        s=L[i]
        m=re.match(r'^\s*PRTI=\s*(\d+)',s)
        if m:
            cur=int(m.group(1)); d['prti'][cur]=dict(indx={},nl={}); sub='indx'; i+=1; continue
        if re.match(r'^C-{5,}',s) or s.startswith('C '): break
        if 'INDX TE=' in s:
            te=[fl(x) for x in s.split('TE=')[1].split()]+[fl(x) for x in L[i+1].split()[2:]]
            d['TE']=np.array(te); i+=2; continue
        m=re.search(r'PRTF=\s*(\d+)',s)
        if m: prtf=int(m.group(1)); sub=None; i+=1; continue
        t=s.split()
        if t and t[0]=='ILREP': sub='ilrep'; i+=1; continue
        if t and t[0]=='INREP': sub='inrep'; i+=1; continue
        if len(t)>=11 and t[0].isdigit() and 'E' in t[1]:
            v=[fl(x) for x in t[1:]]+[fl(x) for x in L[i+1].split()]
            k=int(t[0])
            if sub=='indx': d['prti'][cur]['indx'][k]=np.array(v)
            elif sub=='ilrep' and prtf==1 and k in d.get('nl',{}): d['prti'][cur]['nl'][d['nl'][k]]=np.array(v)
            i+=2; continue
        i+=1
    return d
def interp(TE,v,T):
    return float(np.exp(np.interp(np.log(T),np.log(TE),np.log(np.maximum(v,1e-300)))))
