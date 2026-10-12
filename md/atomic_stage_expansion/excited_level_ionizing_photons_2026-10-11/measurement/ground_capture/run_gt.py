import sys,glob,numpy as np
sys.path.insert(0,'.')
from badgt import read,interp,B
# ion key (MoCHII transition label = recombined ion stage) -> rr file, dr files, ground-term share rule for CA files
CASES={
 's1':('rr/eab12#p_s1ic',['dr/eab12#p_s1ic33','dr/eab12#p_s1ic34'],dict(rr=1.0,dr=1.0)),   # 4S + p -> p^4 3P only (cfp); DR spin-allowed only to 3P (LS)
 's2':('rr/jk12#si_s2ic',['dr/jk12#si_s2ic33','dr/jk12#si_s2ic34'],dict(rr=1/3,dr=None)),
 's3':('rr/saan05#al_s3ic',['dr/nrb00#al_s3ic33','dr/saan00#al_s3ic34'],None),
 'ne3':('rr/nrb05#n_ne3ic',['dr/nrb00#n_ne3ic22','dr/nrb00#n_ne3ic23'],None),
 'ne2':('rr/nrb05#o_ne2ic',['dr/oiz00#o_ne2ic22','dr/oiz00#o_ne2ic23'],None),
 'o2':('rr/nrb05#c_o2ic',['dr/nrb00#c_o2ic22','dr/nrb00#c_o2ic23'],None),
 'o3':('rr/nrb05#b_o3ic',['dr/za00#b_o3ic22','dr/za00#b_o3ic23'],None),
 'n2':('rr/nrb05#b_n2ic',['dr/za00#b_n2ic22','dr/za00#b_n2ic23'],None),
 'c2':('rr/nrb05#be_c2ic',['dr/jc00#be_c2ic12','dr/jc00#be_c2ic23','dr/nrb00#be_c2ic22'],None),
}
def ground_levels(d):
    lv=d['lev']; g=lv[1]
    return [k for k,v in lv.items() if v[0]==g[0] and v[1]==g[1] and v[2]==g[2]]
def parent_ground(d):
    p=d['par']; g=p[1]
    return {k:2*v[3]+1 for k,v in p.items() if v[0]==g[0] and v[1]==g[1] and v[2]==g[2]}
def gt_rate(path,T,rule,kind,prti_avg):
    d=read(f'{B}/{path}.dat.gz')
    pw=parent_ground(d) if prti_avg else {1:1.0}
    tot=0;w=0
    for p,wt in pw.items():
        if p not in d['prti']: continue
        blk=d['prti'][p]
        if d['lev']:
            gl=ground_levels(d)
            r=sum(interp(d['TE'],blk['indx'][k],T) for k in gl if k in blk['indx'])
        else:
            share=rule[kind]
            if kind=='rr':
                r=share*interp(d['TE'],blk['nl'][(3,1)],T)
            else:
                r=(share if share is not None else 0.0)*(interp(d['TE'],blk['indx'][1],T) if 1 in blk['indx'] else 0)
        tot+=wt*r; w+=wt
    return tot/w, d
if __name__=='__main__':
    for T in (5e3,1e4,2e4):
        print(f'T = {T:.0f} K')
        for k,(rr,drs,rule) in CASES.items():
            for avg in (False,True):
                a_rr,d=gt_rate(rr,T,rule,'rr',avg)
                a_dr=sum(gt_rate(f,T,rule,'dr',avg)[0] for f in drs)
                print(f'  {k:4s} {"termavg" if avg else "PRTI=1 "}  RR_gt {a_rr:.4e}  DR_gt {a_dr:.4e}  sum {a_rr+a_dr:.4e}')
