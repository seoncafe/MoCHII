"""Ground-term capture rate alpha_1 [cm^3 s^-1] of the recombined ion, five determinations
(diffuse channel 6 check, 2026-10-11).

  Milne(VFKY96)   Milne relation on the MoCHII PHOTO row (VFKY96 outer-shell fit), g_n/g_+ from
                  data/atomic/ground_term_weights.txt; reproduces the alpha_1 the channel-6 gate prints.
  Milne(TOPbase)  Milne relation on the Opacity Project ground-state cross section with resonances
                  (TOPbase, cdsweb.u-strasbg.fr/topbase/xsections.html), grid shifted so its first point
                  sits on the experimental IP.
  Nahar_unified   state-specific unified (RR+DR, R-matrix) rate into the ground state, NORAD
                  (norad.astronomy.osu.edu/<ion>/<ion>.rrc.txt; Nahar 1995 ApJS 101, 423; 1999 ApJS 120, 131;
                  2000 ApJS 126, 537).
  Badnell_RR_gt   adf48 partial RR into the ground-term levels (IC files), parent levels of the ground
                  term weighted by 2J+1. S I / S II files are configuration-resolved: the ground-term share
                  of the valence-p capture is the LS fractional-parentage share g_T |cfp|^2 / sum
                  (S I: 4S + p -> p^4 3P only, 1; S II: 3P + p -> p^3 4S, 1/3).
  Badnell_DR_gt   adf09 partial DR into the same levels; S I: all DR into 3s2 3p4 counted as 3P (LS spin
                  selection from the 4S parent); S II: DR into 4S taken as 0 (LS: 3p2(1D,1S) cores carry
                  no 4S parentage, doublet resonances).
Env: BADNELL_DIR (adf48 rr/ and adf09 dr/ files, list in ../badnell_files_used.txt),
     REFDATA_DIR (TOPbase <ion>.txt, NORAD <ion>_<ion>.rrc.txt, ne3_ne3.px.gd.txt, NeIII_ALS.dat).
"""
import sys, numpy as np
sys.path.insert(0,'.')
import op_milne as OM, run_gt as RG, norad_rrc as NR
from sig import rows, vfky
NOR={'s2':'s2','s3':'s3','o2':'o2','o3':'o3','n2':'n2','c2':'c2'}
def nahar(n,T):
    Tg,o=NR.read(f'{OM.TB}/{n}_{n}.rrc.txt')
    k,(be,v)=sorted(o.items(),key=lambda x:x[1][0])[0]
    return float(np.exp(np.interp(np.log(T),np.log(Tg),np.log(v))))
print('ion  T      Milne(VFKY96)  Milne(TOPbase)  Nahar_unified  Badnell_RR_gt  Badnell_RR+DR_gt  VF/Nahar  VF/RR  VF/(RR+DR)')
for n in ['s1','s2','s3','ne3','o2','o3','n2','c2']:
    el,q=OM.key(n); r=rows(el)[q]; Eth=r['ph'][0]; gl,gu=OM.G[(el,q)],OM.G[(el,q+1)]
    Eb,Ep,s=OM.op(n); Ee=Ep*OM.RY-Ep[0]*OM.RY+Eth
    rr,drs,rule=RG.CASES[n]
    for T in (5e3,1e4,2e4):
        E2=Eth+np.linspace(0,40*8.617e-5*2e4,40001)
        av=OM.milne(E2,vfky(E2,r['ph']),Eth,T,gl,gu); ao=OM.milne(Ee,s,Eth,T,gl,gu)
        a_rr=RG.gt_rate(rr,T,rule,'rr',True)[0]; a_dr=sum(RG.gt_rate(f,T,rule,'dr',True)[0] for f in drs)
        an=nahar(n,T) if n in NOR else float('nan')
        print(f'{n:4s} {T:6.0f} {av:12.4e} {ao:14.4e} {an:14.4e} {a_rr:14.4e} {a_rr+a_dr:16.4e}   {av/an:7.3f} {av/a_rr:6.3f} {av/(a_rr+a_dr):8.3f}')
