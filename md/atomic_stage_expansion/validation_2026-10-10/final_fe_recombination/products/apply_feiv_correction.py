"""Apply the verified Fe3+ ground-level RR/DR correction only."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]
sys.path.insert(0,str(ROOT/'tools/fitting'))
from zhang_fe_recombination import dr_fit,rr_fit,read_ground_dr_grid,ARCHIVE
from atomic_output import atomic_write_text
work=Path(__file__).parent
production=ROOT/'data/atomic/element_fe.txt'
before=production.read_text()
saved=work/'element_fe.txt.before'
if saved.exists(): raise RuntimeError('Correction already prepared; inspect preserved products')
saved.write_text(before)
def rates_at_transition(text,q):
    part=text.split(f'TRANSITION {q}\n',1)[1].split('TRANSITION ',1)[0]
    return [s for s in part.splitlines() if s.startswith(('RR ','RR2 ','DR '))]
selected=rates_at_transition((work/'generated/element_fe.txt').read_text(),3)
old=rates_at_transition(before,3)
assert len(selected)==len(old)==2
new=before
for a,b in zip(old,selected): new=new.replace(a,b,1)
note='# Fe3+ ground-level RR/DR: Zhang et al. 2026 ApJS 285:50, Tables 6/7, S1.\n# Valid 180--9e7 K; zero density; no metastable-parent capture or DR suppression.\n# Analytic extrapolation outside range is uncertified. DR fit differs from ADF09\n# by up to 4.718 percent on the 18 inspected nodes within the published range.\n'
new=new.replace('TRANSITION 3\n','TRANSITION 3\n'+note,1)
for q in (1,2): assert rates_at_transition(new,q)==rates_at_transition(before,q)
restored=new.replace(note,'')
for a,b in zip(selected,old): restored=restored.replace(a,b,1)
assert restored==before
atomic_write_text(production,new)
# Full recombination coefficient curves: same density, initial parent and domain.
fig,axes=plt.subplots(2,2,figsize=(9,7),constrained_layout=True)
for q,ax in zip(range(3,7),axes.flat):
    ts=np.geomspace(20*q*q,1e7*q*q,601)
    ax.loglog(ts,[dr_fit(q,float(t)) for t in ts],label='Published S1 DR')
    ax.loglog(ts,[rr_fit(q,float(t)) for t in ts],label='Published RR',ls='--')
    grid=read_ground_dr_grid(ARCHIVE,q)
    valid=[(t,a) for t,a in grid if 20*q*q<=t<=1e7*q*q]
    ax.loglog(*zip(*valid),'o',ms=3,label='Author ADF09 ground total')
    ax.set_title('Recombining Fe'+str(q)+'+')
    ax.set_xlabel('Electron temperature (K)')
    ax.set_ylabel(r'Coefficient (cm$^3$ s$^{-1}$)')
    ax.grid(alpha=.25); ax.legend(fontsize=7)
fig.savefig(work/'verified_fe_recombination.pdf'); plt.close(fig)
a,eta=map(float,old[0].split()[1:])
terms=old[1].split(); n=int(terms[1]); c=list(map(float,terms[2:2+n])); e=list(map(float,terms[2+n:]))
rows=[]
for t in [1000.,10000.,100000.]:
    old_rr=a*(t/1e4)**(-eta)
    old_dr=sum(x*np.exp(-y/t) for x,y in zip(c,e))/t**1.5
    rr,dr=rr_fit(3,t),dr_fit(3,t)
    rows.append({'temperature_K':t,'old_rr':old_rr,'old_dr':old_dr,'new_rr':rr,'new_dr':dr,'total_ratio':(rr+dr)/(old_rr+old_dr)})
result={'passed':True,'changed_transition':3,'production_stages':4,'unchanged_other_records':True,'before_sha256':hashlib.sha256(before.encode()).hexdigest(),'after_sha256':hashlib.sha256(new.encode()).hexdigest(),'rates_cm3_s':rows,'stale_outputs':'Previous Fe IV recombination-dependent fractions and cooling require recalculation for reuse.'}
(work/'production_correction.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
