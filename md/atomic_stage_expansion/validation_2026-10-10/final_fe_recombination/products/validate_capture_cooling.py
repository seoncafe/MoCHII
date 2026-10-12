"""Check captured electron cooling using analytic Maxwellian energy moments."""
import sys,json,math,subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[5]
sys.path.insert(0,str(ROOT/'tools/fitting'))
from validate_added_charge_transfer import compile_physics_check
from zhang_fe_recombination import RR_FITS,DR_FITS,rr_fit
work=Path(__file__).parent
exe=compile_physics_check(str((work/'check_fe_capture_cooling.f90').relative_to(ROOT)),work)
subprocess.run([str(exe),str(work/'atomic'),str(work/'capture_cooling.txt')],cwd=ROOT,check=True)
samples=np.loadtxt(work/'capture_cooling.txt');assert samples.shape==(301,9)
errors=[]
for t,computed,*fractions in samples:
 expected=0.
 for q in range(1,7):
  a,b,t0,t1,c,t2=RR_FITS.get(q,(1e-11,.5,1e4,1e4,0,1))
  bt=b+c*math.exp(-t2/t);db=c*math.exp(-t2/t)*t2/t
  s0=math.sqrt(t/t0);s1=math.sqrt(t/t1)
  rr=a/(s0*(1+s0)**(1-bt)*(1+s1)**(1+bt))
  slope=-.5-(1-bt)*s0/(2*(1+s0))-(1+bt)*s1/(2*(1+s1))+db*math.log((1+s0)/(1+s1))
  cs,es=DR_FITS.get(q,((),()))
  moment=sum(ci*ei*math.exp(-ei/t) for ci,ei in zip(cs,es))/t**1.5
  beta=1.380649e-16*(t*rr*(1.5+slope)+moment)
  expected+=10.*100.*1e-4*fractions[q]*beta
 errors.append(abs(computed/expected-1))
assert max(errors)<1e-6,max(errors)
r={'passed':True,'samples':len(errors),'maximum_relative_error':max(errors),'temperature_K':[1000,100000],'reference':'Analytic Maxwellian RR and resonance-energy-weighted DR moments','excluded':'CI and CX disabled in coefficient fixture only'}
(work/'capture_cooling_result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
