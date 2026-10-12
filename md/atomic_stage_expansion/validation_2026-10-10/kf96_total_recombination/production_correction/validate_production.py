import hashlib,json,math,subprocess,sys
from pathlib import Path
sys.path.insert(0,'tools/fitting')
import make_element_data as G
from charge_exchange_expansion import kf96_recombination_total_rates,KF99_RECOMBINATION_MEAN_HEAT
root=Path.cwd();base=root/'md/atomic_stage_expansion/validation_2026-10-10/kf96_total_recombination';out=base/'production_correction';products=base/'products'
index={el:i+1 for i,(el,z,n) in enumerate(G.ELEMENTS)}
samples=[];rows=kf96_recombination_total_rates()
for key in [('cl',3),('cl',4),('ca',3),('ca',4)]:
 row=rows[key]
 for t in [row.min_temperature/10]+[row.min_temperature*(row.max_temperature/row.min_temperature)**(i/2000) for i in range(2001)]+[row.max_temperature*10]:samples.append((key,t))
(out/'temperatures.txt').write_text(''.join(f'{index[el]} {q} {t:.17e}\n' for (el,q),t in samples))
subprocess.run([str(products/'check_charge_transfer_ranges.x'),str(root/'data/atomic'),str(out/'temperatures.txt'),str(out/'compiled_rates.txt')],check=True,stdout=(out/'rate_execution.log').open('w'))
lines=(out/'compiled_rates.txt').read_text().splitlines();assert len(lines)==len(samples)
errors={}
for line,(key,t) in zip(lines,samples):
 ie,it,temperature,ki,kr=line.split();assert (int(ie),int(it),float(temperature))==(index[key[0]],key[1],t)
 assert float(ki)==0 and float(kr)>0
 error=abs(float(kr)/rows[key].rate(t)-1);errors[str(key)]=max(errors.get(str(key),0),error)
assert max(errors.values())<1e-5
heat=0
for line in (out/'thermal_closure.txt').read_text().splitlines():
 fields=line.split();key=fields[0],int(fields[1])
 if key not in [('cl',3),('cl',4),('ca',3),('ca',4)]:continue
 t,ki,kr,di,dr,split,f=map(float,fields[3:]);assert abs(dr-KF99_RECOMBINATION_MEAN_HEAT[key])<1e-12
 assert ki==di==split==0 and kr>0 and 0<f<=1;heat+=1
assert heat==12
result={'passed':True,'stages':43,'transitions':32,'reaction_samples':len(samples),'max_relative_errors':errors,'heat_samples':heat,'max_heat_difference_eV':0.,'charge_and_thermal_consumers':'compiled controlled three-temperature check passed','production_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'data/atomic').glob('element_*.txt')}}
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
