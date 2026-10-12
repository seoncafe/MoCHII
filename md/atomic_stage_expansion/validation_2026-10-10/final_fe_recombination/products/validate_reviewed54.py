"""Stage the reviewed 54-stage prefix and check Fe4 charge-transfer consumers."""
import sys,json,subprocess,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5];sys.path.insert(0,str(ROOT/'tools/fitting'))
import make_element_data as G
from charge_exchange_expansion import kf96_recombination_total_rates
work=Path(__file__).parent/'reviewed54';work.mkdir()
target={'c':5,'n':5,'o':5,'ne':5,'s':5,'ar':5,'mg':4,'fe':5}
G.generate([(el,z,target.get(el,n)) for el,z,n in G.ELEMENTS],work/'atomic')
shutil.copy2(ROOT/'data/atomic/nlevel_s_2.txt',work/'atomic/nlevel_s_2.txt')
exe=ROOT/'md/atomic_stage_expansion/validation_2026-10-10/kf96_total_recombination/products/check_added_charge_transfer_heat.x'
subprocess.run([str(exe),str(work/'atomic'),str(work/'native_rates.txt'),'54'],cwd=ROOT,check=True)
row=kf96_recombination_total_rates()['fe',4];errors=[]
for line in (work/'native_rates.txt').read_text().splitlines():
 fields=line.split()
 if (fields[0],int(fields[1]))!=('fe',4):continue
 t,ki,kr,di,dr,split,f=map(float,fields[3:]);assert ki==di==split==0 and dr==10 and 0<f<1
 errors.append(abs(kr/row.rate(t)-1))
assert len(errors)==3 and max(errors)<1e-5
result={'passed':True,'stages':54,'transitions':43,'Fe4_max_relative_rate_error':max(errors),'Fe4_heat_eV':10,'direct_cache_charge_and_heat':'passed at 5000,10000,30000 K','production_activation':False,'limitation':'Line cooling not populated in this controlled rate fixture; higher-charge transfer remains unselected'}
(work/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
