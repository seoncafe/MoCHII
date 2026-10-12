"""Measure O IV source radiation-frequency consistency without modifying data."""
from pathlib import Path
import hashlib,json,math,re,sys
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'tools/fitting'))
import chianti_cooling as C
HERE=Path(__file__).resolve().parent

def main():
 stem=Path(C.ion_dir('o',4))/'o_4'
 raw=C.read_elvlc(str(stem)+'.elvlc')
 selected,radiation,collisions=C.bound_excitation_atom('o',4)
 parity={}
 for line in Path(str(stem)+'.elvlc').read_text().splitlines():
  if line.strip()=='-1':break
  index=int(line[:7]);orbitals=re.findall(r'(\d+)([spdfgh])(\d*)',line[7:37])
  if not orbitals:raise ValueError('Unparsed configuration')
  parity[index]=sum((int(n) if n else 1)*'spdfgh'.index(o) for _,o,n in orbitals)%2
 e1=[];unmatched=[];non_e1=[]
 for lower,upper,wavelength,gf,rate in radiation:
  if gf<=0 or parity[lower]==parity[upper]:
   non_e1.append([lower,upper,gf,rate]);continue
  gap=selected[upper][1]-selected[lower][1]
  if gap<=0 or rate<=0:raise ValueError('Invalid radiative direction')
  factor=(abs(wavelength)*gap/1e8)**3
  record={'lower':lower,'upper':upper,'wavelength_angstrom':wavelength,'source_gf':gf,'source_A_s_inverse':rate,'selected_gap_cm_inverse':gap,'conditional_A_frequency_scaling_factor':factor}
  if upper>20:
   inferred=1e8*math.sqrt(rate*raw[upper][0]/(6.6703e15*gf))
   candidates={'theoretical':raw[upper][2]-raw[lower][2],'observed':raw[upper][1]-raw[lower][1],'wavelength':1e8/abs(wavelength)}
   errors={label:abs(value/inferred-1) for label,value in candidates.items() if value>0}
   label=min(errors,key=errors.get);record.update(inferred_Einstein_frequency_cm_inverse=inferred,minimum_available_frequency_relative_difference=errors[label],closest_frequency_source=label)
   if errors[label]>.01:unmatched.append(record)
  e1.append(record)
 result={'passed':True,'physical_model_accepted':False,'level_count':len(selected),'radiative_count':len(radiation),'collision_count':len(collisions),'source_nodes':sum(len(row['xs']) for row in collisions),'opposite_parity_positive_gf_count':len(e1),'other_radiative_count':len(non_e1),'higher_E1_count':sum(row['upper']>20 for row in e1),'higher_E1_frequency_mismatches_above_one_percent':len(unmatched),'largest_conditional_wavelength_scaling_examples':sorted(e1,key=lambda row:abs(row['conditional_A_frequency_scaling_factor']-1),reverse=True)[:8],'unmatched_higher_E1':unmatched,'interpretation':'A frequency-consistency diagnostic, not an adopted correction or atomic accuracy estimate. Low upper states 3--20 have a different A source; the original higher-state gf/A calculation and its conversion must be identified.','source_sha256':{str(Path(str(stem)+'.'+extension)):hashlib.sha256(Path(str(stem)+'.'+extension).read_bytes()).hexdigest() for extension in ('elvlc','wgfa','scups')}}
 target=HERE/'result.json'
 if target.exists():raise FileExistsError('Preserve completed diagnostic')
 target.write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if not isinstance(v,(list,dict))})
if __name__=='__main__':main()
