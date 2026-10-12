from pathlib import Path
import json,hashlib,re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
root=Path.cwd();base=root/'md/atomic_stage_expansion/validation_2026-10-11/feii_density_cooling_checked'
state=json.loads((base/'controller/state.json').read_text());result=json.loads((base/'native/result.json').read_text())
assert state['status']=='completed' and result['passed']
for group in (state,result):
 for name,digest in group['source_sha256'].items():
  assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest,name
lines=(base/'native/registered.txt').read_text().splitlines()
rows=np.array([[float(re.sub(r'(?<=[0-9.])([+-]\d{3})$',r'E\1',v)) for v in s.split()[1:]] for s in lines if not s.startswith('#')])
oldpath=root/'md/atomic_stage_expansion/validation_2026-10-11/feii_forbidden_effect/products/with_forbidden.txt'
old=np.loadtxt(oldpath)
errors=[]
for row in rows[:20]:
 candidates=old[np.isclose(old[:,0],row[0],rtol=1e-12)&np.isclose(old[:,1],row[1],rtol=1e-12,atol=0)]
 assert len(candidates)==1
 errors.append(abs(row[2]/candidates[0,2]-1))
assert max(errors)<1e-8
pdf=base/'matching_density_profiles.pdf'
with PdfPages(pdf) as pages:
 for j,ne in enumerate([1e-6,100.,1e6,1e10]):
  curve=old[j::4];grid=rows[j:20:4];middle=rows[20+j::4]
  fig,axes=plt.subplots(2,1,figsize=(7.8,7),constrained_layout=True)
  axes[0].loglog(curve[:,0],curve[:,2],label='Previous direct 401-point profile')
  axes[0].loglog(grid[:,0],grid[:,2],'o',ms=5,label='Optimized direct solve, 5 points')
  axes[0].loglog(grid[:,0],grid[:,4],'x',ms=7,label='Actual matching-fit application')
  axes[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]');axes[0].legend(fontsize=8)
  axes[0].set_title(r'Fe II: 142 levels, 40 preserved forbidden transitions; $n_e = '+f'{ne:g}'+r'$ cm$^{-3}$')
  axes[1].semilogx(grid[:,0],100*(grid[:,4]/grid[:,2]-1),'o-',label='Application error at regular density')
  axes[1].semilogx(middle[:,0],100*(middle[:,4]/middle[:,2]-1),'x--',label=r'Application error at nearby $n_e = '+f'{middle[0,1]:.4g}'+r'$ cm$^{-3}$')
  axes[1].semilogx(middle[:,0],100*(middle[:,3]/middle[:,2]-1),'s:',label='Table error at nearby density')
  axes[1].axhline(0,color='gray',lw=.5);axes[1].set_xlabel(r'$T$ [K]');axes[1].set_ylabel('Relative difference [percent]');axes[1].legend(fontsize=8)
  pages.savefig(fig);plt.close(fig)
review=dict(controller_completed=True,source_and_artifact_hashes_verified=True,native_samples=40,table_cells=48461,previous_direct_comparisons=20,maximum_optimized_previous_direct_relative_difference=max(errors),maximum_relative_thermal_direct_energy_difference=result['maximum_relative_thermal_direct_energy_difference'],maximum_relative_table_interpolation_difference=result['maximum_relative_table_interpolation_difference'],figure_pages=4,figure_pages_visually_inspected=0,production_modified=False,physical_default_accepted=False,atomic_accuracy_verified=False,additional_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (oldpath,pdf)})
(base/'root_review.json').write_text(json.dumps(review,indent=2)+'\n')
print(json.dumps({k:v for k,v in review.items() if k!='additional_hashes'},indent=2))
