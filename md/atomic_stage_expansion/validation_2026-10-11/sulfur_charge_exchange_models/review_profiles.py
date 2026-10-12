from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
b=Path('md/atomic_stage_expansion/validation_2026-10-11/sulfur_charge_exchange_models').resolve()
s=json.loads((b/'controller/state.json').read_text());r=json.loads((b/'native/result.json').read_text())
assert s['status']=='completed' and r['passed']
for group in (s,r):
 for p,h in group['source_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
v=np.loadtxt(b/'native/mocassin_v2.txt');e=np.loadtxt(b/'native/installed.txt')
assert np.array_equal(v[:,:3],e[:,:3])
pdf=b/'charge_transfer_profiles.pdf'
with PdfPages(pdf) as out:
 for i,ne in enumerate([1e-6,100,1e6,1e10]):
  a=v[i::4];c=e[i::4];T=a[:,0]
  fig,ax=plt.subplots(3,1,figsize=(8,9.5),constrained_layout=True)
  ax[0].loglog(T,a[:,3],label='Supplied v2 5-level product atom');ax[0].loglog(T,c[:,3],'--',label='Installed 70-level product atom')
  ax[0].set_ylabel(r'$k_{\rm reverse}$ [cm$^3$ s$^{-1}$]');ax[0].legend(fontsize=8)
  ax[0].set_title(r'S$^+$ + H charge transfer: $n_e = '+f'{ne:g}'+r'$ cm$^{-3}$')
  ax[1].semilogx(T,a[:,5],label='5 levels');ax[1].semilogx(T,c[:,5],'--',label='70 levels');ax[1].set_ylabel('Mean reverse energy [eV]');ax[1].legend(fontsize=8)
  for j,label in enumerate([r'$^4S$',r'$^2D$',r'$^2P$']):
   line,=ax[2].loglog(T,a[:,7+j],label=label+' (5 levels)')
   ax[2].loglog(T,c[:,7+j],'--',color=line.get_color(),label=label+' (70 levels)')
  ax[2].set_ylabel('Product-term population');ax[2].set_xlabel(r'$T$ [K]');ax[2].legend(fontsize=7,ncol=2)
  out.savefig(fig);plt.close(fig)
idx=int(np.argmax(abs(e[:,3]/v[:,3]-1)))
review=dict(controller_completed=True,sources_and_artifacts_verified=True,native_samples=808,maximum_rate_difference_temperature_K=float(v[idx,0]),maximum_rate_difference_density_cm_inverse_3=float(v[idx,1]),figure_pages=4,figure_pages_visually_inspected=0,production_physics_modified=False,public_model_selection_implemented=False,atomic_accuracy_verified=False,figure_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest())
(b/'root_review.json').write_text(json.dumps(review,indent=2)+'\n')
print(json.dumps(review,indent=2))
