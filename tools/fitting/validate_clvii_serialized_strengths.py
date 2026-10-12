"""Audit serialized Cl VII strengths and the cold native type-2 limit."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import chianti_cooling as C
import fit_nlevel as N
from atomic_output import atomic_write_text


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work-dir',type=Path,required=True)
    args=p.parse_args();work=args.work_dir.resolve();out=work/'products/ClVII'
    atom=out/'nlevel_cl_7.txt';lines=atom.read_text().splitlines()
    start=next(i for i,s in enumerate(lines) if s.startswith('NUPS'))
    rows=lines[start+1:];raw,rad,coll=C.read_excitation_atom('cl',7)
    if len(rows)!=len(coll):raise ValueError('Incomplete serialized collisions')
    temperature=np.geomspace(1e3,1e5,1601);errors=[];cold_errors=[]
    for line,tr in zip(rows,coll):
        f=line.split();l,u,tt=int(f[0]),int(f[1]),int(f[2]);count=int(f[8])
        if (l,u,tt)!=(tr['ll'],tr['ul'],2):raise ValueError('Incorrect collision identity')
        stored=dict(tr,de=float(f[4]),cups=float(f[3]),ttype=tt)
        cf=np.array(f[9:9+count],dtype=float);lo,hi=float(f[5]),float(f[6]);logf=int(f[7])
        measured=N.eval_ups_fit(stored,lo,hi,cf,logf,temperature)
        errors.append(float(np.max(np.abs(measured/C.upsilon(tr,temperature)-1))))
        cold=N.eval_ups_fit(stored,lo,hi,cf,logf,[1000.,2000.,5000.,9799.])
        cold_errors.append(float(np.max(np.abs(cold/cold[0]-1))))
    result=dict(passed=max(errors)<N.ERR_TARGET and max(cold_errors)==0.,transitions=len(rows),
        temperature_samples=1601,maximum_relative_serialized_strength_error=max(errors),
        maximum_relative_cold_strength_variation=max(cold_errors),
        product_sha256=hashlib.sha256(atom.read_bytes()).hexdigest(),
        verification='Serialized evaluator follows src/nlevel_mod.f90 ups_eval; full native rate matrices independently verified by atom_result.json')
    atomic_write_text(work/'serialized_strength_result.json',json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise ValueError(result)
    # Inspect complete excited-level populations; ground normalization is
    # checked in the native driver. Each density has all 31 excited levels.
    data=np.load(out/'full_density_profiles.npz');native=np.loadtxt(out/'atom_native.txt')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(7,8),sharex=True,constrained_layout=True)
    for j,ne in enumerate([1e-6,100.,1e6]):
        pop=data[f'populations_full_density_{j}'];t=data['temperature_K']
        for level in range(1,32):
            line=axes[j].loglog(t,pop[:,level],lw=.65)[0]
            axes[j].loglog(t,native[j::3,4+level],'--',lw=.65,color=line.get_color())
        axes[j].set_ylabel(r'$p_u$');axes[j].set_title(r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$')
    axes[0].text(.03,.06,'Solid: source; dashed: native; all excited states',transform=axes[0].transAxes,fontsize=8)
    axes[-1].set_xlabel(r'$T$ [K]');fig.savefig(out/'ClVII_population_profiles.pdf');plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
