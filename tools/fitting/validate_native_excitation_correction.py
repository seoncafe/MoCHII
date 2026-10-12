"""Check corrected collision interpolation for an existing native atom.

Compare the unchanged installed level set with its direct source matrix.
This validates interpolation and representation, not truncation against the
full source atom or other omitted physical processes.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import chianti_cooling as C
import fit_nlevel as N
from atomic_output import atomic_write_text
from validate_added_charge_transfer import compile_physics_check
from validate_excitation_density_profiles import full_excitation_balance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ion',choices=['CI','SIII','SI'],required=True)
    p.add_argument('--work-dir',type=Path,required=True)
    args=p.parse_args();el,st,n=('s',1,5) if args.ion=='SI' else N.IONS[args.ion]
    out=args.work_dir.resolve()/args.ion;out.mkdir(parents=True,exist_ok=True)
    path,nlev,nrad,ncoll,error=N.build_ion(args.ion,el,st,n,output_dir=out)
    if error>=N.ERR_TARGET:raise ValueError(f'Unaccepted strength fit: {error}')
    levels,rad,coll=C.bound_excitation_atom(el,st)
    ids=sorted(levels)[:n];lev={i:levels[i] for i in ids}
    rr=[r for r in rad if r[0] in lev and r[1] in lev]
    cc=[t for t in coll if t['ll'] in lev and t['ul'] in lev]
    temperature=np.geomspace(1e3,1e5,801);density=[1e-6,100.,1e6,1e10]
    source=full_excitation_balance(lev,rr,cc,temperature,density)
    exe=compile_physics_check('tests/ion_stages/check_bound_excitation_atom.f90',out)
    native=out/'atom_native.txt'
    cooling=Path('data/atomic')/f'cooling_{el}_{st}.txt'
    if not cooling.exists():
        cooling=out/f'cooling_{el}_{st}.txt'
        if not cooling.exists():
            raise FileNotFoundError('Prepare the staged cooling fit before native validation')
    subprocess.run([str(exe),str(cooling),str(path),str(native)],check=True)
    rows=np.loadtxt(native)
    if rows.shape!=(3204,4+n):raise ValueError('Incomplete native profile')
    errors=[];pop_errors=[];saved={'temperature_K':temperature}
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(9,7),sharex=True,constrained_layout=True)
    for j,ne in enumerate(density):
        full,pop=source[ne];measured=rows[j::4,3]
        difference=measured/full-1
        errors.append(float(np.max(np.abs(difference))))
        pop_errors.append(float(np.max(np.abs(rows[j::4,4:]-pop))))
        saved[f'cooling_density_{j}']=full;saved[f'populations_density_{j}']=pop
        ax=axes.flat[j];ax.semilogx(temperature,100*difference)
        ax.set_title(r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$');ax.set_ylabel(r'Difference [\%]');ax.set_xlabel(r'$T$ [K]')
    fig.suptitle(args.ion+' native / corrected source, installed level set')
    fig.savefig(out/'native_differences.pdf');plt.close(fig)
    if max(errors)>.005:raise ValueError(f'Unaccepted native cooling fit: {errors}')
    np.savez_compressed(out/'source_profiles.npz',**saved)
    result=dict(passed=True,levels=nlev,radiative_rows=nrad,collision_rows=ncoll,
                points=3204,maximum_relative_cooling_differences=errors,
                maximum_absolute_population_differences=pop_errors,
                maximum_transition_fit_error=error,
                scope='Installed level set, four electron densities; full-atom truncation is not tested')
    atomic_write_text(out/'result.json',json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
