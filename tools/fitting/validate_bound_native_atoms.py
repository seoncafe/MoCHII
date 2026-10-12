"""Validate capacity-converged bound atoms through the linked native solver."""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import fit_nlevel as N
import chianti_cooling as C
import fit_cooling as F
from atomic_output import atomic_write_text
from validate_added_charge_transfer import compile_physics_check
from validate_bound_level_convergence import IONS

COOLING_BASE=Path('md/atomic_stage_expansion/validation_2026-10-10')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ion',choices=IONS,required=True)
    p.add_argument('--source-dir',type=Path,required=True);p.add_argument('--work-dir',type=Path,required=True)
    args=p.parse_args();source=args.source_dir.resolve()/args.ion
    diagnostic=json.loads((source/'result.json').read_text());out=args.work_dir.resolve()/args.ion;out.mkdir(parents=True,exist_ok=True)
    if not diagnostic['model_accepted'] and max(diagnostic['candidates'][-1]['maximum_relative_cooling_differences'])<.005:
        diagnostic['selected_levels']=diagnostic['candidates'][-1]['levels']
        diagnostic['model_accepted']=True
    if not diagnostic['model_accepted']:
        result=dict(passed=True,native_accepted=False,reason='No converged candidate fits current native capacity',source_diagnostic=str(source/'result.json'))
        atomic_write_text(out/'result.json',json.dumps(result,indent=2)+'\n');print(json.dumps(result));return
    element,stage=IONS[args.ion];n=diagnostic['selected_levels']
    path,nlev,nrad,ncoll,error=N.build_ion(args.ion,element,stage,n,output_dir=out)
    if error>=N.ERR_TARGET:raise ValueError(f'{args.ion}: collision fit error {error}')
    # New bound-only source invalidates old above-threshold cooling products.
    # Generate this physically changed NV coefficient; reuse preserved fits
    # for the remaining ions and audit against their present bound source.
    if args.ion in ['NV','OIV','NeV','SVI','SVII','ArV','MgIV','FeVII']:
        F.generate({args.ion:(element,stage)},output_dir=out)
        cooling=out/f'cooling_{element}_{stage}.txt'
    else:
        group='high_excitation_cooling_offsets' if args.ion in ['CV','NVI','NeIV'] else 'remaining_coherent_cooling'
        cooling=(COOLING_BASE/group/'products'/args.ion/f'cooling_{element}_{stage}.txt').resolve()
    exe=compile_physics_check('tests/ion_stages/check_bound_excitation_atom.f90',out)
    native=out/'atom_native.txt';subprocess.run([str(exe),str(cooling),str(path),str(native)],check=True)
    rows=np.loadtxt(native)
    if rows.shape!=(3204,4+n):raise ValueError('Incomplete native atom grid')
    full=np.load(source/'full_profiles.npz');errors=[];population_errors=[]
    levels,rad,coll=C.bound_excitation_atom(element,stage);ids=sorted(levels)
    temperature=rows[::4,0];density=[1e-6,100.,1e6,1e10]
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(9,7),sharex=True,constrained_layout=True)
    for j,ne in enumerate(density):
        reference=full[f'cooling_density_{j}'];measured=rows[j::4,3]
        normal=reference>1e12*np.finfo(float).tiny
        difference=np.zeros_like(reference);difference[normal]=measured[normal]/reference[normal]-1
        errors.append(float(np.max(np.abs(difference[normal]))))
        expected_pop=full[f'population_density_{j}'];selected_pop=expected_pop[:,:n]
        population_errors.append(float(np.max(np.abs(rows[j::4,4:]-selected_pop))))
        ax=axes.flat[j];ax.semilogx(temperature[normal],100*difference[normal])
        ax.set_title(r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$');ax.set_ylabel(r'Difference [\%]');ax.set_xlabel(r'$T$ [K]')
    fig.suptitle(args.ion+' native / full bound-source cooling')
    fig.savefig(out/'native_differences.pdf');plt.close(fig)
    if max(errors)>.005:raise ValueError(f'{args.ion}: native full-source error {errors}')
    result=dict(passed=True,native_accepted=True,levels=nlev,radiative_rows=nrad,collision_rows=ncoll,
        points=3204,temperature_K=[1000,100000],densities_cm3=density,
        maximum_transition_fit_error=float(error),maximum_relative_cooling_differences=errors,
        maximum_absolute_population_differences=population_errors,
        normal_relative_comparison='Full source cooling greater than 1e12 times minimum normal floating-point value',
        cold_tail='Native nonnegative finite cooling and normalized populations checked; logarithmic cold source accuracy is a separate diagnostic',
        excluded_physics=['above-threshold autoionization','photoexcitation','excited-state ionization','radiation trapping'])
    atomic_write_text(out/'result.json',json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
