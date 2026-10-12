"""Staged 32-bound-level Cl VII atom with adopted cold constant Upsilon.

Historical record (2026-10-11 staging of the Cl VII atom). It cannot run on
the current tree: its native driver tests/ion_stages/check_excitation_atom.f90
was removed; the installed Cl VII atom is data/atomic/excitation/chianti/
nlevel_cl_7.txt, checked by tests/ion_stages/check_excitation_atoms.py.

Compare native rates/populations with the full independent statistical matrix
using source tables at actual T. Conditions match: Maxwellian electrons,
spontaneous emission, no photoexcitation, trapping, or excited ionization.
The source cold approximation is shared, so this is a numerical verification.
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
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    args=parser.parse_args();work=args.work_dir.resolve();out=work/'products/ClVII'
    path,nl,nr,nc,error=N.build_ion('ClVII','cl',7,32,output_dir=out)
    fit_result=dict(passed=bool(error<N.ERR_TARGET),levels=nl,radiative_rows=nr,collision_rows=nc,
                    maximum_transition_fit_error=float(error))
    atomic_write_text(work/'transition_fit_result.json',json.dumps(fit_result,indent=2)+'\n')
    if error>=N.ERR_TARGET or (nl,nr,nc)!=(32,496,496):
        raise ValueError(f'Cl VII source strength fit failed: {fit_result}')
    (work/'native_atom').mkdir(parents=True,exist_ok=True)
    exe=compile_physics_check('tests/ion_stages/check_excitation_atom.f90',work/'native_atom')
    dest=out/'atom_native.txt'
    subprocess.run([str(exe),str(out/'cooling_cl_7.txt'),str(path),str(dest)],check=True)
    rows=np.loadtxt(dest)
    if rows.shape!=(2403,36):raise ValueError('Incomplete native Cl VII profile')
    temperature=rows[::3,0];density=[1e-6,100.,1e6]
    raw,rad,coll=C.read_excitation_atom('cl',7)
    levels=C.selected_energies('cl',7,raw,rad,coll)
    source=full_excitation_balance(levels,rad,coll,temperature,density)
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(6,6),sharex=True,constrained_layout=True)
    errors=[];population_errors=[];saved={'temperature_K':temperature}
    for j,ne in enumerate(density):
        full,pop=source[ne];measured=rows[j::3,3];difference=measured/full-1
        if np.any(full<=0) or not np.all(np.isfinite(full)):raise ValueError('Invalid full source cooling')
        errors.append(float(np.max(np.abs(difference))))
        population_errors.append(float(np.max(np.abs(rows[j::3,4:]-pop))))
        label=r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$'
        line=axes[0].loglog(temperature,full,label=label)[0]
        axes[0].loglog(temperature,measured,'--',color=line.get_color())
        axes[1].semilogx(temperature,100*difference,color=line.get_color())
        saved[f'cooling_full_density_{j}']=full;saved[f'populations_full_density_{j}']=pop
    axes[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
    axes[0].set_title('Cl VII (solid: source; dashed: native 32 states)')
    axes[0].legend(fontsize=8);axes[1].set_xlabel(r'$T$ [K]');axes[1].set_ylabel(r'Difference [\%]')
    fig.savefig(out/'ClVII_density_profiles.pdf');plt.close(fig)
    np.savez_compressed(out/'full_density_profiles.npz',**saved)
    result=dict(passed=max(errors)<.005,points=2403,levels=nl,densities_cm3=density,
                temperatures_K=[1000,100000],maximum_relative_cooling_differences=errors,
                maximum_absolute_population_differences=population_errors,
                maximum_transition_fit_error=float(error),
                cold_approximation='Source Upsilon fixed below 9800 K; actual T in Maxwellian rates',
                physical_accuracy='Low-temperature resonances unverified; above-threshold core states excluded')
    atomic_write_text(work/'atom_result.json',json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise ValueError(f'Cl VII native excitation mismatch: {result}')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
