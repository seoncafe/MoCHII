"""Check a selected N IV atom against the full supplied collision model.

Historical record (2026-10-11 staging of the N IV atom). Its native step
cannot run on the current tree: the driver tests/ion_stages/
check_niv_staging.f90 was removed. The installed N IV atoms are checked by
tests/ion_stages/check_excitation_atoms.py.

This is stationary electron-collision and spontaneous-radiation physics only.
The temperature interval is numerical; resonance accuracy remains uncertain.
"""
import argparse
import json
from pathlib import Path
import subprocess
from unittest.mock import patch

import numpy as np
import chianti_cooling as C
import fit_nlevel as N
from validate_added_charge_transfer import compile_physics_check
from validate_stage_cooling import coefficients
from fit_cooling import evaluate

ROOT=Path(__file__).resolve().parents[2]


def raw_atom():
    stem=Path(C.ion_dir('n',4))/'n_4'
    raw=C.read_elvlc(str(stem)+'.elvlc')
    radiation=C.read_wgfa(str(stem)+'.wgfa')
    collisions=C.read_scups(str(stem)+'.scups')
    return C.selected_energies('n',4,raw,radiation,collisions),radiation,collisions


def truncated_cooling(n,temperatures):
    levels,radiation,collisions=raw_atom()
    selected={i:r for i,r in levels.items() if i<=n}
    rad=[r for r in radiation if r[0]<=n and r[1]<=n]
    coll=[r for r in collisions if r['ll']<=n and r['ul']<=n]
    # Source identity was checked in full above; only the model cutoff changes.
    with patch.object(C,'selected_energies',return_value=selected), \
         patch.object(C,'read_wgfa',return_value=rad), \
         patch.object(C,'read_scups',return_value=coll):
        return C.cooling_low_density_limit('n',4,temperatures)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step',choices=['curves','nlevel','native'])
    parser.add_argument('--work-dir',required=True,type=Path)
    args=parser.parse_args(); work=args.work_dir.resolve()
    work.mkdir(parents=True,exist_ok=True)
    result={'passed':True,'step':args.step}
    if args.step=='curves':
        temperatures=np.logspace(3,5,801)
        full=C.cooling_low_density_limit('n',4,temperatures)
        reduced={n:truncated_cooling(n,temperatures) for n in (29,87)}
        assert np.all(full>0) and np.all(np.isfinite(full))
        differences={str(n):float(np.max(np.abs(lam/full-1))) for n,lam in reduced.items()}
        assert differences['29']<.002 and differences['87']<.0002,differences
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,1,figsize=(6,6),sharex=True,constrained_layout=True)
        axes[0].loglog(temperatures,full,label='136 states')
        for n,lam in reduced.items():
            axes[0].loglog(temperatures,lam,'--',label=str(n)+' states')
            axes[1].semilogx(temperatures,100*(lam/full-1),label=str(n)+' states')
        axes[0].set_ylabel(r'Cooling (erg cm$^3$ s$^{-1}$)')
        axes[1].set_ylabel(r'Difference (\%)');axes[1].set_xlabel('Electron temperature (K)')
        for ax in axes: ax.legend();ax.grid(alpha=.25)
        fig.savefig(work/'niv_level_cooling.pdf');plt.close(fig)
        np.savetxt(work/'level_curves.txt',np.column_stack([temperatures,full,reduced[29],reduced[87]]),
                   header='T_K Lambda136 Lambda29 Lambda87 [erg cm3 s-1]')
        result.update(points=801,max_relative_differences=differences,
                      comparison='Same source, energies, initial ion and stationary limit; cutoff only changes',
                      limitation='Source resonance uncertainty and complete target convergence are not certified')
    elif args.step=='nlevel':
        path,nl,nr,nu,error=N.build_ion('NIV','n',4,29,output_dir=work)
        assert (nl,nr,nu)==(29,170,406)
        assert error<N.ERR_TARGET,error
        result.update(levels=nl,radiation=nr,collisions=nu,max_transition_fit_error=error)
    else:
        exe=compile_physics_check('tests/ion_stages/check_niv_staging.f90',work)
        output=work/'native_rates.txt'
        subprocess.run([str(exe),str(work),str(output)],cwd=ROOT,check=True)
        rows=np.loadtxt(output);assert rows.shape==(2403,33)
        a,t=coefficients(work/'cooling_n_4.txt')
        fit_error=float(np.max(np.abs(rows[:,2]/evaluate(a,t,rows[:,0])-1)))
        assert fit_error<1e-10,fit_error
        low=rows[::3];source=truncated_cooling(29,low[:,0])
        low_error=float(np.max(np.abs(low[:,3]/source-1)))
        assert low_error<.005,low_error
        # Independent finite-density matrix uses source spline strengths.
        import sys
        sys.path.insert(0,str(ROOT))
        from tests.fitting.test_cooling_stationary_limit import finite_density_cooling
        levels,rad,coll=raw_atom(); levels={i:r for i,r in levels.items() if i<=29}
        rad=[r for r in rad if r[0]<=29 and r[1]<=29]
        coll=[r for r in coll if r['ll']<=29 and r['ul']<=29]
        errors=[]
        for temperature in (1000.,10000.,100000.):
            for density in (1e-6,100.,1e6):
                row=rows[np.argmin(abs(rows[:,0]-temperature)+abs(rows[:,1]-density))]
                expected=finite_density_cooling(levels,rad,coll,temperature,density)
                errors.append(abs(row[3]/expected-1))
        assert max(errors)<.005,max(errors)
        result.update(rows=2403,temperatures=801,density_cm3=[1e-6,100,1e6],
                      max_fit_reader_error=fit_error,max_low_density_error=low_error,
                      max_independent_finite_density_error=max(errors),
                      native_population_checks='All finite, nonnegative and normalized')
    (work/(args.step+'_result.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':
    main()
