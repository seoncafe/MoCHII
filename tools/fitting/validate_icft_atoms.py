"""Validate selected 29-state O V and Cl VI excitation atoms in staging.

Historical record (2026-10-11 staging of the O V and Cl VI atoms). The
'native' step cannot run on the current tree: its driver
tests/ion_stages/check_excitation_atom.f90 was removed, and it imports
test_cooling_stationary_limit from tests/fitting, which is not on this
directory's import path. The installed atoms are checked by
tests/ion_stages/check_excitation_atoms.py.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch
import numpy as np
import chianti_cooling as C
import fit_nlevel as N
from atomic_output import atomic_write_text
from validate_added_charge_transfer import compile_physics_check

ROOT=Path(__file__).resolve().parents[2]

def atom(el,stage,n=None):
    p=Path(C.ion_dir(el,stage))/f'{el}_{stage}'
    lev,rad,coll=C.bound_excitation_atom(el,stage)
    if n:
        lev={i:r for i,r in lev.items() if i<=n}
        rad=[r for r in rad if r[0]<=n and r[1]<=n]
        coll=[r for r in coll if r['ll']<=n and r['ul']<=n]
    return lev,rad,coll


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('step',choices=['nlevel','native'])
    p.add_argument('--work-dir',required=True,type=Path)
    args=p.parse_args();work=args.work_dir.resolve();results={}
    exe=None
    if args.step=='native':exe=compile_physics_check('tests/ion_stages/check_excitation_atom.f90',work)
    for name,el,stage in [('OV','o',5),('ClVI','cl',6)]:
        out=work/'products'/name;out.mkdir(parents=True,exist_ok=True)
        if args.step=='nlevel':
            path,nl,nr,nc,error=N.build_ion(name,el,stage,29,output_dir=out)
            if error>=N.ERR_TARGET or nl!=29 or nc!=406:raise ValueError('ICFT collision fit failed')
            results[name]=dict(levels=nl,radiative_rows=nr,collision_rows=nc,maximum_transition_fit_error=error)
        else:
            dest=out/'atom_native.txt'
            subprocess.run([str(exe),str(out/f'cooling_{el}_{stage}.txt'),str(out/f'nlevel_{el}_{stage}.txt'),str(dest)],check=True)
            rows=np.loadtxt(dest)
            if rows.shape!=(2403,33):raise ValueError('Incomplete excitation atom profile')
            lev,rad,coll=atom(el,stage,29)
            with patch.object(C,'selected_energies',return_value=lev),patch.object(C,'read_wgfa',return_value=rad),patch.object(C,'read_scups',return_value=coll):
                source=C.cooling_low_density_limit(el,stage,rows[::3,0])
            low_error=float(np.max(np.abs(rows[::3,3]/source-1)))
            sys.path.insert(0,str(ROOT/'tests/fitting'))
            from test_cooling_stationary_limit import finite_density_cooling
            errors=[]
            for tindex in (0,400,800):
                for dindex in range(3):
                    row=rows[3*tindex+dindex]
                    expected=finite_density_cooling(lev,rad,coll,row[0],row[1])
                    errors.append(abs(row[3]/expected-1))
            if low_error>.005 or max(errors)>.005:raise ValueError('ICFT excitation balance mismatch')
            results[name]=dict(points=2403,maximum_low_density_difference=low_error,maximum_independent_finite_density_difference=max(errors),densities_cm3=[1e-6,100,1e6])
    atomic_write_text(work/(args.step+'_atoms_result.json'),json.dumps(dict(passed=True,ions=results),indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
