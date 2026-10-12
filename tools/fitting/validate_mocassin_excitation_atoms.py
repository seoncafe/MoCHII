"""Validate converted v2 source tables and all level/cooling consumers."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from mocassin_v2_atom_inventory import MOCASSIN, read_atom
from validate_source_excitation_atom import compile_driver, grid_file, DENSITY
from validate_excitation_density_profiles import full_excitation_balance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate-dir',required=True,type=Path)
    p.add_argument('--work-dir',required=True,type=Path)
    a=p.parse_args();work=a.work_dir.resolve()
    if work.exists(): raise FileExistsError('Preserve previous validation')
    work.mkdir(parents=True)
    manifest=json.loads((a.candidate_dir/'manifest.json').read_text())
    exe=compile_driver(work)
    result={'passed':False,'atoms':[], 'scope':'Supplied finite-level excitation models only; no matched spatial plasma or absolute atomic accuracy claim'}
    for row in manifest['atoms']:
        name=f"{row['element']}_{row['stage']}";directory=work/name;directory.mkdir()
        source=read_atom(MOCASSIN/row['relative_file'],include_tables=True)
        temperature=np.asarray(source['temperatures_K'])
        atom=Path(row['path'])
        assert hashlib.sha256(atom.read_bytes()).hexdigest()==row['native_sha256']
        grid_file(directory/'nodes.in',1,temperature)
        subprocess.run([str(exe),str(atom),str(directory/'nodes.in'),str(directory/'nodes.txt')],check=True,timeout=120)
        measured=np.loadtxt(directory/'nodes.txt',ndmin=2)[:,1:]
        levels={r['index']:(r['statistical_weight'],r['energy_cm_inverse']) for r in source['levels']}
        collisions=[]
        for r in source['collision_tables']:
            if not any(r['upsilon']):continue
            lo,up=sorted((r['lower'],r['upper']),key=lambda k:(levels[k][1],k))
            collisions.append(dict(ll=lo,ul=up,temperature_grid=temperature,upsilon_grid=np.asarray(r['upsilon'])))
        collisions.sort(key=lambda r:(r['ll'],r['ul']))
        reference=np.asarray([r['upsilon_grid'] for r in collisions]).T
        node_error=float(np.max(np.abs(measured/reference-1)))
        if node_error>1e-12:raise ValueError(name+' collision nodes differ')
        radiation=[(r['lower'],r['upper'],0.,0.,r['A_s_inverse'])
                   for r in source['spontaneous_probabilities'] if r['A_s_inverse']>0]
        # Full temperature profiles in the numerical interval, including the
        # adopted cold continuation and each original node in this interval.
        ceiling=min(1e5,temperature[-1])
        profile=np.unique(np.concatenate((np.geomspace(1000.,ceiling,17),
                      temperature[(temperature>=1000.)&(temperature<=ceiling)])))
        grid_file(directory/'profile.in',2,profile)
        subprocess.run([str(exe),str(atom),str(directory/'profile.in'),str(directory/'profile.txt')],check=True,timeout=180)
        native=np.loadtxt(directory/'profile.txt',ndmin=2)
        independent=full_excitation_balance(levels,radiation,collisions,profile,DENSITY)
        cool_error=pop_error=0.
        for j,ne in enumerate(DENSITY):
            power,pop=independent[ne];rows=native[j::len(DENSITY)]
            if np.any(power<=0):raise ValueError(name+' has nonpositive cooling')
            cool_error=max(cool_error,float(np.max(np.abs(rows[:,2]/power-1))))
            pop_error=max(pop_error,float(np.max(np.abs(rows[:,3:]-pop))))
        if cool_error>1e-8 or pop_error>1e-8:raise ValueError(name+' excitation balance differs')
        result['atoms'].append(dict(element=row['element'],stage=row['stage'],nlev=row['nlev'],nrad=row['nrad'],nups=row['nups'],
            temperatures=len(profile),node_relative_error=node_error,cooling_relative_error=cool_error,population_absolute_error=pop_error,
            native_sha256=row['native_sha256']))
        (work/'progress.json').write_text(json.dumps(result,indent=2)+'\n')
        print(f'PASS: {name}: {row["nlev"]} levels, {len(profile)} temperatures, four densities',flush=True)
    result['passed']=True
    result['source_sha256']={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in
        (a.candidate_dir/'manifest.json',Path(__file__),exe)}
    (work/'checks_result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
