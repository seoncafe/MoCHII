"""Compare cooling interpolation at fixed density against full source matrices.

Keep the existing 401 temperature nodes. Compare linear and shape-preserving
cubic interpolation of ln Lambda in -1/T, preserving the Boltzmann limit.
No production module or atomic file is modified.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from mocassin_excitation_profiles import source_atom
from validate_excitation_density_profiles import full_excitation_balance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    selection=parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--preparation', type=Path)
    selection.add_argument('--candidate-manifest', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve interpolation assessments')
    manifest=args.preparation or args.candidate_manifest
    records = json.loads(manifest.read_text())['atoms']
    regular_nodes = 10.**np.linspace(2.,6.,401)
    # Every midpoint in the source-qualified 1000--100000 K comparison range.
    regular_query = 10.**np.linspace(3.005,4.995,200)
    densities = [1e-6,100.,1e6,1e10]
    results = []
    for record in records:
        atom = source_atom(record)
        ceiling=min(float(tr['temperature_grid'][-1]) for tr in atom[2])
        nodes=regular_nodes[regular_nodes <= ceiling]
        if nodes[-1] < min(ceiling,regular_nodes[-1]):
            nodes=np.append(nodes,min(ceiling,regular_nodes[-1]))
        query=regular_query[regular_query <= ceiling]
        if not len(query):
            raise ValueError('Source ceiling below the comparison interval')
        grid = full_excitation_balance(*atom,nodes,densities)
        direct = full_excitation_balance(*atom,query,densities)
        curves = {'linear':{},'pchip':{}}
        for ne in densities:
            power = grid[ne][0]
            positive = power > 0
            x,y = -1/nodes[positive],np.log(power[positive])
            if not np.all(((-1/query)>=x[0])&((-1/query)<=x[-1])):
                raise ValueError('Query crosses a zero endpoint')
            curves['linear'][ne] = np.exp(np.interp(-1/query,x,y))
            curves['pchip'][ne] = np.exp(PchipInterpolator(x,y,extrapolate=False)(-1/query))
        for method,values in curves.items():
            errors=[]
            for ne in densities:
                cooling_error=abs(values[ne]/direct[ne][0]-1)
                density_error=abs((values[ne]/values[1e-6])/(direct[ne][0]/direct[1e-6][0])-1)
                errors.append(dict(electron_density_cm_inverse3=ne,
                    maximum_relative_cooling_difference=float(max(cooling_error)),
                    maximum_relative_density_ratio_difference=float(max(density_error))))
            results.append(dict(element=record['element'],stage=record['stage'],method=method,
                source_ceiling_K=ceiling,query_temperatures=len(query),rows=errors))
    root=Path(__file__).resolve().parents[2]
    sources=[manifest,Path(__file__),root/'tools/fitting/mocassin_excitation_profiles.py',
             root/'tools/fitting/validate_excitation_density_profiles.py']
    sources.extend(Path(r.get('selected_path',r['path'])) for r in records)
    result=dict(production_modified=False,temperature_nodes=401,query_temperatures=200,
                densities_cm_inverse3=densities,source_models=len(records),
                scope='Fixed-density interpolation and density ratio only; no density interpolation or spatial model',
                results=results,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for method in ('linear','pchip'):
        worst=max(row['maximum_relative_density_ratio_difference'] for r in results
                  if r['method']==method for row in r['rows'])
        print(method,'maximum density-ratio relative difference:',worst)


if __name__=='__main__':
    main()
