"""Prepare low-density cooling fits for supplied finite-level atoms.

The fits are made to the stationary state at n_e = 1e-6 cm^-3, the reference
column of the runtime cooling tables (chianti_cooling.reference_log_cooling),
not to the strict n_e -> 0 limit (decision of 2026-10-11; see log_cooling).

Candidate products are separate from production. Signed exponential terms
are a mathematical expansion of positive total cooling, not emission rates.
Accuracy applies only within the stated numerical temperature interval.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.special import logsumexp
import chianti_cooling as C
from mocassin_excitation_profiles import source_atom


def log_cooling(levels,radiation,collisions,T):
    """Natural log of the fit target: the stationary balance at
    n_e = 1e-6 cm^-3 (C.REFERENCE_DENSITY), the reference column the runtime
    divides by (nlevel_cooling_mod::nlevel_cooling_apply). Not the strict
    n_e -> 0 limit, in which a placeholder spontaneous rate (e.g. CHIANTI
    A = 1e-30 s^-1 for a one-photon-forbidden J = 0 -> 0 decay) radiates
    every excitation into a level that is collisionally de-excited at any
    physical density. Below 1e-280 erg cm^3 s^-1 the logarithmic strict
    ground-state sum takes over (C.reference_log_cooling)."""
    return C.reference_log_cooling(levels,radiation,collisions,T)[0]


def strict_log_cooling(levels,radiation,collisions,T):
    """Natural log of the strict n_e -> 0 limit (former fit target)."""
    closed,_=C.radiative_cascade_probabilities(levels,collisions,radiation)
    ground=min(levels,key=lambda k:levels[k][1])
    if closed==[ground]:
        return C.log_ground_state_cooling_from_levels(levels,radiation,collisions,T)
    source=C.cooling_low_density_from_levels(levels,radiation,collisions,T)
    if np.any(source<=0):raise ValueError('Closed-state cold tail needs logarithmic balance')
    return np.log(source)


def exponential_fit(levels,radiation,collisions,maximum):
    temperature=np.geomspace(1000.,maximum,401)
    dense=np.geomspace(1000.,maximum,1601)
    log_source=log_cooling(levels,radiation,collisions,temperature)
    dense_source=log_cooling(levels,radiation,collisions,dense)
    closed,probability=C.radiative_cascade_probabilities(levels,collisions,radiation)
    end_energy=np.array([levels[i][1] for i in closed])
    thresholds=[]
    for tr in collisions:
        lo,up=C.excitation_levels(levels,tr)
        if lo in closed and levels[up][1]-np.dot(probability[up],end_energy)>0:
            energy=levels[up][1]-min(e for g,e in levels.values())
            if energy>0:thresholds.append(energy*C.HC_ERG_CM/C.K_B_ERG)
    if not thresholds:raise ValueError('No photon-emitting excitation from a radiative end state')
    # Excited radiative end states also have Boltzmann-suppressed populations
    # in the cold limit. Factor the radiating target's energy above the ground,
    # rather than its small gap above an already excited end state. This is
    # numerical conditioning of the same total-cooling function.
    minimum=min(thresholds)
    target=log_source+.5*np.log(temperature)+minimum/temperature
    scale=float(target.max());normalized=np.exp(target-scale)
    best=None
    for count in (12,18,24,32,48,64,96,128):
        for spacing in (.00125,.0025,.005,.01,.02,.04,.08):
            offsets=np.arange(count)*spacing*maximum
            basis=np.exp(-offsets[None,:]/temperature[:,None])
            relative=basis/normalized[:,None]
            constraints=np.vstack((
                np.column_stack((relative,-relative,-np.ones(len(temperature)))),
                np.column_stack((-relative,relative,-np.ones(len(temperature)))),
                np.column_stack((basis*(1e-5-1),basis*(1e-5+1),np.zeros(len(temperature))))))
            rhs=np.r_[np.ones(len(temperature)),-np.ones(len(temperature)),np.zeros(len(temperature))]
            solution=linprog(np.r_[np.zeros(2*count),1.],A_ub=constraints,b_ub=rhs,
                            bounds=[(0,None)]*(2*count+1),method='highs')
            if not solution.success:continue
            coefficient=(solution.x[:count]-solution.x[count:2*count])*np.exp(scale)
            excitation=minimum+offsets
            parts=coefficient[:,None]*np.exp(-(excitation[:,None]-minimum)/dense)
            total=parts.sum(axis=0)
            if np.any(total<=0):continue
            condition=float(np.max(np.abs(parts).sum(axis=0)/total))
            log_fit=np.log(total)-minimum/dense-.5*np.log(dense)
            error=float(np.max(np.abs(np.expm1(log_fit-dense_source))))
            if condition<1e5 and error<.005:
                keep=coefficient!=0
                return coefficient[keep],excitation[keep],error,condition,dense,np.exp(dense_source)
            if best is None or error<best:best=error
    raise ValueError(f'No accepted exponential fit; best relative error={best}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--ion',action='append',metavar='EL_STAGE',
                   help='Prepare only these ions (e.g. n_2); repeat for several')
    p.add_argument('--retry-manifest',type=Path,
                   help='Prepare only unaccepted fits; preserve completed earlier products')
    a=p.parse_args();out=a.output_dir.resolve()
    if out.exists():raise FileExistsError('Preserve previous fit products')
    out.mkdir(parents=True)
    manifest=json.loads((a.candidate_dir/'manifest.json').read_text())
    if a.retry_manifest:
        retry={r['ion'] for r in json.loads(a.retry_manifest.read_text())['requires_review']}
        manifest['atoms']=[r for r in manifest['atoms'] if f"{r['element']}_{r['stage']}" in retry]
    if a.ion:
        manifest['atoms']=[r for r in manifest['atoms'] if f"{r['element']}_{r['stage']}" in set(a.ion)]
    result={'passed':False,'fits':[],'requires_review':[],'production_installed':False}
    for row in manifest['atoms']:
        name=f"{row['element']}_{row['stage']}"
        try:
            levels,radiation,collisions=source_atom(row)
            maximum=min(1e5,row['temperature_K'][-1])
            amplitude,excitation,error,condition,T,source=exponential_fit(levels,radiation,collisions,maximum)
        except ValueError as error:
            result['requires_review'].append({'ion':name,'reason':str(error)})
            print(f'Review required: {name}: {error}',flush=True)
        else:
            text=[f'# Supplied MOCASSIN v2 finite-level cooling model: {name}',
                  f'# Source atom SHA-256: {row["native_sha256"]}',
                  '# Stationary electron-collision/spontaneous-radiation balance at n_e = 1e-06 cm^-3 '
                  '(the reference density of the cooling tables), not the strict n_e -> 0 limit.',
                  '# Log-log PCHIP strengths; adopted constant strengths below first source node.',
                  '# Signed coefficients describe total positive cooling, not individual radiation.',
                  f'# Numerical fit interval: 1000..{maximum:g} K; no extrapolation accuracy claim.',
                  '# Enforce this interval in native consumers before installation; this file is a candidate.',
                  f'# Maximum relative error on 1601 temperatures: {error:.17g}.',
                  '# Atomic accuracy and original publication correspondence remain unverified.',str(len(amplitude))]
            text.extend(f'{v:.17g} {t:.17g}' for v,t in zip(amplitude,excitation))
            path=out/f'cooling_{name}.txt';path.write_text('\n'.join(text)+'\n')
            np.savez_compressed(out/(name+'.npz'),temperature_K=T,source_cooling=source,
                                amplitude=amplitude,excitation_K=excitation)
            result['fits'].append({'ion':name,'terms':len(amplitude),'maximum_relative_error':error,
                                    'maximum_cancellation':condition,'temperature_maximum_K':maximum,
                                    'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
            print(f'Prepared: {name}: {len(amplitude)} terms; error={error:.6g}',flush=True)
        (out/'progress.json').write_text(json.dumps(result,indent=2)+'\n')
    result['passed']=not result['requires_review']
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
