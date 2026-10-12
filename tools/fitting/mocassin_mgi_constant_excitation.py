"""Prepare the v2 single-temperature Mg I approximation for explicit review.

Constant Upsilon above 10000 K is the inspected v2 assumption, not additional
published collision data. This candidate is not installed and its hot use
requires the user's physical-model decision. Preserve the approved cold data.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from mocassin_v2_atom_inventory import MOCASSIN,read_atom
from validate_source_excitation_atom import compile_driver,grid_file,DENSITY
from validate_excitation_density_profiles import full_excitation_balance

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    work=parser.parse_args().work_dir.resolve()
    if work.exists():raise FileExistsError('Preserve earlier candidates')
    work.mkdir(parents=True)
    source=MOCASSIN/'data/mgi.dat';a=read_atom(source,include_tables=True)
    cold=ROOT/'md/atomic_stage_expansion/validation_2026-10-11/mgi_cold_source/products'
    cold_result=json.loads((cold/'result.json').read_text())
    cold_atom=cold/'nlevel_mg_1.txt'
    if not cold_result['passed'] or a['temperatures_K']!=[10000.]:
        raise ValueError('Cold result or source node changed')
    for p in (source,cold_atom):
        if hashlib.sha256(p.read_bytes()).hexdigest()!=cold_result['source_sha256'][str(p)]:
            raise ValueError('Approved cold source changed')
    # Preserve all cold-model states, A values, and source strengths. Only add
    # the explicitly synthetic hot representation node; it is not new data.
    rows=cold_atom.read_text().splitlines()
    rows=rows[next(i for i,s in enumerate(rows) if s.startswith('NLEVEL_VERSION')):]
    output=['# Mg I single-node constant-Upsilon model: hot extension pending user acceptance.',
        '# Only 10000 K is a supplied collision temperature; 1000 and 100000 K nodes are synthetic.',
        '# Model interval: 1000 <= T <= 100000 K. Atomic accuracy is not certified.',
        '# Below the source node follows the approved cold rule; hot constancy follows the v2 assumption.',
        '# Above 100000 K reject evaluation. Do not infer a physical accuracy range from this model interval.',
        '# Retained radiative rates and energies are the supplied five-level source approximation.']
    cursor=0
    while cursor<len(rows):
        line=rows[cursor]
        if line=='NODES 20':line='NODES 30'
        if ' TEMPERATURE_TABLE ' in line:
            fields=line.split()
            if fields[-1]!='2':raise ValueError('Unexpected cold table length')
            first,second=rows[cursor+1].split(),rows[cursor+2].split()
            if first[0]!='1000' or second[0]!='10000' or first[1]!=second[1]:
                raise ValueError('Cold strengths are not the expected constant pair')
            fields[-1]='3';output.extend((' '.join(fields),rows[cursor+1],rows[cursor+2],'100000 '+second[1]))
            cursor+=3
        else:
            output.append(line);cursor+=1
    native=work/'nlevel_mg_1.txt';native.write_text('\n'.join(output)+'\n')
    exe=compile_driver(work)
    temperature=np.geomspace(10000.,100000.,101)
    grid_file(work/'hot.in',2,temperature)
    hot=work/'hot_profiles.txt'
    subprocess.run([str(exe),str(native),str(work/'hot.in'),str(hot)],check=True,timeout=120)
    actual=np.loadtxt(hot)
    if actual.shape!=(404,8) or not np.isfinite(actual).all():raise ValueError('Incomplete hot profile')
    levels={r['index']:(r['statistical_weight'],r['energy_cm_inverse']) for r in a['levels']}
    radiation=[(r['lower'],r['upper'],0.,0.,r['A_s_inverse']) for r in a['spontaneous_probabilities'] if r['A_s_inverse']>0]
    collisions=[dict(ll=r['lower'],ul=r['upper'],temperature_grid=np.array([1000.,10000.,100000.]),
                    upsilon_grid=np.full(3,r['upsilon'][0])) for r in a['collision_tables']]
    independent=full_excitation_balance(levels,radiation,collisions,temperature,DENSITY)
    cooling_error=population_error=0.
    for j,ne in enumerate(DENSITY):
        power,pop=independent[ne]
        cooling_error=max(cooling_error,float(np.max(abs(actual[j::4,2]/power-1))))
        population_error=max(population_error,float(np.max(abs(actual[j::4,3:]-pop))))
    if cooling_error>1e-8 or population_error>1e-8:raise ValueError('Native hot model differs from matrix')
    previous=np.loadtxt(cold/'profile.txt')
    if previous.shape!=(404,8) or not np.allclose(previous[-4:],actual[:4],rtol=1e-11,atol=0):
        raise ValueError('Hot model does not join the preserved cold source-node result')
    joined=np.vstack((previous,actual[4:]))
    np.savetxt(work/'combined_profiles.txt',joined,fmt='%.17e')
    grid_file(work/'above_model.in',1,np.array([100000.1]))
    rejected=subprocess.run([str(exe),str(native),str(work/'above_model.in'),str(work/'above_model.txt')],
                            capture_output=True,text=True,timeout=30)
    (work/'above_model.log').write_text(rejected.stdout+rejected.stderr)
    if not rejected.returncode or 'temperature exceeds collision table' not in rejected.stdout+rejected.stderr:
        raise ValueError('Model ceiling not enforced')
    fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
    for j,ne in enumerate(DENSITY):
        row=joined[j::4];label=r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$'
        axes.flat[0].loglog(row[:,0],row[:,2],label=label)
        axes.flat[0].loglog(temperature,independent[ne][0],'k:',lw=.5)
        for k in range(5):axes.flat[j+1].loglog(row[:,0],np.maximum(row[:,3+k],1e-30))
        axes.flat[j+1].set_title(label);axes.flat[j+1].set_ylabel('Level populations')
        axes.flat[5].semilogx(temperature,actual[j::4,2]/independent[ne][0]-1,label=label)
    axes.flat[0].legend(fontsize=7);axes.flat[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
    axes.flat[5].set_ylabel('Relative hot-model matrix difference')
    for ax in axes.flat:
        ax.set_xlabel(r'$T$ [K]');ax.axvline(10000.,color='.6',ls='--',lw=.6)
    fig.suptitle('Mg I constant strength: approved cold / proposed hot approximation')
    fig.savefig(work/'constant_excitation_profiles.pdf');plt.close(fig)
    result=dict(passed=True,production_installed=False,hot_model_accepted=False,
        atomic_accuracy_verified=False,source_temperature_K=10000.,model_interval_K=[1000.,100000.],
        new_hot_samples=404,reused_cold_samples=404,combined_profile_samples=804,
        maximum_hot_cooling_relative_difference=cooling_error,
        maximum_hot_population_absolute_difference=population_error,model_ceiling_rejected=True,
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       (source,cold_atom,cold/'result.json',cold/'profile.txt',native,Path(__file__),exe,hot)})
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: 404 hot candidate samples; 404 cold samples preserved; no production installation or hot acceptance.')


if __name__=='__main__':main()
