"""Isolate existing Fe II forbidden A values in the supplied 142-level model.

Historical record (2026-10-11 Fe II forbidden-transition comparison). It
cannot run on the current tree: it reads the root data/atomic/
nlevel_fe_2.txt, removed when the metal atoms moved to data/atomic/
excitation/<model>/.

This is a controlled radiation-model diagnostic, not full Fe II atomic
completion. Both calculations use identical energies, collision curves,
electron densities, and temperatures. The added A values are estimates from
the installed 13-level model; its energy identities are checked explicitly.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from validate_source_excitation_atom import compile_driver,grid_file,DENSITY
from mocassin_v2_atom_inventory import read_atom,MOCASSIN

ROOT=Path(__file__).resolve().parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work-dir',type=Path,required=True)
    a=p.parse_args();work=a.work_dir.resolve()
    if work.exists():raise FileExistsError('Preserve previous comparisons')
    work.mkdir(parents=True)
    baseline=ROOT/'md/atomic_stage_expansion/phase1_77_mocassin_v2/candidate_atoms/nlevel_fe_2.txt'
    installed=ROOT/'data/atomic/nlevel_fe_2.txt'
    records=[line.split() for line in installed.read_text().splitlines() if line.strip() and not line.startswith('#')]
    li=next(i for i,r in enumerate(records) if r[0]=='NLEV');n=int(records[li][1]);levels=records[li+1:li+1+n]
    ri=next(i for i,r in enumerate(records) if r[0]=='NRAD');nr=int(records[ri][1]);radiation=records[ri+1:ri+1+nr]
    source=read_atom(MOCASSIN/'data/feii.dat',include_tables=True)
    for row in levels:
        i,g,e=int(row[0]),float(row[1]),float(row[2]);original=source['levels'][i-1]
        if g!=original['statistical_weight'] or abs(e-original['energy_cm_inverse'])>.01:
            raise ValueError('First 13 Fe II identities require a separate physical mapping')
    text=baseline.read_text().splitlines();where=next(i for i,r in enumerate(text) if r.startswith('NRAD '))
    old=int(text[where].split()[1]);old_rows=text[where+1:where+1+old]
    existing={tuple(map(int,line.split()[:2])) for line in old_rows}
    if any(tuple(map(int,row[:2])) in existing for row in radiation):
        raise ValueError('Radiation diagnostic must add missing pairs only')
    updated=[' '.join(row[:3]) for row in radiation]
    text[where:where+1+old]=[f'NRAD {old+nr}',*old_rows,*updated]
    enriched=work/'nlevel_fe_2_with_forbidden.txt'
    enriched.write_text('# Diagnostic adds identified installed 13-level Fe II A estimates; not full atomic acceptance.\n'+ '\n'.join(text)+'\n')
    driver=compile_driver(work);T=np.geomspace(1000.,100000.,401);grid_file(work/'grid.in',2,T)
    results={}
    for name,atom in [('allowed',baseline),('with_forbidden',enriched)]:
        output=work/(name+'.txt')
        subprocess.run([str(driver),str(atom),str(work/'grid.in'),str(output)],check=True,timeout=600)
        data=np.loadtxt(output)
        if data.shape!=(1604,145):raise ValueError('Incomplete 142-level profile')
        results[name]=data
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    measured=[]
    for j,ne in enumerate(DENSITY):
        old=results['allowed'][j::4,2];new=results['with_forbidden'][j::4,2]
        ratio=new/old
        ax=axes.flat[j];ax.loglog(T,old,label='Allowed source A');ax.loglog(T,new,label='Add 40 installed forbidden A')
        ax.set_title(r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$');ax.set_xlabel(r'$T$ [K]')
        ax.set_ylabel(r'$Lambda$ [erg cm$^3$ s$^{-1}$]');ax.legend(fontsize=8)
        measured.append({'electron_density_cm_inverse_3':ne,'ratio_at_1000_K':float(ratio[0]),
                          'ratio_at_10000_K':float(ratio[200]),'ratio_at_100000_K':float(ratio[-1]),
                          'maximum_ratio':float(ratio.max())})
    fig.suptitle('Fe II: identical 142 levels and collision data; changed radiation only')
    fig.savefig(work/'forbidden_cooling_comparison.pdf');plt.close(fig)
    result={'passed':True,'samples_each_model':1604,'added_radiative_pairs':nr,'nlev':142,
        'maximum_identity_energy_difference_cm_inverse':max(abs(float(r[2])-source['levels'][int(r[0])-1]['energy_cm_inverse']) for r in levels),
        'density_comparisons':measured,'production_installed':False,
        'limits':['Only 40 currently installed forbidden A values were added.','Additional missing radiative probabilities and absolute atomic accuracy remain unverified.','Imposed fixed-ion excitation conditions; not a spatial nebular model.'],
        'source_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (baseline,installed,enriched,driver)},
        'visual_review_completed':False}
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(measured,indent=2),flush=True)


if __name__=='__main__':main()
