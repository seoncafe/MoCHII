"""Measure full bound-level cooling convergence within existing native capacity.

Source audit only: diagnostic completion is distinct from model acceptance.
No photoexcitation, trapping, autoionization or excited-state ionization.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import chianti_cooling as C
from atomic_output import atomic_write_text
from validate_excitation_density_profiles import full_excitation_balance

IONS={'CV':('c',5),'NV':('n',5),'NVI':('n',6),'NeIV':('ne',4),
      'ArVI':('ar',6),'ArVII':('ar',7),'FeV':('fe',5),'FeVI':('fe',6),'CaVI':('ca',6),
      'OIV':('o',4),'NeV':('ne',5),'SVI':('s',6),'SVII':('s',7),
      'ArV':('ar',5),'MgIV':('mg',4),'FeVII':('fe',7),'FeIV':('fe',4),'SI':('s',1)}


def relative_difference(value,reference):
    normal=reference>64*np.finfo(float).tiny
    if np.any(~np.isfinite(value)) or np.any(value<0):raise ValueError('Invalid cooling')
    return float(np.max(np.abs(value[normal]/reference[normal]-1))) if normal.any() else 0.


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--ion',choices=IONS,required=True)
    parser.add_argument('--capacity-study',action='store_true',
        help='Measure larger source atoms for a storage proposal; does not change native limits')
    parser.add_argument('--work-dir',type=Path,required=True);args=parser.parse_args()
    element,stage=IONS[args.ion];out=args.work_dir.resolve()/args.ion;out.mkdir(parents=True,exist_ok=True)
    levels,rad,coll=C.bound_excitation_atom(element,stage)
    for l,u,w,g,a in rad:
        if a>0 and levels[u][1]<=levels[l][1]:raise ValueError('Source radiation direction inconsistent')
    temperature=np.geomspace(1e3,1e5,801);densities=[1e-6,100.,1e6,1e10]
    full=full_excitation_balance(levels,rad,coll,temperature,densities)
    records=[];saved={'temperature_K':temperature}
    for j,ne in enumerate(densities):
        cooling,pop=full[ne]
        saved[f'cooling_density_{j}']=cooling;saved[f'population_density_{j}']=pop
    candidate_ids=sorted(levels);selected=None;profiles={}
    candidates=[5,10,15,20,25,29,30,32,34,35,40]
    if args.capacity_study:candidates += [45,50,60,70,80,96,110,130,150,170,190,len(levels)]
    sizes=sorted(set(min(n,len(levels)) for n in candidates))
    for n in sizes:
        retained=set(candidate_ids[:n]);lev={i:r for i,r in levels.items() if i in retained}
        rr=[r for r in rad if r[0] in retained and r[1] in retained]
        cc=[tr for tr in coll if tr['ll'] in retained and tr['ul'] in retained]
        within_capacity=n<=40 and sum(r[4]>0 for r in rr)<=600 and len(cc)<=600
        if not args.capacity_study and not within_capacity:continue
        measured=full_excitation_balance(lev,rr,cc,temperature,densities)
        differences=[relative_difference(measured[ne][0],full[ne][0]) for ne in densities]
        records.append(dict(levels=n,radiative_rows=sum(r[4]>0 for r in rr),collision_rows=len(cc),within_native_capacity=within_capacity,
                            maximum_relative_cooling_differences=differences))
        profiles[n]=measured
        if selected is None and max(differences)<.003:selected=n
    if not records:raise ValueError('No candidate within native capacity')
    displayed=selected or records[-1]['levels']
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(9,7),sharex=True,constrained_layout=True)
    for j,ne in enumerate(densities):
        ax=axes.flat[j];source=full[ne][0];ax.loglog(temperature,source,label=f'Full: {len(levels)} states')
        ax.loglog(temperature,profiles[displayed][ne][0],'--',label=f'Selected: {displayed} states')
        ax.set_title(r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$');ax.set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
        ax.set_xlabel(r'$T$ [K]');ax.legend(fontsize=8)
    fig.suptitle(args.ion+' bound-level cooling');fig.savefig(out/'cooling_profiles.pdf');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(9,7),sharex=True,constrained_layout=True)
    for j,ne in enumerate(densities):
        ax=axes.flat[j];reference=full[ne][0];normal=reference>64*np.finfo(float).tiny
        for record in records:
            n=record['levels'];ax.semilogx(temperature[normal],100*(profiles[n][ne][0][normal]/reference[normal]-1),label=str(n))
        ax.set_title(r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$');ax.set_ylabel(r'Difference [\%]');ax.set_xlabel(r'$T$ [K]')
        ax.legend(title='States',fontsize=7,ncol=3)
    fig.suptitle(args.ion+' cutoff convergence');fig.savefig(out/'cutoff_differences.pdf');plt.close(fig)
    np.savez_compressed(out/'full_profiles.npz',**saved)
    raw_levels=C.read_excitation_atom(element,stage)[0]
    stem=Path(C.ion_dir(element,stage))/f'{element}_{stage}'
    source_paths=[Path(str(stem)+suffix) for suffix in ['.elvlc','.wgfa','.scups']]
    selected_record=next((r for r in records if r['levels']==selected),None)
    result=dict(passed=True,model_accepted=selected_record is not None and selected_record['within_native_capacity'],selected_levels=selected,
        capacity_study=args.capacity_study,source_convergence_accepted=selected is not None,
        full_bound_collision_connected_levels=len(levels),raw_levels=len(raw_levels),
        temperature_samples=801,temperature_K=[1000,100000],densities_cm3=densities,
        candidates=records,source_cutoff_relative_target=.003,
        native_capacity=dict(levels=40,radiative_rows=600,collision_rows=600),
        source_sha256={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths},
        comparison_domain='Relative source-cutoff comparison where full cooling exceeds 64 times minimum normal float; colder tails require logarithmic coronal diagnostics before final acceptance',
        source_accuracy='Bound states only; adopted energy vector; scattering resonance and omitted-process accuracy are not certified',
        acceptance='passed means diagnostic completed; model_accepted separately records capacity-constrained cooling convergence')
    atomic_write_text(out/'result.json',json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
