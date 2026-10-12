"""Stage nine KF96/KF99 reactions and check the actual linked calculations.

Historical record (2026-10-11 staging of nine KF96/KF99 reactions). Its
main() cannot run on the current tree: it copies the root data/atomic/
nlevel_s_2.txt, removed when the metal atoms moved to data/atomic/
excitation/<model>/. compile_physics_check, imported by other scripts here,
is unaffected.

No production file is replaced. Mean heat is the published low-density,
ground-reactant approximation, not a determination of product branching.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import make_element_data as G
from charge_exchange_expansion import (kf96_added_rates, KF99_ADDED_HEAT,
    kf96_recombination_total_rates, KF99_RECOMBINATION_MEAN_HEAT)

ROOT=Path(__file__).resolve().parents[2]
TARGET={'c':5,'n':4,'o':5,'ne':5,'s':5,'ar':5,'mg':4}


def compile_physics_check(source,work):
    prefix=Path(os.environ.get('HDF5_PREFIX','/data/opt/hdf5_intel'))
    sedust=ROOT/'SEDust/sed/lib'
    objects=sorted(p for p in (ROOT/'src').glob('*.o') if p.name!='main.o')
    if not objects or not (ROOT/'src/species_mod.o').exists():
        raise RuntimeError('Build the native debug objects first')
    exe=work/(Path(source).stem+'.x')
    subprocess.run(['mpiifort','-cpp','-DMPI','-DHDF5','-O0','-no-ipo','-check',
        'all,noarg_temp_created','-fpe0','-traceback','-qopenmp','-Isrc',
        '-module',str(work),'-o',str(exe),str(ROOT/source),*map(str,objects),
        '-lcfitsio','-L/usr/local/lib',str(sedust/'libsedust.a'),
        str(prefix/'lib/libhdf5_fortran.a'),str(prefix/'lib/libhdf5.a'),
        '-lsz','-ldl','-lz','-lm','-qopenmp'],cwd=ROOT,check=True)
    return exe


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step',choices=['prepare','fortran','plot'])
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument("--include-reviewed-totals",action="store_true",
                        help="Include original N4 and existing Cl/Ca total rates; Fe DR remains blocked")
    args=parser.parse_args()
    target=TARGET.copy()
    reactions=kf96_added_rates()
    heats=KF99_ADDED_HEAT.copy()
    if args.include_reviewed_totals:
        target["n"]=5
        reactions.update({key:value for key,value in kf96_recombination_total_rates().items()
                          if key[0]!="fe"})
        heats.update({key:value for key,value in KF99_RECOMBINATION_MEAN_HEAT.items()
                      if key[0]!="fe"})
    expected_stages=sum(target.get(el,n) for el,z,n in G.ELEMENTS)
    work=args.work_dir.resolve(); work.mkdir(parents=True,exist_ok=True)
    atomic=work/'atomic'
    result={'passed':False,'step':args.step}
    production={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (ROOT/'data/atomic').glob('element_*.txt')}
    if args.step=='prepare':
        if atomic.exists(): raise RuntimeError('Preserve prior staging; choose a new work directory')
        roster=[(el,z,target.get(el,n)) for el,z,n in G.ELEMENTS]
        G.generate(roster,atomic)
        shutil.copy2(ROOT/'data/atomic/nlevel_s_2.txt',atomic/'nlevel_s_2.txt')
        result.update(stages=sum(n for el,z,n in roster),transitions=sum(n-1 for el,z,n in roster),
                      production_sha256=production,scope='atomic-rate tests only; missing line cooling')
    elif args.step=='fortran':
        if production!=json.loads((work/'prepare_result.json').read_text())['production_sha256']:
            raise RuntimeError('Production element data changed during staging')
        exe=compile_physics_check('tests/ion_stages/check_added_charge_transfer_heat.f90',work)
        subprocess.run([str(exe),str(atomic),str(work/'thermal_closure.txt'),str(expected_stages)],cwd=ROOT,check=True)
        seen=set(); heat_error=0.
        for line in (work/'thermal_closure.txt').read_text().splitlines():
            fields=line.split(); key=fields[0],int(fields[1])
            if key not in heats: continue
            temperature,ki,kr,di,dr,split,fraction=map(float,fields[3:])
            heat_error=max(heat_error,abs(dr-heats[key]))
            if ki!=0 or di!=0 or split!=0 or not 0<fraction<=1 or kr<=0:
                raise RuntimeError('Invalid selected staged reaction: '+line)
            seen.add((key,int(fields[2])))
        if len(seen)!=3*len(heats) or heat_error>1e-12: raise RuntimeError('Incomplete mean-heat checks')
        index={el:i+1 for i,(el,z,n) in enumerate(G.ELEMENTS)}
        samples=[]
        for key,reaction in reactions.items():
            ts=[reaction.min_temperature/10]+[reaction.min_temperature*(
                reaction.max_temperature/reaction.min_temperature)**(j/2000)
                for j in range(2001)]+[reaction.max_temperature*10]
            for temperature in ts:
                samples.append((key,temperature))
        (work/'temperatures.txt').write_text(''.join(
            f'{index[el]} {q} {t:.17e}\n' for (el,q),t in samples))
        exe=compile_physics_check('tests/ion_stages/check_charge_transfer_ranges.f90',work)
        subprocess.run([str(exe),str(atomic),str(work/'temperatures.txt'),
                        str(work/'compiled_rates.txt')],cwd=ROOT,check=True)
        output=(work/'compiled_rates.txt').read_text().splitlines()
        if len(output)!=len(samples): raise RuntimeError('Incomplete compiled rate curves')
        errors={}; curves={}
        for line,(key,t) in zip(output,samples):
            ie,it,temperature,ki,kr=line.split(); kr=float(kr)
            if (int(ie),int(it),float(temperature))!=(index[key[0]],key[1],t) or float(ki)!=0:
                raise RuntimeError('Changed reaction sample identity')
            reaction=reactions[key]
            t4=min(max(t,reaction.min_temperature),reaction.max_temperature)/1e4
            a,b,c,d=reaction.coefficients
            reference=a*1e-9*t4**b*(1+c*math.exp(d*t4))
            if not math.isfinite(kr) or kr<=0: raise RuntimeError('Invalid compiled rate')
            label=key[0]+str(key[1]); error=abs(kr/reference-1)
            errors[label]=max(errors.get(label,0),error)
            curve=curves.setdefault(label,{'temperature_K':[],'compiled_cm3_s':[],
                     'published_cm3_s':[],'mean_heat_eV':heats[key],
                     'min_temperature_K':reaction.min_temperature,'max_temperature_K':reaction.max_temperature})
            curve['temperature_K'].append(t);curve['compiled_cm3_s'].append(kr)
            curve['published_cm3_s'].append(reference)
        if max(errors.values())>1e-5: raise RuntimeError('Added rate interpolation accuracy failed')
        (work/'curves.json').write_text(json.dumps(curves,indent=2)+'\n')
        result.update(reaction_samples=len(samples),max_relative_errors=errors,
                      heat_samples=len(seen),max_heat_difference_eV=heat_error,
                      production_preserved=True)
    else:
        import matplotlib.pyplot as plt
        from matplotlib.ticker import NullFormatter
        curves=json.loads((work/'curves.json').read_text())
        ncol=3
        nrow=(len(curves)+ncol-1)//ncol
        fig,axes=plt.subplots(nrow,ncol,figsize=(10,2.8*nrow),constrained_layout=True)
        for ax,(label,curve) in zip(axes.flat,curves.items()):
            ts=curve['temperature_K'][1:-1]
            scale=curve['mean_heat_eV']*1.602176634e-12
            ax.loglog(ts,[v*scale for v in curve['compiled_cm3_s'][1:-1]],
                      label='Compiled',linewidth=1.4)
            ax.loglog(ts,[v*scale for v in curve['published_cm3_s'][1:-1]],
                      label='KF96 rate x KF99 mean heat',linestyle='--',linewidth=.8)
            ax.set_title(label+'; mean heat '+f"{curve['mean_heat_eV']:.2f}"+' eV')
            ax.set_xlabel('Temperature (K)')
            ax.set_ylabel(r'$k\Delta E$ (erg cm$^3$ s$^{-1}$)')
            # Dense minor labels overlap for the Ne3 interval (one decade).
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.grid(alpha=.25)
        for ax in list(axes.flat)[len(curves):]: ax.set_visible(False)
        axes.flat[0].legend(fontsize=6)
        fig.suptitle('Charge-transfer thermal coefficients: selected temperature intervals\n'
                     'Ground reactants; low-density mean-heat approximation')
        output=work/'added_charge_transfer_heat.pdf'
        fig.savefig(output);plt.close(fig)
        if not output.exists() or output.stat().st_size<1000: raise RuntimeError('Missing comparison figure')
        result['figure']=str(output)
    current={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
             for p in (ROOT/'data/atomic').glob('element_*.txt')}
    if current!=production: raise RuntimeError('Staging altered production data')
    result['passed']=True
    (work/(args.step+'_result.json')).write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__': main()
