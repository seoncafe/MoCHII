"""Generate the 64-stage roster and check adopted high-charge capture physics.

Historical record (2026-10-11 64-stage roster; the default is now stages
I-VII, make_element_data.ELEMENTS). It cannot run on the current tree: it
copies the root data/atomic/nlevel_s_2.txt, removed when the metal atoms
moved to data/atomic/excitation/<model>/.

This validates atomic rate and energy consumers, not a transported nebular
model. Line cooling and nonadjacent Auger physics remain separate work.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import make_element_data as G
from charge_exchange_capture import langevin_recombination_rates, avoided_crossing_mean_heat
from validate_added_charge_transfer import compile_physics_check

ROOT=Path(__file__).resolve().parents[2]
TARGET={'c':5,'n':6,'o':5,'ne':5,'s':7,'ar':7,'mg':4,'fe':7,'si':5,'cl':7,'ca':6}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step',choices=['prepare','native'])
    parser.add_argument('--work-dir',type=Path,required=True)
    args=parser.parse_args();work=args.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
    atomic=work/'atomic';rates=langevin_recombination_rates()
    production={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (ROOT/'data/atomic').glob('element_*.txt')}
    if args.step=='prepare':
        if atomic.exists():raise RuntimeError('Preserve existing products; use a new directory')
        roster=[(el,z,TARGET[el]) for el,z,n in G.ELEMENTS]
        assert sum(n for el,z,n in roster)==64
        G.generate(roster,atomic)
        shutil.copy2(ROOT/'data/atomic/nlevel_s_2.txt',atomic/'nlevel_s_2.txt')
        result={'passed':True,'stages':64,'transitions':53,'selected_high_charge_reactions':len(rates),
                'rate_model':'H(1s) polarization capture with unit electron-transfer probability and reduced mass',
                'mean_heat_model':'KF99 avoided crossing with exit charge: 2.86*(q-1) eV for initial charge q>4',
                'temperature_accuracy_interval':None,'production_activation':False}
    else:
        exe=compile_physics_check('tests/ion_stages/check_added_charge_transfer_heat.f90',work)
        subprocess.run([str(exe),str(atomic),str(work/'thermal_rates.txt'),'64'],cwd=ROOT,check=True)
        heat_samples=0;heat_errors=[]
        for line in (work/'thermal_rates.txt').read_text().splitlines():
            p=line.split();key=p[0],int(p[1])
            if key not in rates:continue
            temperature,ki,kr,di,dr,split,frac=map(float,p[3:])
            assert ki==di==split==0 and kr>0 and 0<frac<1
            heat_errors.append(abs(dr-avoided_crossing_mean_heat(key[1])))
            heat_samples+=1
        assert heat_samples==30 and max(heat_errors)<1e-12
        index={el:i+1 for i,(el,z,n) in enumerate(G.ELEMENTS)}
        samples=[(key,10.**(1+6*j/100)) for key in rates for j in range(101)]
        (work/'temperatures.txt').write_text(''.join(
            f'{index[el]} {q} {t:.17e}\n' for (el,q),t in samples))
        exe=compile_physics_check('tests/ion_stages/check_charge_transfer_ranges.f90',work)
        subprocess.run([str(exe),str(atomic),str(work/'temperatures.txt'),str(work/'native_rates.txt')],
                       cwd=ROOT,check=True)
        errors=[]
        lines=(work/'native_rates.txt').read_text().splitlines()
        assert len(lines)==len(samples)
        for line,(key,t) in zip(lines,samples):
            ie,it,temperature,ki,kr=line.split()
            assert (int(ie),int(it),float(temperature))==(index[key[0]],key[1],t)
            assert float(ki)==0
            errors.append(abs(float(kr)/rates[key]-1))
        assert max(errors)<1e-6
        result={'passed':True,'stages':64,'transitions':53,'rate_samples':len(errors),
                'maximum_relative_rate_error':max(errors),'heat_samples':heat_samples,
                'maximum_heat_error_eV':max(heat_errors),
                'consumer_checks':'Direct/cached fractions, electron charge, reaction sums and kinetic heat passed',
                'sampled_temperature_K':[10,1e7],
                'temperature_accuracy_interval':None,
                'temperature_sampling_note':'Checks numerical temperature independence, not physical accuracy',
                'production_activation':False}
    assert production=={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (ROOT/'data/atomic').glob('element_*.txt')}
    (work/(args.step+'_result.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':
    main()
