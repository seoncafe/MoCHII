#!/usr/bin/env python3
"""Prepare and execute the bounded continuous-energy MC/RQMC workflow."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = ROOT/'results/qmc_continuous_2026-10-18'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,data):
    path.write_text(json.dumps(data,indent=2)+'\n')


def prepare():
    OUT.mkdir(exist_ok=False)
    steps = [dict(id='build',command=[sys.executable,str(__file__),'build'],
                  timeout_seconds=900,check=dict(path=str(OUT/'build.json'),equals={'ok':True}))]
    for power in range(17,21):
        for mcseed,qseed in zip((11,22,33,44),(101,202,303,404)):
            for mode,seed in (('mc',mcseed),('qmc',qseed)):
                tag=f'{mode}_p{power}_{seed}'
                original=HERE/'v3'/f'{tag}.in'
                text=original.read_text()
                text=re.sub(r"par%out_file\s*=\s*'[^']*'",f"par%out_file = '{OUT/tag}.h5'",text)
                text=text.replace('/\n', ''' par%gas_niter = 0
 par%te_fixed = 1.0e4
 par%use_metals = .false.
 par%metal_ne = .false.
 par%add_fuv = .false.
 par%ion_dust_scatter = .false.
 par%use_sec_ion = .false.
 par%heii_diffuse = .false.
 par%hei_diffuse = .false.
 par%use_master_slave = .false.
 par%ion_peel = .false.
 par%data_dir = 'DATA_ROOT'
/
'''.replace('DATA_ROOT',str(ROOT/'data')))
                (OUT/f'{tag}.in').write_text(text)
                steps.append(dict(id=tag,command=[sys.executable,str(__file__),'run',tag],
                    timeout_seconds=300,check=dict(path=str(OUT/f'{tag}.json'),equals={'ok':True})))
    steps.append(dict(id='analysis',command=[sys.executable,str(HERE/'check_continuous.py'),str(OUT)],
        timeout_seconds=600,check=dict(path=str(OUT/'summary.json'),equals={'ok':True})))
    sources=list((ROOT/'src').glob('*.f90'))+list((ROOT/'src').glob('*.inc'))
    sources += [ROOT/'Makefile',ROOT/'tests/g0_gamma/amr_sphere_L5.fits',Path(__file__),
        HERE/'check_continuous.py',HERE/'continuous_reference.py',ROOT/'tools/python/mochii_rates.py']
    sources += list(OUT.glob('*.in'))
    save(OUT/'plan.json',dict(cwd=str(ROOT),source_files=list(map(str,sources)),
        environment={'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1',
                     'I_MPI_FABRICS':'shm','I_MPI_PIN':'1'},steps=steps))
    print(OUT/'plan.json')


def build():
    subprocess.run(['make','-j4','MoCHII.x'],cwd=ROOT,check=True,timeout=850)
    shutil.copy2(ROOT/'MoCHII.x',OUT/'MoCHII.x')
    save(OUT/'build.json',dict(ok=True,binary_sha256=digest(OUT/'MoCHII.x'),
        build_configuration=(ROOT/'src/.build_configuration').read_text(),
        source_sha256={str(p.relative_to(ROOT)):digest(p) for p in (ROOT/'src').glob('*.f90')}))


def run(tag):
    output=OUT/f'{tag}_rates.h5'
    if output.exists() or (OUT/f'{tag}.json').exists():
        raise RuntimeError('Refusing to overwrite an existing calculation')
    if digest(OUT/'MoCHII.x') != json.loads((OUT/'build.json').read_text())['binary_sha256']:
        raise RuntimeError('Executable changed')
    command=['mpirun','-np','8',str(OUT/'MoCHII.x'),str(OUT/f'{tag}.in')]
    start=time.perf_counter()
    with (OUT/f'{tag}.log').open('w') as stream:
        subprocess.run(command,cwd=OUT,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=280)
    elapsed=time.perf_counter()-start
    log=(OUT/f'{tag}.log').read_text()
    assert '>>> STOP' in log and 'continuous band total (sampler)' in log
    assert 'diffuse field:' not in log
    with h5py.File(output) as f:
        data={k:f[k]['data'][:] for k in ('Gamma_HI','LeafXYZ','LeafSize','T_e','x_HI','x_HeI','n_e')}
    assert data['Gamma_HI'].size == 32768
    assert all(np.isfinite(a).all() for a in data.values())
    assert np.all(data['Gamma_HI']>=0)
    for key,target in [('T_e',1e4),('x_HI',.1),('x_HeI',1.)]:
        assert np.allclose(data[key],target,rtol=0,atol=1e-12)
    radius=np.linalg.norm(data['LeafXYZ'],axis=0)
    assert np.allclose(data['n_e'][radius<.85],.9)
    match=re.search(r'Total Execution Time\s*:\s*([0-9.]+)',log)
    save(OUT/f'{tag}.json',dict(ok=True,command=command,elapsed_seconds=elapsed,
        program_seconds=float(match[1])*60,input_sha256=digest(OUT/f'{tag}.in'),
        output_sha256=digest(output),binary_sha256=digest(OUT/'MoCHII.x')))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','build','run'])
    parser.add_argument('tag',nargs='?');args=parser.parse_args()
    {'prepare':prepare,'build':build,'run':lambda:run(args.tag)}[args.action]()
