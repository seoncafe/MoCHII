"""Validate expanded active atoms in thermal and photon-energy consumers."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from run_checks import production_objects

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools/fitting'))
from mocassin_excitation_profiles import source_atom
from validate_excitation_density_profiles import full_excitation_balance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--installation', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--diagnostic-temperature-nodes', type=int, choices=(401,801,1601),
                        default=401, help='Refine only a diagnostic copy of the cooling module')
    parser.add_argument('--native-timeout',type=int,default=300)
    parser.add_argument('--reuse-verified-profiles',type=Path,
                        help='Reuse only complete matching profiles from a terminal diagnostic')
    args = parser.parse_args()
    if not 1 <= args.native_timeout <= 1800:
        raise ValueError('Require a finite native limit no greater than 1800 seconds')
    work = args.work_dir.resolve()
    if work.exists():
        raise FileExistsError('Preserve earlier validation')
    work.mkdir(parents=True)
    installed = json.loads(args.installation.read_text())
    if not installed['passed'] or installed['expanded_atoms'] != 6:
        raise ValueError('Incomplete expansion')
    driver = ROOT/'tests/ion_stages/check_shared_atom_cooling.f90'
    exe = work/'shared_cooling.x'
    hdf = Path('/data/opt/hdf5_intel/lib')
    objects = [p for p in production_objects() if p.name != 'main.o']
    sources = [args.installation, Path(__file__), driver]
    if args.diagnostic_temperature_nodes != 401:
        original_module = ROOT/'src/nlevel_cooling_mod.f90'
        text = original_module.read_text()
        count = args.diagnostic_temperature_nodes
        if text.count('NT = 401') != 1 or text.count('DT_L  = 0.01_wp') != 1:
            raise ValueError('Production grid declaration changed')
        text = text.replace('NT = 401',f'NT = {count}').replace(
            'DT_L  = 0.01_wp',f'DT_L  = {4.0/(count-1):.17g}_wp')
        diagnostic = work/'nlevel_cooling_mod.f90'
        diagnostic.write_text(text)
        obj = work/'nlevel_cooling_mod.o'
        subprocess.run(['mpiifort','-cpp','-DMPI','-DHDF5','-O0','-no-ipo',
            '-check','all,noarg_temp_created','-fpe0','-traceback','-qopenmp',
            '-Isrc','-module',str(work),'-c',str(diagnostic),'-o',str(obj)],
            cwd=ROOT,check=True)
        objects = [p for p in objects if p.name != 'nlevel_cooling_mod.o']+[obj]
        sources.extend((original_module,diagnostic,obj))
    subprocess.run(['mpiifort', '-cpp', '-DMPI', '-DHDF5', '-O0', '-no-ipo',
        '-check', 'all,noarg_temp_created', '-fpe0', '-traceback', '-qopenmp',
        '-Isrc', '-module', str(work), '-o', str(exe), str(driver), *map(str,objects),
        '-lcfitsio', str(ROOT/'SEDust/sed/lib/libsedust.a'),
        str(hdf/'libhdf5_fortran.a'), str(hdf/'libhdf5.a'),
        '-lsz', '-ldl', '-lz', '-lm', '-qopenmp'], cwd=ROOT, check=True)
    results = []
    sources.append(exe)
    reuse = args.reuse_verified_profiles.resolve() if args.reuse_verified_profiles else None
    if reuse:
        saved_state = reuse.parent/'controller/state.json'
        state = json.loads(saved_state.read_text())
        if state['status'] not in ('completed','failed'):
            raise ValueError('Do not reuse a live or uncertain diagnostic')
        for p in (driver,ROOT/'src/species_mod.f90',ROOT/'src/nlevel_mod.f90',
                  ROOT/'src/nlevel_cooling_mod.f90',ROOT/'src/define.f90',
                  ROOT/'tools/fitting/mocassin_excitation_profiles.py',
                  ROOT/'tools/fitting/validate_excitation_density_profiles.py'):
            if state['source_sha256'].get(str(p)) != hashlib.sha256(p.read_bytes()).hexdigest():
                raise ValueError('Native physics or reference source changed')
        if args.diagnostic_temperature_nodes != 401 and (reuse/'nlevel_cooling_mod.f90').read_bytes() != diagnostic.read_bytes():
            raise ValueError('Diagnostic grid source changed')
        sources.append(saved_state)
    reused_profiles = 0
    densities = [1e-6,100.,1e6,1e10]
    for record in installed['atoms']:
        atom = Path(record['selected_path'])
        if hashlib.sha256(atom.read_bytes()).hexdigest() != record['native_sha256']:
            raise ValueError('Installed model changed')
        element,stage = record['element'],record['stage']
        name = element+'_'+str(stage)
        fixture = work/(name+'_atomic')
        fixture.mkdir()
        fit = ROOT/'data/atomic'/('cooling_'+name+'.txt')
        rate = ROOT/'data/atomic'/('element_'+element+'.txt')
        for p in (atom,fit,rate):
            shutil.copyfile(p,fixture/p.name)
            sources.append(p)
        if element == 's':
            # S II levels are needed by sulfur product-term charge transfer.
            # Zero ONLY its fixture cooling fit to isolate the S IV consumer.
            dependency = ROOT/'data/atomic/nlevel_s_2.txt'
            shutil.copyfile(dependency,fixture/dependency.name)
            (fixture/'cooling_s_2.txt').write_text('# Diagnostic S II cooling disabled; charge-transfer atom retained.\n1\n0 0\n')
            sources.append(dependency)
        fit_rows = [s.split() for s in fit.read_text().splitlines()
                    if s.strip() and not s.lstrip().startswith('#')]
        coefficients = np.array(fit_rows[1:], dtype=float)
        if coefficients.shape != (int(fit_rows[0][0]),2):
            raise ValueError('Incomplete cooling fit')
        reference = None
        for mode in ('local_ne','low_density'):
            inp = work/(name+'_'+mode+'.in')
            inp.write_text('&parameters\n'+(' par%cooling_model="low_density"\n' if mode=='low_density' else '')+'/\n')
            output = work/(name+'_'+mode+'.txt')
            log = work/(name+'_'+mode+'.log')
            old_output = reuse/output.name if reuse else None
            old_log = reuse/log.name if reuse else None
            complete = (reuse and old_output.exists() and old_log.exists()
                        and 'PASS: shared-ion thermal coefficient' in old_log.read_text())
            if complete:
                if (reuse/inp.name).read_bytes() != inp.read_bytes():
                    raise ValueError('Input choices changed')
                for p in fixture.iterdir():
                    if (reuse/fixture.name/p.name).read_bytes() != p.read_bytes():
                        raise ValueError('Fixture atom, fit, or rate data changed')
                shutil.copyfile(old_output,output)
                shutil.copyfile(old_log,log)
                sources.extend((old_output,old_log))
                reused_profiles += 1
            else:
                run = subprocess.run([str(exe),str(fixture),str(inp),element,str(stage),str(output)],
                                     cwd=ROOT,env=dict(os.environ,OMP_NUM_THREADS='1'),
                                     capture_output=True,text=True,timeout=args.native_timeout)
                log.write_text(run.stdout+run.stderr)
                if run.returncode:
                    raise RuntimeError(name+' '+mode+' failed; inspect preserved log')
            rows = np.loadtxt(output)
            if rows.shape != (132,6) or not np.isfinite(rows).all():
                raise ValueError('Incomplete cooling and emissivity profile')
            temperature = rows[::4,0]
            if reference is None:
                reference = full_excitation_balance(*source_atom(record),temperature,densities)
            fitted = np.sum(coefficients[:,0,None]*np.exp(-coefficients[:,1,None]/temperature),axis=0)/np.sqrt(temperature)
            thermal_error = emitted_error = 0.
            for j,ne in enumerate(densities):
                power = reference[ne][0]
                expected = fitted if mode=='low_density' else fitted*power/reference[1e-6][0]
                thermal_error=max(thermal_error,float(np.max(abs(rows[j::4,4]/expected-1))))
                emitted_error=max(emitted_error,float(np.max(abs(rows[j::4,5]/power-1))))
            if thermal_error > .005 or emitted_error > 1e-8:
                raise ValueError('Shared-ion consumer differs from independent source matrix')
            results.append(dict(element=element,stage=stage,nlev=record['nlev'],mode=mode,
                samples=132,maximum_relative_thermal_difference=thermal_error,
                maximum_relative_emitted_energy_difference=emitted_error))
    result = dict(passed=True,expanded_atoms=6,thermal_samples=1584,
        emitted_energy_samples=1584,zero_density_checks=24,cases=results,
        production_installed=installed['production_installed'],
        diagnostic_temperature_nodes=args.diagnostic_temperature_nodes,
        reused_complete_profiles=reused_profiles,
        scope='Existing rate chains, selected atoms and existing fits in isolated fixtures; fixed ne and photo rates; no spatial plasma or atomic accuracy claim',
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: six expanded shared atoms, 1584 thermal and emitted-energy samples')


if __name__ == '__main__':
    main()
