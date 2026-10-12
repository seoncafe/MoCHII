"""Validate source-limited cooling through actual rebuilt production objects."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
SOURCE_CHECKS = HERE.parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    phase = argparse.ArgumentParser()
    phase.add_argument('phase', choices=['build', 'native'])
    phase = phase.parse_args().phase
    objects = [ROOT/'src'/name for name in ('define.o', 'nlevel_mod.o', 'nlevel_cooling_mod.o')]
    if phase == 'build':
        subprocess.run(['make', 'DEBUG=1', '-j4'], cwd=ROOT, check=True, timeout=600)
        assert objects[-1].stat().st_mtime_ns >= (ROOT/'src/nlevel_cooling_mod.f90').stat().st_mtime_ns
        assert (ROOT/'MoCHII.x').stat().st_size > 0
        result = {'passed': True, 'products': {str(p): digest(p) for p in objects+[ROOT/'MoCHII.x']}}
    else:
        directory = HERE/'native'
        directory.mkdir(exist_ok=False)
        flags = ['mpiifort', '-O0', '-no-ipo', '-no-ftz', '-check',
                 'all,noarg_temp_created', '-fpe0', '-traceback', '-I'+str(ROOT/'src'),
                 '-module', str(directory)]
        driver = directory/'check_registry.x'
        subprocess.run(flags+['-o', str(driver), str(SOURCE_CHECKS/'check_registry.f90')]+
                       list(map(str,objects)), cwd=ROOT, check=True, timeout=60)
        spec = importlib.util.spec_from_file_location('source_registry_checks', SOURCE_CHECKS/'run_checks.py')
        checks = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(checks)
        checks.BUILD = directory
        checks.RUN = directory
        checks.ATOMIC = SOURCE_CHECKS/'corrected/build/atomic'
        checks.run = lambda command: subprocess.run(command, cwd=directory, text=True,
                            capture_output=True, check=True, timeout=60)
        result = {'passed': True, 'objects': {str(p):digest(p) for p in objects},
                  'driver_sha256':digest(driver), 'gradual_underflow':'Focused tests use -no-ftz; production defaults unchanged'}
        for mode, name, stage in [('argon','ar',2),('toy','xx',1),('near','zz',1),('unlimited','yy',1)]:
            result[mode] = checks.inspect(mode,name,stage)
        for name, mode in [('ar','reject'),('xx','reject'),('zz','reject_near')]:
            result[name+'_'+mode] = checks.rejection(name,2 if name=='ar' else 1,mode,
                                                   'temperature exceeds source collision table')
        for mode in ['nanT','infT','zeroT','nanNe','infNe','zeroNe']:
            result[mode] = checks.rejection('yy',1,mode,'must be positive and finite')
        small = directory/'check_small.x'
        subprocess.run(flags+['-o',str(small),str(ROOT/'tests/ion_stages/check_small_excitation_cooling.f90')]+
                       list(map(str,objects)),cwd=ROOT,check=True,timeout=60)
        fixtures=ROOT/'md/atomic_stage_expansion/validation_2026-10-10/cooling_subnormal_candidate'
        result['subnormal_and_zero']={}
        for name in ['ratio','partial','zero_ratio','zero_partial']:
            output=directory/(name+'.txt')
            completed=subprocess.run([str(small),str(fixtures/name),'test','1',str(output)],
                                     cwd=ROOT,check=True,capture_output=True,text=True,timeout=60)
            (directory/(name+'_stdout.txt')).write_text(completed.stdout)
            rows=[line for line in output.read_text().splitlines() if not line.startswith('#')]
            assert len(rows)==48
            result['subnormal_and_zero'][name]={'rows':len(rows),'sha256':digest(output)}
    (HERE/(phase+'_result.json')).write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
