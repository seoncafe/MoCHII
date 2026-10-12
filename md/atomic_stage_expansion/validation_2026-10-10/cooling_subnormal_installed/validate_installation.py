"""Build and check the installed absolute excitation-cooling correction."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['build', 'native'])
    phase = parser.parse_args().phase
    if phase == 'build':
        command = ['make', 'DEBUG=1', '-j4']
        subprocess.run(command, cwd=ROOT, check=True, timeout=600)
        products = [ROOT/'MoCHII.x', ROOT/'src/nlevel_mod.o',
                    ROOT/'src/nlevel_cooling_mod.o', ROOT/'src/cooling_mod.o',
                    ROOT/'src/species_mod.o']
        for path in products:
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f'Missing build product: {path}')
        if (ROOT/'src/nlevel_cooling_mod.o').stat().st_mtime_ns < (
                ROOT/'src/nlevel_cooling_mod.f90').stat().st_mtime_ns:
            raise ValueError('The installed cooling source was not rebuilt')
        result = {'passed': True, 'command': command,
                  'products': {str(path): digest(path) for path in products}}
    else:
        directory = HERE/'native'
        directory.mkdir(exist_ok=False)
        flags = ['mpiifort', '-O0', '-no-ipo', '-no-ftz', '-check',
                 'all,noarg_temp_created', '-fpe0', '-traceback',
                 '-Isrc', '-module', str(directory)]
        command = flags + ['-o', str(directory/'check_small.x'),
                           str(ROOT/'tests/ion_stages/check_small_excitation_cooling.f90'),
                           str(ROOT/'src/define.o'), str(ROOT/'src/nlevel_mod.o'),
                           str(ROOT/'src/nlevel_cooling_mod.o')]
        subprocess.run(command, cwd=ROOT, check=True, timeout=60)
        fixtures = ROOT/'md/atomic_stage_expansion/validation_2026-10-10/cooling_subnormal_candidate'
        records = {}
        for name in ['ratio', 'partial', 'zero_ratio', 'zero_partial']:
            output = directory/(name+'.txt')
            run = [str(directory/'check_small.x'), str(fixtures/name),
                   'test', '1', str(output)]
            completed = subprocess.run(run, cwd=ROOT, check=True, timeout=60,
                                       capture_output=True, text=True)
            (directory/(name+'_stdout.txt')).write_text(completed.stdout)
            rows = [line.split() for line in output.read_text().splitlines()
                    if not line.startswith('#')]
            if len(rows) != 48:
                raise ValueError('Incomplete subnormal/zero fixture')
            records[name] = {'rows': len(rows), 'command': run,
                             'output_sha256': digest(output)}
        result = {'passed': True, 'compile_command': command, 'fixtures': records,
                  'gradual_underflow': 'Enabled explicitly with -no-ftz for these tests',
                  'scope': 'Actual rebuilt production objects; public cooling entry points',
                  'objects': {str(path): digest(path) for path in
                              [ROOT/'src/nlevel_mod.o', ROOT/'src/nlevel_cooling_mod.o']},
                  'executable_sha256': digest(directory/'check_small.x')}
    (HERE/(phase+'_result.json')).write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
