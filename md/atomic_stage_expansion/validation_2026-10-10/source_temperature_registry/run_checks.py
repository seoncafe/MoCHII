"""Build and check source-limited and legacy cooling tables in isolation."""

from pathlib import Path
import hashlib
import json
import math
import shutil
import subprocess
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RUN = HERE / 'corrected'
BUILD = RUN / 'build'
ATOMIC = BUILD / 'atomic'
FLAGS = ['mpifort', '-O1', '-g', '-fcheck=all',
         '-ffpe-trap=invalid,zero,overflow', '-ffree-line-length-none']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command, *, cwd=BUILD):
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True,
                          timeout=120, check=True)


def prepare():
    RUN.mkdir(exist_ok=True)
    BUILD.mkdir(exist_ok=True)
    ATOMIC.mkdir(exist_ok=True)
    for name in ('define.f90', 'nlevel_mod.f90', 'nlevel_cooling_mod.f90'):
        shutil.copy2(ROOT / 'src' / name, BUILD / name)
        run(FLAGS + ['-c', name])
    run(FLAGS + [str(HERE / 'check_registry.f90'), 'define.o', 'nlevel_mod.o',
                 'nlevel_cooling_mod.o', '-o', 'check_registry.x'])
    ar = ROOT / 'md/atomic_stage_expansion/validation_2026-10-10/fine_structure_native/ArII/nlevel_ar_2.txt'
    shutil.copy2(ar, ATOMIC / 'nlevel_ar_2.txt')
    (ATOMIC / 'nlevel_xx_1.txt').write_text('''\
NLEVEL_VERSION 2
NLEV 2
1 1 0
2 3 1000
NRAD 1
1 2 0.01
NODES 3
NUPS 1
1 2 TEMPERATURE_TABLE 0 0 0 LOG_PCHIP_T 3
100 1
1000 2
158123 3
''')
    (ATOMIC / 'nlevel_yy_1.txt').write_text('''\
NLEV 2
1 1 0
2 3 1000
NRAD 1
1 2 0.01
NUPS 1
1 2 2 1 0.02 0 1 0 1 0.3
''')
    (ATOMIC / 'nlevel_zz_1.txt').write_text('''\
NLEVEL_VERSION 2
NLEV 2
1 1 0
2 3 1000
NRAD 1
1 2 0.01
NODES 3
NUPS 1
1 2 TEMPERATURE_TABLE 0 0 0 LOG_PCHIP_T 3
100 1
1000 2
999900 3
''')


def inspect(mode, name, stage):
    command = [str(BUILD / 'check_registry.x'), str(ATOMIC), name,
               str(stage), mode]
    start = time.monotonic()
    result = run(command)
    (RUN / f'{mode}_stdout.txt').write_text(result.stdout)
    rows = []
    for line in result.stdout.splitlines():
        fields = line.split()
        if fields and fields[0] == 'ROW':
            _, _, index, T, ne, direct, table = fields
            values = [float(T), float(ne), float(direct), float(table)]
            assert all(math.isfinite(v) for v in values)
            assert values[2] >= 0 and values[3] >= 0
            error = abs(values[3] / values[2] - 1) if values[2] > 0 else abs(values[3])
            rows.append({'index': int(index), 'T_K': values[0],
                         'ne_cm3': values[1], 'direct': values[2],
                         'table': values[3], 'relative_error': error})
    assert len(rows) == 24
    exact = [r for r in rows if r['index'] in (1, 2, 3, 6)]
    assert max(r['relative_error'] for r in exact) < 1e-9
    assert max(r['relative_error'] for r in rows) < 0.005
    return {'elapsed_seconds': time.monotonic() - start,
            'max_relative_error': max(r['relative_error'] for r in rows),
            'max_exact_node_error': max(r['relative_error'] for r in exact),
            'rows': rows}


def rejection(name, stage, mode, expected):
    command = [str(BUILD / 'check_registry.x'), str(ATOMIC), name, str(stage), mode]
    proc = subprocess.run(command, cwd=BUILD, text=True, capture_output=True,
                          timeout=120)
    (RUN / f'{name}_{mode}_stdout.txt').write_text(proc.stdout + proc.stderr)
    assert proc.returncode != 0
    assert expected in proc.stderr
    assert 'SIGFPE' not in proc.stderr
    return {'exit_code': proc.returncode, 'expected_diagnostic': True}


if __name__ == '__main__':
    prepare()
    data = {'source_sha256': {name: sha(ROOT / 'src' / name)
                              for name in ('define.f90', 'nlevel_mod.f90',
                                           'nlevel_cooling_mod.f90')},
            'executable_sha256': sha(BUILD / 'check_registry.x')}
    for mode, name, stage in [('argon', 'ar', 2), ('toy', 'xx', 1),
                              ('near', 'zz', 1), ('unlimited', 'yy', 1)]:
        data[mode] = inspect(mode, name, stage)
    data['above_ceiling'] = rejection('xx', 1, 'reject', 'temperature exceeds source collision table')
    data['argon_above_ceiling'] = rejection('ar', 2, 'reject', 'temperature exceeds source collision table')
    data['near_above_ceiling'] = rejection('zz', 1, 'reject_near', 'temperature exceeds source collision table')
    for mode, expected in [('nanT', 'temperature must be positive and finite'),
                           ('infT', 'temperature must be positive and finite'),
                           ('zeroT', 'temperature must be positive and finite'),
                           ('nanNe', 'electron density must be positive and finite'),
                           ('infNe', 'electron density must be positive and finite'),
                           ('zeroNe', 'electron density must be positive and finite')]:
        data[mode] = rejection('xx', 1, mode, expected)
    (RUN / 'result.json').write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps({k: v['max_relative_error'] for k, v in data.items()
                      if isinstance(v, dict) and 'max_relative_error' in v}))
