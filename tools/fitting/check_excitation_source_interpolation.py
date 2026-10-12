"""Focused compiled checks of version-2 source-node interpolation and parsing.

The numerical interval is 1,000--100,000 K. Synthetic strengths test the
serializer/parser and interpolation formulas; they make no physical claim about
an ion. The temperature-table fixture adopts constant Upsilon below 9,800 K.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

import chianti_cooling as C
from atomic_output import atomic_write_text
from validate_source_excitation_atom import ROOT, grid_file


TMIN, TMAX = 1.0e3, 1.0e5
RELATIVE_LIMIT = 1.0e-10
TIMEOUT_SECONDS = 60


def native_executable(directory):
    executable = directory / 'check_source_excitation_atom.x'
    command = ['mpiifort', '-O0', '-no-ipo', '-check', 'all,noarg_temp_created',
               '-fpe0', '-traceback', '-Isrc', '-module', str(directory),
               '-o', str(executable),
               str(ROOT / 'tests/ion_stages/check_source_excitation_atom.f90'),
               str(ROOT / 'src/define.o'), str(ROOT / 'src/nlevel_mod.o')]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                            timeout=TIMEOUT_SECONDS)
    (directory / 'compile.log').write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f'Native driver compilation failed: {result.returncode}')
    return executable


def atom_text(rows, version=2, declared_nodes=None):
    """Construct a two-level atom; a row is (header, [(x, y), ...])."""
    parts = []
    if version == 2:
        parts.append('NLEVEL_VERSION 2')
    parts += ['NLEV 2', '1 1 0', '2 3 10000', 'NRAD 1', '1 2 1']
    if version == 2:
        actual_nodes = sum(len(nodes) for _, nodes in rows)
        parts.append(f'NODES {actual_nodes if declared_nodes is None else declared_nodes}')
    parts.append(f'NUPS {len(rows)}')
    for header, nodes in rows:
        parts.append(header)
        parts.extend(f'{x} {y}' for x, y in nodes)
    return '\n'.join(parts) + '\n'


def collision_fixture(ttype, x, y, cups=1.1, de=1.0, interpolation='CUBIC'):
    return (f'1 2 BT_SOURCE {ttype} {cups:.17g} {de:.17g} {interpolation} {len(x)}',
            list(zip((format(v, '.17g') for v in x),
                     (format(v, '.17g') for v in y))))


def source_transition(ttype, x, y, cups=1.1, de=1.0):
    return dict(ll=1, ul=2, ttype=ttype, cups=cups, de=de,
                xs=np.asarray(x), ys=np.asarray(y))


def _node_temperatures(transition):
    x = transition['xs']
    c = transition['cups']
    de = transition['de']
    if transition['ttype'] in (1, 4):
        valid = x[(x > 0) & (x < 1)]
        et = np.exp(np.log(c) / (1 - valid)) - c
    else:
        valid = x[(x > 0) & (x < 1)]
        et = c * valid / (1 - valid)
    return et * de / C.KB_OVER_RY


def temperature_grid(transition=None, table_nodes=None):
    base = np.geomspace(TMIN, TMAX, 801)
    nodes = (np.asarray(table_nodes) if table_nodes is not None else
             _node_temperatures(transition))
    neighborhoods = np.ravel(nodes[:, None] * np.array([1 - 1e-9, 1, 1 + 1e-9]))
    return np.unique(np.r_[base, neighborhoods[(neighborhoods >= TMIN) &
                                                   (neighborhoods <= TMAX)],
                            9800 * (1 - 1e-9), 9800., 9800 * (1 + 1e-9)])


def type2_spline_strength(x, y, temperature, interpolation):
    """Independent type-2 expectation for an explicitly chosen spline."""
    et = C.KB_OVER_RY * temperature
    scaled = np.clip(et / (et + 1.1), x[0], x[-1])
    spline = (CubicSpline(x, y) if interpolation == 'CUBIC'
              else PchipInterpolator(x, y))
    return spline(scaled)


def execute_fixture(directory, executable, name, contents, temperature,
                    expected=None, reject=False):
    atom = directory / f'{name}.txt'
    grid = directory / f'{name}.grid'
    output = directory / f'{name}.native'
    atom.write_text(contents, encoding='ascii')
    grid_file(grid, 1, temperature)
    started = time.monotonic()
    result = subprocess.run([str(executable), str(atom), str(grid), str(output)],
                            cwd=ROOT, text=True, capture_output=True,
                            timeout=TIMEOUT_SECONDS)
    (directory / f'{name}.log').write_text(result.stdout + result.stderr)
    elapsed = time.monotonic() - started
    if reject:
        if result.returncode == 0:
            raise AssertionError(f'Malformed atom accepted: {name}')
        return dict(name=name, rejected=True, exit_code=result.returncode,
                    seconds=elapsed)
    if result.returncode:
        raise AssertionError(f'Valid atom rejected: {name}; see {name}.log')
    values = np.loadtxt(output)
    if values.shape != (len(temperature), 2):
        raise AssertionError(f'Incomplete native output: {name}')
    if not np.allclose(values[:, 0], temperature, rtol=1e-13, atol=0):
        raise AssertionError(f'Temperature grid changed: {name}')
    actual = values[:, 1]
    if not np.all(np.isfinite(actual)) or np.any(actual < 0):
        raise AssertionError(f'Invalid native strength: {name}')
    zero = expected == 0
    if np.any(actual[zero] != 0):
        raise AssertionError(f'Source zeros changed: {name}')
    positive = ~zero
    error = float(np.max(np.abs(actual[positive] / expected[positive] - 1)))
    if error > RELATIVE_LIMIT:
        raise AssertionError(f'Source interpolation differs: {name}: {error}')
    return dict(name=name, passed=True, points=len(temperature),
                exact_zero_points=int(zero.sum()), maximum_relative_error=error,
                seconds=elapsed)


def checks(directory, executable):
    results = []
    x = np.array([0., .07, .29, .62, .87, 1.])
    y = np.array([1., 1.12, 1.18, 1.15, 1.21, 1.25])
    irregular_x = np.array([0., .035, .21, .235, .61, .9, 1.])
    irregular_y = np.array([1., 1.13, 1.1, 1.16, 1.2, 1.17, 1.22])
    for ttype in range(1, 7):
        xx, yy = (irregular_x, irregular_y) if ttype in (2, 3, 6) else (x, y)
        if ttype == 6:
            yy = np.log10(yy) - 1.2  # Negative logarithmic ordinates are valid.
        transition = source_transition(ttype, xx, yy)
        if isinstance(C.scaled_collision_strength_interpolator(transition), PchipInterpolator):
            raise AssertionError(f'Expected the source cubic for BT type {ttype}')
        grid = temperature_grid(transition)
        expected = C.upsilon(transition, grid)
        fixture = collision_fixture(ttype, xx, yy)
        results.append(execute_fixture(directory, executable, f'bt_type_{ttype}',
                                       atom_text([fixture]), grid, expected))

    # Three exact-zero source nodes followed by positive strength force PCHIP.
    zero_x = np.array([0., .025, .08, .2, .6, 1.])
    zero_y = np.array([0., 0., 0., .5, .7, .75])
    zero_transition = source_transition(2, zero_x, zero_y)
    if not isinstance(C.scaled_collision_strength_interpolator(zero_transition), PchipInterpolator):
        raise AssertionError('Source zero interval did not select PCHIP')
    grid = temperature_grid(zero_transition)
    expected = C.upsilon(zero_transition, grid)
    results.append(execute_fixture(directory, executable, 'pchip_zero',
                                   atom_text([collision_fixture(2, zero_x, zero_y,
                                                                interpolation='PCHIP')]),
                                   grid, expected))

    # Two and three nodes use distinct native endpoint formulas. Unequal
    # spacing in the three-node fixtures checks the parabolic cubic branch.
    for node_count, node_x, node_y in (
        (2, np.array([.015, .61]), np.array([1.3, 1.9])),
        (3, np.array([.015, .12, .61]), np.array([1.3, 1.9, 1.25])),
    ):
        transition = source_transition(2, node_x, node_y)
        grid = temperature_grid(transition)
        for interpolation in ('CUBIC', 'PCHIP'):
            name = f'{interpolation.lower()}_{node_count}_nodes'
            expected = type2_spline_strength(node_x, node_y, grid, interpolation)
            fixture = collision_fixture(2, node_x, node_y,
                                        interpolation=interpolation)
            results.append(execute_fixture(directory, executable, name,
                                           atom_text([fixture]), grid, expected))

    table_t = np.array([9800., 12000., 17500., 31000., 62000., 100000., 9.8e7])
    table_y = np.array([.14, .17, .11, .25, .23, .31, .38])
    table = dict(temperature_grid=table_t, upsilon_grid=table_y)
    grid = temperature_grid(table_nodes=table_t)
    expected = C.upsilon(table, grid)
    table_fixture = ('1 2 TEMPERATURE_TABLE 0 0 0 LOG_PCHIP_T 7',
                     list(zip((format(v, '.17g') for v in table_t),
                              (format(v, '.17g') for v in table_y))))
    result = execute_fixture(directory, executable, 'log_pchip_temperature',
                             atom_text([table_fixture]), grid, expected)
    cold = grid < table_t[0]
    native = np.loadtxt(directory / 'log_pchip_temperature.native')[:, 1]
    if not np.all(native[cold] == table_y[0]):
        raise AssertionError('Cold temperature-table continuation changed')
    result['constant_cold_points'] = int(cold.sum())
    results.append(result)

    # Independent Chebyshev T0/T1/T2 formula and type-4 descaling.
    legacy_grid = np.geomspace(TMIN, TMAX, 801)
    et = C.KB_OVER_RY * legacy_grid / .4
    st = 1 - np.log(1.1) / np.log(et + 1.1)
    u = 2 * np.clip(st, 0, 1) - 1
    expected = (1.2 + .1 * u - .03 * (2 * u * u - 1)) * np.log(et + 1.1)
    legacy = '\n'.join(['NLEV 2', '1 1 0', '2 3 10000', 'NRAD 1',
                        '1 2 1', 'NUPS 1',
                        '1 2 4 1.1 0.4 0 1 0 3 1.2 0.1 -0.03']) + '\n'
    results.append(execute_fixture(directory, executable, 'legacy_chebyshev',
                                   legacy, legacy_grid, expected))

    base = collision_fixture(2, x, y)
    malformed = {
        'unordered_nodes': atom_text([(base[0], [base[1][0], base[1][2],
                                                      base[1][1], *base[1][3:]])]),
        'duplicate_nodes': atom_text([(base[0], [base[1][0], base[1][1],
                                                      base[1][1], *base[1][3:]])]),
        'nodes_count_short': atom_text([base], declared_nodes=len(x)-1),
        'nodes_count_long': atom_text([base], declared_nodes=len(x)+1),
        'duplicate_collision_pair': atom_text([base, base]),
        'nan_node': atom_text([(base[0], [*base[1][:-1], (base[1][-1][0], 'NaN')])]),
        'negative_ordinate': atom_text([(base[0], [*base[1][:-1],
                                                     (base[1][-1][0], '-0.1')])]),
        'truncated_source': '\n'.join(atom_text([base]).splitlines()[:-1]) + '\n',
    }
    cubic_x = np.array([0., .5, 1.])
    cubic_y = np.array([0., 0., 1.])
    cubic = CubicSpline(cubic_x, cubic_y)
    if not np.any(cubic(np.linspace(0, .5, 21)) < 0):
        raise AssertionError('Malformed cubic fixture has no negative extremum')
    malformed['negative_cubic_extremum'] = atom_text([
        collision_fixture(2, cubic_x, cubic_y, interpolation='CUBIC')])
    for name, content in malformed.items():
        results.append(execute_fixture(directory, executable, name, content,
                                       np.array([1e4]), reject=True))
    results.append(execute_fixture(directory, executable,
                                   'temperature_above_table', atom_text([table_fixture]),
                                   np.array([9.8e7 * (1 + 1e-6)]), reject=True))
    for name, temperature in (('zero_temperature', 0.),
                              ('nonfinite_temperature', np.nan)):
        results.append(execute_fixture(directory, executable, name,
                                       atom_text([base]), np.array([temperature]),
                                       reject=True))
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path,
                        default=ROOT / 'md/atomic_stage_expansion/validation_2026-10-10/source_interpolation_checks')
    args = parser.parse_args(argv)
    directory = args.work_dir.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / 'result.json').exists():
        raise FileExistsError('Preserve prior result; select a new work directory')
    executable = native_executable(directory)
    rows = checks(directory, executable)
    valid = [row for row in rows if row.get('passed')]
    rejected = [row for row in rows if row.get('rejected')]
    result = {'passed': True, 'valid_fixtures': len(valid),
              'rejected_fixtures': len(rejected),
              'evaluated_temperatures': sum(row['points'] for row in valid),
              'exact_zero_points': sum(row['exact_zero_points'] for row in valid),
              'maximum_relative_error': max(row['maximum_relative_error'] for row in valid),
              'native_object_sha256': sha256((ROOT / 'src/nlevel_mod.o').read_bytes()).hexdigest(),
              'driver_sha256': sha256((ROOT / 'tests/ion_stages/check_source_excitation_atom.f90').read_bytes()).hexdigest(),
              'fixtures': rows}
    atomic_write_text(directory / 'result.json', json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
