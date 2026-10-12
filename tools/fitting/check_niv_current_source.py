"""Check the current N IV source against 29-state and full-bound native atoms."""

import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys

import numpy as np

import chianti_cooling as C
from atomic_output import atomic_write_text
from validate_excitation_density_profiles import full_excitation_balance
from validate_source_excitation_atom import DENSITY, source_subset
from write_excitation_atom import _selection_code_files, _source_files


ROOT = Path(__file__).resolve().parents[2]
MAX_BYTES = 8 * 1024**3


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def memory_limit():
    resource.setrlimit(resource.RLIMIT_AS, (MAX_BYTES, MAX_BYTES))


def full_source(directory):
    memory_limit()
    directory.mkdir(parents=True, exist_ok=True)
    temperature = np.geomspace(1e3, 1e5, 801)
    levels, radiation, collisions = C.bound_excitation_atom('n', 4)
    if len(levels) < 29:
        raise ValueError('The bound N IV source has fewer than 29 levels')
    selected, selected_radiation, selected_collisions = source_subset(
        'n', 4, sorted(levels)[:29])
    full = full_excitation_balance(levels, radiation, collisions,
                                   temperature, DENSITY)
    restricted = full_excitation_balance(selected, selected_radiation,
                                         selected_collisions, temperature, DENSITY)
    data = {'temperature_K': temperature}
    maxima = []
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True,
                                constrained_layout=True)
    difference_figure, difference_axes = plt.subplots(2, 2, figsize=(9, 7),
                                                       sharex=True, constrained_layout=True)
    for index, density in enumerate(DENSITY):
        reference, full_population = full[density]
        chosen, selected_population = restricted[density]
        if (not np.all(np.isfinite(reference)) or np.any(reference <= 0)
            or not np.all(np.isfinite(chosen)) or np.any(chosen < 0)):
            raise ValueError(f'Invalid N IV source cooling at density {density:g}')
        relative = chosen / reference - 1
        maxima.append(float(np.max(np.abs(relative))))
        data[f'cooling_density_{index}'] = reference
        data[f'cooling_29_density_{index}'] = chosen
        data[f'populations_full_density_{index}'] = full_population
        data[f'populations_29_density_{index}'] = selected_population
        label = r'$n_e = ' + f'{density:g}' + r'$ cm$^{-3}$'
        axis = axes.flat[index]
        axis.loglog(temperature, reference, label='Full bound source')
        axis.loglog(temperature, chosen, '--', label='29 selected states')
        axis.set_title(label)
        axis.set_xlabel(r'$T$ [K]')
        axis.set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
        axis.legend(fontsize=8)
        axis = difference_axes.flat[index]
        axis.semilogx(temperature, 100 * relative)
        axis.set_title(label)
        axis.set_xlabel(r'$T$ [K]')
        axis.set_ylabel(r'Difference [\%]')
    figure.savefig(directory/'cooling_profiles.pdf')
    difference_figure.savefig(directory/'cutoff_differences.pdf')
    plt.close(figure)
    plt.close(difference_figure)
    np.savez_compressed(directory/'full_source_profiles.npz', **data)
    sources = list(_source_files('n', 4)) + list(_selection_code_files('n', 4))
    sources += [Path(__file__), ROOT/'tools/fitting/chianti_cooling.py',
                ROOT/'tools/fitting/validate_excitation_density_profiles.py',
                ROOT/'tools/fitting/validate_source_excitation_atom.py']
    result = dict(passed=True, full_levels=len(levels), selected_levels=29,
                  temperature_points=len(temperature),
                  temperature_K=[float(temperature[0]), float(temperature[-1])],
                  densities_cm3=DENSITY, maximum_relative_cutoff_errors=maxima,
                  cutoff_target=0.003, selected_29_meets_cutoff=max(maxima) <= 0.003,
                  selected_native_levels=29 if max(maxima) <= 0.003 else len(levels),
                  sources={str(path): digest(path) for path in sources},
                  profiles_sha256=digest(directory/'full_source_profiles.npz'))
    atomic_write_text(directory/'result.json', json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


def native(directory, driver):
    result = json.loads((directory/'fullsource/result.json').read_text())
    if not result['passed']:
        raise ValueError('Full source calculation did not pass')
    source_profiles = directory/'fullsource/full_source_profiles.npz'
    native_dir = directory/'native'
    native_dir.mkdir(parents=True, exist_ok=True)
    levels = [29]
    if result['selected_native_levels'] != 29:
        levels.append(result['selected_native_levels'])
    for count in levels:
        work = native_dir/f'level{count}'
        work.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(ROOT/'tools/fitting/validate_source_excitation_atom.py'),
                   '--element', 'n', '--stage', '4', '--levels', str(count),
                   '--work-dir', str(work), '--points', '801',
                   '--full-source-profiles', str(source_profiles), '--driver', str(driver)]
        with (work/'validation.log').open('w') as log:
            try:
                subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=750 if count == 29 else 1050,
                               preexec_fn=memory_limit)
            except subprocess.CalledProcessError:
                if count != 29 or result['selected_29_meets_cutoff']:
                    raise
                if 'Full source cutoff target failed' not in (work/'validation.log').read_text():
                    raise
        if count != 29 and not (work/'result.json').exists():
            raise ValueError('Full-bound native validation is incomplete')
    accepted = native_dir/f"level{result['selected_native_levels']}"/'result.json'
    if not accepted.exists():
        raise ValueError('Selected native atom did not pass')
    accepted_result = json.loads(accepted.read_text())
    if not accepted_result['passed'] or not accepted_result['density_ratio_cooling']:
        raise ValueError('Selected native atom lacks measured density-ratio acceptance')
    report = dict(passed=True, selected_levels=result['selected_native_levels'],
                  selected_result=str(accepted), selected_sha256=digest(accepted),
                  driver_sha256=digest(driver), failed_29_retained=not result['selected_29_meets_cutoff'])
    atomic_write_text(native_dir/'result.json', json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


def registered(directory, driver):
    selected = json.loads((directory/'native/result.json').read_text())
    if not selected['passed']:
        raise ValueError('Native source atom is not accepted')
    source_dir = directory/'native'/f"level{selected['selected_levels']}"
    work = directory/'registered'
    command = [sys.executable, str(ROOT/'tools/fitting/validate_registered_excitation_atom.py'),
               '--source-dir', str(source_dir), '--element', 'n', '--stage', '4',
               '--driver', str(driver), '--work-dir', str(work), '--timeout', '1500']
    subprocess.run(command, cwd=ROOT, check=True, timeout=1560,
                   preexec_fn=memory_limit)
    measured = json.loads((work/'result.json').read_text())
    if not measured['passed'] or measured['native_peak_rss_KiB'] > 8*1024*1024:
        raise ValueError('Registered native table failed or exceeded 8 GiB')
    print(json.dumps(measured, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['fullsource', 'native', 'registered'])
    parser.add_argument('--work-dir', required=True, type=Path)
    parser.add_argument('--driver', type=Path)
    args = parser.parse_args()
    directory = args.work_dir.resolve()
    if args.mode == 'fullsource':
        full_source(directory/'fullsource')
    elif args.mode == 'native':
        if args.driver is None:
            parser.error('native mode requires --driver')
        native(directory, args.driver.resolve())
    else:
        if args.driver is None:
            parser.error('registered mode requires --driver')
        registered(directory, args.driver.resolve())


if __name__ == '__main__':
    main()
