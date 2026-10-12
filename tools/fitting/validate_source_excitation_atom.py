"""Compare a compiled source-node atom with independent excitation balance."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import chianti_cooling as C
from atomic_output import atomic_write_text
from write_excitation_atom import write_excitation_atom
from validate_excitation_density_profiles import full_excitation_balance

ROOT = Path(__file__).resolve().parents[2]
DENSITY = [1e-6, 100., 1e6, 1e10]


def compile_driver(directory):
    executable = directory/'check_source_excitation_atom.x'
    subprocess.run(['mpiifort', '-O0', '-no-ipo', '-check', 'all,noarg_temp_created',
                    '-fpe0', '-traceback', '-Isrc', '-module', str(directory),
                    '-o', str(executable),
                    str(ROOT/'tests/ion_stages/check_source_excitation_atom.f90'),
                    str(ROOT/'src/define.o'), str(ROOT/'src/nlevel_mod.o')],
                   cwd=ROOT, check=True)
    return executable


def grid_file(path, mode, temperature):
    path.write_text(f'{mode} {len(temperature)}\n'+
                    ' '.join(format(t, '.17g') for t in temperature)+'\n')


def source_subset(element, stage, ids):
    bound, radiation, collisions = C.bound_excitation_atom(element, stage)
    levels = {i: bound[i] for i in ids}
    radiation = [r for r in radiation if r[0] in levels and r[1] in levels]
    collisions = [r for r in collisions if r['ll'] in levels and r['ul'] in levels]
    index = {i: j for j, i in enumerate(ids)}
    def pair(row):
        l, u = C.excitation_levels(levels, row)
        return index[l], index[u]
    collisions.sort(key=pair)
    return levels, radiation, collisions


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--element', required=True)
    p.add_argument('--stage', type=int, required=True)
    p.add_argument('--levels', type=int, required=True)
    p.add_argument('--work-dir', type=Path, required=True)
    p.add_argument('--points', type=int, default=801)
    p.add_argument('--full-source-profiles', type=Path)
    p.add_argument('--driver', type=Path,
                   help='Previously compiled driver; preserve its isolated object files')
    p.add_argument('--temperature-min', type=float, default=1e3)
    p.add_argument('--temperature-max', type=float, default=1e5)
    args = p.parse_args()
    args.element = args.element.lower()
    if (not np.isfinite(args.temperature_min) or not np.isfinite(args.temperature_max)
        or not 0 < args.temperature_min < args.temperature_max or args.points < 2):
        raise ValueError('Require finite, positive, ordered temperatures and at least two points')
    directory = args.work_dir.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    metadata = write_excitation_atom(args.element, args.stage, args.levels, directory)
    levels, radiation, collisions = source_subset(args.element, args.stage,
                                                metadata['selected_ids'])
    executable = args.driver.resolve() if args.driver else compile_driver(directory)
    temperature = np.geomspace(args.temperature_min, args.temperature_max, args.points)
    source_nodes = np.concatenate([r['temperature_grid'] for r in collisions
                                   if 'temperature_grid' in r]) if any(
                                       'temperature_grid' in r for r in collisions) else np.array([])
    source_nodes = source_nodes[(source_nodes >= args.temperature_min)
                                & (source_nodes <= args.temperature_max)]
    strength_temperature = np.unique(np.r_[np.geomspace(args.temperature_min,
                                                       args.temperature_max, 201), source_nodes])
    grid_file(directory/'strength_grid.txt', 1, strength_temperature)
    grid_file(directory/'profile_grid.txt', 2, temperature)
    strength_output = directory/'native_strengths.txt'
    profile_output = directory/'native_profiles.txt'
    if strength_output.exists() or profile_output.exists():
        raise RuntimeError('Preserve prior products; use a new validation directory')
    subprocess.run([str(executable), metadata['path'], str(directory/'strength_grid.txt'),
                    str(strength_output)], check=True)
    measured = np.loadtxt(strength_output)
    source_strength = np.column_stack([C.upsilon(r, strength_temperature)
                                       for r in collisions])
    if measured.shape != (len(strength_temperature), len(collisions)+1):
        raise ValueError('Incomplete native strengths')
    native_strength = measured[:, 1:]
    zero = source_strength == 0
    if np.any(native_strength[zero] != 0):
        raise ValueError('Source zero interval was not preserved exactly')
    positive = ~zero
    strength_error = float(np.max(np.abs(native_strength[positive]/source_strength[positive]-1)))
    if not np.all(np.isfinite(native_strength)) or np.any(native_strength < 0) or strength_error > 1e-9:
        raise ValueError(f'Native collision source differs: {strength_error}')
    started = time.monotonic()
    subprocess.run([str(executable), metadata['path'], str(directory/'profile_grid.txt'),
                    str(profile_output)], check=True)
    native_seconds = time.monotonic()-started
    rows = np.loadtxt(profile_output)
    if rows.shape != (4*args.points, 3+args.levels):
        raise ValueError('Incomplete native populations')
    if not np.allclose(rows[:, 0], np.repeat(temperature, 4), rtol=1e-13, atol=0):
        raise ValueError('Native temperature grid differs')
    if not np.array_equal(rows[:, 1], np.tile(DENSITY, args.points)):
        raise ValueError('Native electron-density grid differs')
    direct = full_excitation_balance(levels, radiation, collisions, temperature, DENSITY)
    full = np.load(args.full_source_profiles) if args.full_source_profiles else None
    if full is not None and not np.array_equal(full['temperature_K'], temperature):
        raise ValueError('Full source temperature grid differs')
    errors, populations, cutoff = [], [], []
    saved = {'temperature_K': temperature}
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, constrained_layout=True)
    difference_fig, difference_axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True,
                                                 constrained_layout=True)
    for j, ne in enumerate(DENSITY):
        expected, pop = direct[ne]
        actual = rows[j::4, 2]
        if np.any(actual[expected == 0] != 0):
            raise ValueError('Zero source cooling was not preserved')
        normal = expected > 64*np.finfo(float).tiny
        relative = np.zeros_like(expected)
        relative[normal] = actual[normal]/expected[normal]-1
        errors.append(float(np.max(np.abs(relative[normal]))))
        populations.append(float(np.max(np.abs(rows[j::4, 3:]-pop))))
        if full is not None:
            reference = full[f'cooling_density_{j}']
            if reference.shape != expected.shape:
                raise ValueError('Full source grid differs')
            good = reference > 64*np.finfo(float).tiny
            cutoff.append(float(np.max(np.abs(actual[good]/reference[good]-1))))
        saved[f'cooling_source_{j}'] = expected
        saved[f'cooling_native_{j}'] = actual
        saved[f'populations_source_{j}'] = pop
        saved[f'populations_native_{j}'] = rows[j::4, 3:]
        title = r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$'
        ax = axes.flat[j]
        ax.loglog(temperature, expected, label='Selected source')
        ax.loglog(temperature, actual, '--', label='Compiled atom')
        if full is not None:
            ax.loglog(temperature, reference, ':', label='Full source')
        ax.set_title(title); ax.set_xlabel(r'$T$ [K]')
        ax.set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]'); ax.legend(fontsize=8)
        ax = difference_axes.flat[j]
        ax.semilogx(temperature[normal], 100*relative[normal])
        ax.set_title(title); ax.set_xlabel(r'$T$ [K]'); ax.set_ylabel(r'Difference [\%]')
    if max(errors) > 1e-8 or max(populations) > 1e-10:
        raise ValueError(f'Compiled source balance differs: {errors}, {populations}')
    if cutoff and max(cutoff) > .003:
        raise ValueError(f'Full source cutoff target failed: {cutoff}')
    if cutoff and not metadata['density_ratio_cooling']:
        path = Path(metadata['path'])
        # This flag requires measured convergence, not just structural generation.
        contents = path.read_text() + (
            '# Density ratio accepted against full source at 1000--100000 K; '
            'ne = 1e-6, 1e2, 1e6, 1e10 cm^-3; maximum relative difference < 0.003\n'
            'DENSITY_RATIO 1\n')
        atomic_write_text(path, contents)
        metadata['density_ratio_cooling'] = True
    fig.savefig(directory/'cooling_profiles.pdf')
    difference_fig.savefig(directory/'native_differences.pdf')
    plt.close(fig); plt.close(difference_fig)
    np.savez_compressed(directory/'source_native_profiles.npz', **saved)
    result = dict(metadata, passed=True, temperature_points=args.points,
                  temperature_K=[args.temperature_min, args.temperature_max],
                  densities_cm3=DENSITY, maximum_relative_strength_error=strength_error,
                  exact_source_zero_samples=int(zero.sum()),
                  maximum_relative_cooling_errors=errors,
                  maximum_absolute_population_errors=populations,
                  full_source_cutoff_errors=cutoff, native_profile_seconds=native_seconds,
                  product_sha256=hashlib.sha256(Path(metadata['path']).read_bytes()).hexdigest(),
                  native_object_sha256=(hashlib.sha256((ROOT/'src/nlevel_mod.o').read_bytes()).hexdigest()
                                        if not args.driver else None),
                  driver_path=str(executable),
                  driver_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                  scope='Source representation and statistical equilibrium; no spatial plasma model')
    atomic_write_text(directory/'result.json', json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
