"""Validate the private Mg III dipole-strength atom at full temperature/density resolution."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'tools/fitting'))
import chianti_cooling as C
from validate_excitation_density_profiles import full_excitation_balance

TEMPERATURE = np.geomspace(1000., 100000., 801)
DENSITY = [1e-6, 100., 1e6, 1e10]
SOURCE = Path.home()/'RT_Codes/CHIANTI/dbase/mg/mg_3/mg_3'
FULL_ATOM = HERE/'nlevel_mg_3_line_strength_candidate.txt'
SMALL_ATOM = HERE/'first27_original_A/nlevel_mg_3.txt'
DRIVER = HERE/'check_source_excitation_atom.x'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, data):
    if path.exists():
        raise FileExistsError(f'Preserve previous result: {path}')
    path.write_text(json.dumps(data, indent=2)+'\n')


def collision_pair(levels, transition):
    return tuple(sorted((transition['ll'], transition['ul']),
                        key=lambda index: (levels[index][1], index)))


def source_atom():
    original = C.read_elvlc(str(SOURCE)+'.elvlc')
    source_radiation = C.read_wgfa(str(SOURCE)+'.wgfa')
    levels, _, collisions = C.bound_excitation_atom('mg', 3)
    radiation = []
    classes = {'first27_original': 0, 'higher_E1_dipole': 0,
               'higher_unidentified_multipole_original': 0}
    for first, second, wavelength, gf, source_rate in source_radiation:
        if first not in levels or second not in levels or source_rate <= 0:
            continue
        lower, upper = sorted((first, second),
                              key=lambda index: (levels[index][1], index))
        physical_gap = levels[upper][1]-levels[lower][1]
        if physical_gap <= 0:
            raise ValueError('Nonpositive radiative photon energy')
        if max(first, second) <= 27:
            rate = source_rate
            classes['first27_original'] += 1
        elif gf > 0:
            source_gap = original[second][2]-original[first][2]
            if source_gap <= 0:
                raise ValueError('Nonpositive theoretical E1 photon energy')
            rate = (source_rate*(physical_gap/source_gap)**3 *
                    original[second][0]/original[upper][0])
            classes['higher_E1_dipole'] += 1
        else:
            rate = source_rate
            classes['higher_unidentified_multipole_original'] += 1
        radiation.append((lower, upper, wavelength, gf, rate))
    if (len(levels), len(radiation), len(collisions)) != (167, 3364, 13861):
        raise ValueError('Full Mg III bound source counts changed')
    if classes != {'first27_original': 117, 'higher_E1_dipole': 3135,
                   'higher_unidentified_multipole_original': 112}:
        raise ValueError('Radiative source classification changed')
    return levels, radiation, collisions


def grid_file(path, mode, temperatures):
    path.write_text(f'{mode} {len(temperatures)}\n' +
                    ' '.join(f'{t:.17g}' for t in temperatures)+'\n')


def pilot():
    levels, radiation, collisions = source_atom()
    temperatures = np.array([1000., 1800., 10000., 100000.])
    source = full_excitation_balance(levels, radiation, collisions,
                                     temperatures, DENSITY)
    minima = []
    for density in DENSITY:
        cooling, populations = source[density]
        if (not np.all(np.isfinite(cooling)) or np.any(cooling < 0)
            or not np.all(np.isfinite(populations))
            or np.min(populations) < -1e-13
            or np.max(abs(populations.sum(axis=1)-1)) > 1e-12):
            raise ValueError('Pilot statistical equilibrium is invalid')
        minima.append(float(np.min(cooling)))
    save_json(HERE/'pilot_step.json', {
        'passed': True, 'temperature_K': temperatures.tolist(),
        'density_cm3': DENSITY, 'full_levels': len(levels),
        'minimum_cooling_erg_cm3_s': min(minima),
    })


def native():
    started = time.monotonic()
    levels, _, collisions = source_atom()
    temperatures = np.array([1000., 1799.9, 1800., 1800.1, 100000.])
    strength_grid = HERE/'profile_strength_grid.txt'
    profile_grid = HERE/'profile_temperature_grid.txt'
    grid_file(strength_grid, 1, temperatures)
    grid_file(profile_grid, 2, TEMPERATURE)
    strengths = HERE/'profile_native_strengths.txt'
    full_output = HERE/'profile_native_full.txt'
    small_output = HERE/'profile_native_first27.txt'
    if any(path.exists() for path in [strengths, full_output, small_output]):
        raise FileExistsError('Preserve existing native profile product')
    subprocess.run([str(DRIVER), str(FULL_ATOM), str(strength_grid), str(strengths)],
                   check=True, timeout=300)
    measured = np.loadtxt(strengths)
    collisions.sort(key=lambda row: collision_pair(levels, row))
    expected = np.column_stack([C.upsilon(row, temperatures) for row in collisions])
    if measured.shape != (len(temperatures), len(collisions)+1):
        raise ValueError('Incomplete native collision strengths')
    if not np.allclose(measured[:, 0], temperatures, rtol=1e-14, atol=0):
        raise ValueError('Native collision-strength temperatures differ')
    positive = expected > 0
    error = float(np.max(np.abs(measured[:, 1:][positive]/expected[positive]-1)))
    if error > 1e-9 or np.any(measured[:, 1:][~positive] != 0):
        raise ValueError(f'Native BT strengths differ: {error}')
    subprocess.run([str(DRIVER), str(FULL_ATOM), str(profile_grid), str(full_output)],
                   check=True, timeout=900)
    subprocess.run([str(DRIVER), str(SMALL_ATOM), str(profile_grid), str(small_output)],
                   check=True, timeout=300)
    for path, number in [(full_output, 167), (small_output, 27)]:
        rows = np.loadtxt(path)
        if rows.shape != (4*len(TEMPERATURE), number+3):
            raise ValueError(f'Incomplete native profile: {path}')
        if not np.allclose(rows[:, 0], np.repeat(TEMPERATURE, 4), rtol=1e-13, atol=0):
            raise ValueError(f'Native temperature grid differs: {path}')
        if not np.array_equal(rows[:, 1], np.tile(DENSITY, len(TEMPERATURE))):
            raise ValueError(f'Native density grid differs: {path}')
        if (not np.all(np.isfinite(rows)) or np.any(rows[:, 2] < 0)
            or np.min(rows[:, 3:]) < -1e-13
            or np.max(np.abs(np.sum(rows[:, 3:], axis=1)-1)) > 1e-12):
            raise ValueError(f'Invalid native cooling or populations: {path}')
    save_json(HERE/'native_step.json', {
        'passed': True, 'scope': 'Compiled source atom and population profile',
        'temperature_points': len(TEMPERATURE), 'density_cm3': DENSITY,
        'full_levels': 167, 'first27_levels': 27,
        'source_strength_samples': int(expected.size),
        'maximum_relative_source_strength_error': error,
        'exact_zero_source_strength_samples': int(np.sum(~positive)),
        'native_seconds': time.monotonic()-started,
        'full_atom_hash': sha256(FULL_ATOM), 'first27_atom_hash': sha256(SMALL_ATOM),
        'driver_hash': sha256(DRIVER), 'native_object_hash': sha256(ROOT/'src/nlevel_mod.o'),
        'full_profile_hash': sha256(full_output), 'first27_profile_hash': sha256(small_output),
    })


def direct():
    started = time.monotonic()
    levels, radiation, collisions = source_atom()
    source = full_excitation_balance(levels, radiation, collisions, TEMPERATURE, DENSITY)
    small_levels = {i: row for i, row in levels.items() if i <= 27}
    small_radiation = [row for row in radiation if row[0] <= 27 and row[1] <= 27]
    small_collisions = [row for row in collisions if row['ll'] <= 27 and row['ul'] <= 27]
    if (len(small_levels), len(small_radiation), len(small_collisions)) != (27, 117, 351):
        raise ValueError('First 27 source subset changed')
    selected = full_excitation_balance(small_levels, small_radiation, small_collisions,
                                       TEMPERATURE, DENSITY)
    payload = {'temperature_K': TEMPERATURE,
               'density_cm3': np.array(DENSITY)}
    minimum = []
    for index, density in enumerate(DENSITY):
        for name, results, n in [('full', source, 167), ('first27', selected, 27)]:
            cooling, populations = results[density]
            if (cooling.shape != (len(TEMPERATURE),)
                or populations.shape != (len(TEMPERATURE), n)
                or not np.all(np.isfinite(cooling)) or np.any(cooling < 0)
                or not np.all(np.isfinite(populations))
                or np.min(populations) < -1e-13
                or np.max(abs(populations.sum(axis=1)-1)) > 1e-12):
                raise ValueError(f'Invalid direct {name} population profile')
            payload[f'cooling_{name}_{index}'] = cooling
            payload[f'populations_{name}_{index}'] = populations
            minimum.append(float(cooling.min()))
    output = HERE/'direct_profiles.npz'
    if output.exists():
        raise FileExistsError('Preserve previous direct population profile')
    np.savez_compressed(output, **payload)
    save_json(HERE/'direct_step.json', {
        'passed': True, 'scope': 'Independent Python statistical equilibrium',
        'temperature_points': len(TEMPERATURE), 'density_cm3': DENSITY,
        'full_levels': len(levels), 'first27_levels': len(small_levels),
        'full_radiative_rows': len(radiation),
        'full_collision_rows': len(collisions),
        'first27_radiative_rows': len(small_radiation),
        'first27_collision_rows': len(small_collisions),
        'minimum_cooling_erg_cm3_s': min(minimum),
        'direct_seconds': time.monotonic()-started,
        'profile_hash': sha256(output),
    })


def compare():
    profile = np.load(HERE/'direct_profiles.npz')
    if not np.array_equal(profile['temperature_K'], TEMPERATURE):
        raise ValueError('Direct temperature grid changed')
    full_native = np.loadtxt(HERE/'profile_native_full.txt')
    small_native = np.loadtxt(HERE/'profile_native_first27.txt')
    native_errors = {'full': [], 'first27': []}
    population_errors = {'full': [], 'first27': []}
    cutoff = []
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, constrained_layout=True)
    difference_fig, difference_axes = plt.subplots(2, 2, figsize=(9, 7),
                                                   sharex=True, constrained_layout=True)
    for index, density in enumerate(DENSITY):
        full = profile[f'cooling_full_{index}']
        selected = profile[f'cooling_first27_{index}']
        for name, rows in [('full', full_native), ('first27', small_native)]:
            measured = rows[index::4, 2]
            expected = profile[f'cooling_{name}_{index}']
            positive = expected > 64*np.finfo(float).tiny
            if np.any(measured[expected == 0] != 0):
                raise ValueError(f'Native {name} zero cooling differs')
            relative = float(np.max(abs(measured[positive]/expected[positive]-1)))
            native_errors[name].append(relative)
            population_errors[name].append(float(np.max(abs(
                rows[index::4, 3:]-profile[f'populations_{name}_{index}']))))
        positive = full > 64*np.finfo(float).tiny
        cutoff.append(float(np.max(abs(selected[positive]/full[positive]-1))))
        ax = axes.flat[index]
        line = ax.loglog(TEMPERATURE, full, label='Corrected full source')[0]
        ax.loglog(TEMPERATURE, full_native[index::4, 2], '--',
                  color=line.get_color(), label='Compiled full source')
        ax.loglog(TEMPERATURE, selected, ':', label='First 27 levels')
        ax.set_title(r'$n_e='+f'{density:g}'+r'$ cm$^{-3}$')
        ax.set_xlabel(r'$T$ [K]')
        ax.set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
        ax.legend(fontsize=8)
        ax = difference_axes.flat[index]
        ax.semilogx(TEMPERATURE[positive],
                    100*(selected[positive]/full[positive]-1))
        ax.axhline(.3, color='black', linewidth=.5)
        ax.axhline(-.3, color='black', linewidth=.5)
        ax.set_title(r'$n_e='+f'{density:g}'+r'$ cm$^{-3}$')
        ax.set_xlabel(r'$T$ [K]')
        ax.set_ylabel(r'First 27 / full - 1 [\%]')
    if max(native_errors['full']+native_errors['first27']) > 1e-8:
        raise ValueError(f'Native cooling differs from direct source: {native_errors}')
    if max(population_errors['full']+population_errors['first27']) > 1e-10:
        raise ValueError(f'Native populations differ from direct source: {population_errors}')
    fig.savefig(HERE/'full_and_first27_cooling_profiles.pdf')
    difference_fig.savefig(HERE/'first27_truncation_differences.pdf')
    plt.close(fig); plt.close(difference_fig)
    save_json(HERE/'comparison_step.json', {
        'passed': True,
        'scope': 'Numerical equivalence and 27-level cutoff within corrected radiative model',
        'temperature_points': len(TEMPERATURE), 'density_cm3': DENSITY,
        'maximum_relative_native_cooling_errors': native_errors,
        'maximum_absolute_native_population_errors': population_errors,
        'first27_full_relative_cooling_errors': cutoff,
        'first27_within_0p3_percent': max(cutoff) < .003,
        'physical_source_accuracy_verified': False,
        'unidentified_multipole_approximation': 'Original higher-level A retained; multipole identity unresolved',
        'full_atom_hash': sha256(FULL_ATOM),
        'first27_atom_hash': sha256(SMALL_ATOM),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', choices=['pilot', 'native', 'direct', 'compare'],
                        required=True)
    args = parser.parse_args()
    {'pilot': pilot, 'native': native, 'direct': direct,
     'compare': compare}[args.step]()


if __name__ == '__main__':
    main()
