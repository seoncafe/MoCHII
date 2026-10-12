"""Conditional Si I E1 length- and velocity-form finite-density comparison.

The original velocity operator omits some relativistic terms when term mixing
is strong. The difference between its rates and the original length-form
rates is a sensitivity diagnostic, not an atomic uncertainty interval.
"""

import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FITTING = ROOT / 'tools' / 'fitting'
sys.path.insert(0, str(FITTING))
SOURCE_DIR = HERE.parent / 'si_i_electric_dipole_source_audit'
PROFILE_DIR = HERE.parent / 'si_i_radiative_completion'
sys.path.insert(0, str(SOURCE_DIR))

from audit_si_i_electric_dipole import (CHI_DIR, DEFAULT_NIST, chianti_levels,
                                        nist_records)
from validate_excitation_density_profiles import full_excitation_balance
from si_i_spontaneous_probabilities import MULTIPOLE_CONSTANT


TEMPERATURE = np.geomspace(1000.0, 100000.0, 401)
DENSITY = np.array([1e-6, 100.0, 1e6, 1e10])
CONSTANT_E1 = 2.02613e18


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_atom():
    summary = json.loads((SOURCE_DIR / 'candidate_v2/candidate_result.json').read_text())
    atom_path = Path(summary['candidate_path'])
    module_path = Path(summary['source_module_path'])
    if (summary.get('passed') is not True or summary['nlev'] != 35
            or summary['nrad'] != 250 or summary['nups'] != 595
            or summary['nodes'] != 8925 or summary['bt_tmin_K'] != 1003.02
            or digest(atom_path) != summary['candidate_sha256']
            or digest(module_path) != summary['source_sha256'][str(module_path)]):
        raise ValueError('Frozen Si I candidate source changed')
    for path, expected in summary['source_sha256'].items():
        if digest(path) != expected:
            raise ValueError(f'Si I candidate input changed: {path}')
    spec = importlib.util.spec_from_file_location('si_i_sensitivity_source', module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return summary, *module.source_atom()


def velocity_radiation(levels, radiation):
    identities = chianti_levels(CHI_DIR / 'si_1.elvlc')[1]
    records, unmapped = nist_records(DEFAULT_NIST, identities)
    if unmapped:
        raise ValueError('Unmapped NIST source state')
    source = {tuple(sorted((lower, upper))): (lower, upper, wl, gf, rate)
              for lower, upper, wl, gf, rate in radiation}
    if len(source) != 250:
        raise ValueError('Si I radiation pairs changed')
    changed = {}
    for record in records:
        if record['multipole'] != 'E1':
            continue
        lower, upper = sorted((record['first_chianti_id'],
                               record['second_chianti_id']),
                              key=lambda index: levels[index][1])
        pair = (lower, upper)
        if pair not in source or pair in changed:
            raise ValueError('Duplicate or missing E1 source pair')
        gap = levels[upper][1] - levels[lower][1]
        if gap <= 0 or record.get('velocity_S_atomic_units', 0) <= 0:
            raise ValueError('Unusable NIST velocity-form strength')
        factor = CONSTANT_E1 / (levels[upper][0] * (1e8 / gap)**3)
        length_rate = factor * record['length_S_atomic_units']
        velocity_rate = factor * record['velocity_S_atomic_units']
        if not np.isfinite(length_rate + velocity_rate) or min(length_rate, velocity_rate) <= 0:
            raise ValueError('Nonfinite E1 gauge rate')
        other_rates = []
        for other in records:
            if (other['multipole'] != 'E1' and
                    tuple(sorted((other['first_chianti_id'], other['second_chianti_id'])))
                    == tuple(sorted(pair))):
                constant, exponent = MULTIPOLE_CONSTANT[other['multipole']]
                other_rates.append(constant * other['length_S_atomic_units'] /
                                   (levels[upper][0] * (1e8 / gap)**exponent))
        other_rate = sum(other_rates)
        if abs(source[tuple(sorted(pair))][4] / (length_rate + other_rate) - 1) > 1e-12:
            raise ValueError('Identified multipoles do not recover the length-form total')
        changed[pair] = (length_rate, velocity_rate, other_rate)
    if len(changed) != 196:
        raise ValueError('Si I E1 gauge coverage changed')
    alternative = []
    for lower, upper, wavelength, gf, rate in radiation:
        pair = (lower, upper)
        if pair in changed:
            _, velocity_rate, other_rate = changed[pair]
            # Sum identified positive contributions without subtracting a
            # dominant E1 probability to recover a weak M2 probability.
            rate = velocity_rate + other_rate
        if not np.isfinite(rate) or rate <= 0:
            raise ValueError('Velocity-form total rate is invalid')
        alternative.append((lower, upper, wavelength, gf, rate))
    if any(levels[upper][1] <= levels[lower][1]
           for lower, upper, _, _, _ in alternative):
        raise ValueError('Nonphysical emission direction')
    return alternative, changed


def baseline_profiles(summary):
    report_path = PROFILE_DIR / 'direct_step.json'
    profile_path = PROFILE_DIR / 'direct_profiles.npz'
    report = json.loads(report_path.read_text())
    if (report.get('passed') is not True
            or report.get('candidate_sha256') != summary['candidate_sha256']
            or report.get('temperature_points') != 401
            or report.get('densities_cm3') != DENSITY.tolist()
            or digest(profile_path) != report['profile_sha256']):
        raise ValueError('Frozen Si I length-form profile changed')
    with np.load(profile_path) as saved:
        if (not np.array_equal(saved['temperature_K'], TEMPERATURE)
                or not np.array_equal(saved['density_cm3'], DENSITY)):
            raise ValueError('Si I reference grid changed')
        return {density: (saved[f'cooling_{index}'].copy(),
                          saved[f'populations_{index}'].copy())
                for index, density in enumerate(DENSITY)}


def check_profiles(profiles):
    for density in DENSITY:
        cooling, populations = profiles[density]
        if (cooling.shape != (401,) or populations.shape != (401, 35)
                or not np.all(np.isfinite(cooling)) or np.any(cooling <= 0)
                or not np.all(np.isfinite(populations))
                or np.min(populations) < -1e-13
                or np.max(np.abs(populations.sum(axis=1) - 1)) > 1e-12):
            raise ValueError('Invalid Si I gauge population or cooling profile')


def make_pdf(path, baseline, velocity):
    with PdfPages(path) as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True,
                                 constrained_layout=True)
        for index, density in enumerate(DENSITY):
            ax = axes.flat[index]
            original = baseline[density][0]
            measured = velocity[density][0]
            ax.loglog(TEMPERATURE, original, label='Length form')
            ax.loglog(TEMPERATURE, measured, '--', label='Velocity form')
            ax.set_title(r'$n_e = ' + f'{density:g}' + r'\;\mathrm{cm}^{-3}$')
            ax.set_xlabel(r'$T\;[\mathrm{K}]$')
            ax.set_ylabel(r'$\Lambda\;[\mathrm{erg\,cm^3\,s^{-1}}]$')
            ax.legend(fontsize=8)
        pdf.savefig(fig)
        plt.close(fig)
        fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True,
                                 constrained_layout=True)
        for index, density in enumerate(DENSITY):
            ax = axes.flat[index]
            relative = 100 * (velocity[density][0] / baseline[density][0] - 1)
            ax.semilogx(TEMPERATURE, relative)
            ax.set_title(r'$n_e = ' + f'{density:g}' + r'\;\mathrm{cm}^{-3}$')
            ax.set_xlabel(r'$T\;[\mathrm{K}]$')
            ax.set_ylabel(r'$100(\Lambda_v/\Lambda_l-1)$')
        pdf.savefig(fig)
        plt.close(fig)
        for density in DENSITY:
            fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                     constrained_layout=True)
            original = baseline[density][1]
            measured = velocity[density][1]
            for state in range(35):
                axes[0].semilogy(TEMPERATURE,
                                 np.maximum(original[:, state], 1e-320),
                                 linewidth=0.6, alpha=0.65)
            image = axes[1].pcolormesh(TEMPERATURE, np.arange(1, 36),
                                        (measured - original).T,
                                        shading='auto', cmap='coolwarm')
            axes[1].set_xscale('log')
            axes[0].set_ylabel('Length-form population')
            axes[1].set_ylabel('State index')
            axes[1].set_xlabel(r'$T\;[\mathrm{K}]$')
            fig.colorbar(image, ax=axes[1], label='Velocity minus length population')
            fig.suptitle(r'$n_e = ' + f'{density:g}' + r'\;\mathrm{cm}^{-3}$')
            pdf.savefig(fig)
            plt.close(fig)


def run():
    started = time.monotonic()
    result_path = HERE / 'result.json'
    profile_path = HERE / 'full_profiles.npz'
    figure_path = HERE / 'full_profiles.pdf'
    if any(path.exists() for path in (result_path, profile_path, figure_path)):
        raise FileExistsError('Preserve previous Si I gauge sensitivity outputs')
    summary, levels, radiation, collisions = source_atom()
    changed_radiation, identified = velocity_radiation(levels, radiation)
    original = baseline_profiles(summary)
    measured = full_excitation_balance(levels, changed_radiation,
                                       collisions, TEMPERATURE, DENSITY)
    check_profiles(original)
    check_profiles(measured)
    payload = {'temperature_K': TEMPERATURE, 'density_cm3': DENSITY}
    differences = []
    for index, density in enumerate(DENSITY):
        baseline_cooling, baseline_pop = original[density]
        measured_cooling, measured_pop = measured[density]
        payload[f'length_cooling_{index}'] = baseline_cooling
        payload[f'velocity_cooling_{index}'] = measured_cooling
        payload[f'length_populations_{index}'] = baseline_pop
        payload[f'velocity_populations_{index}'] = measured_pop
        relative = measured_cooling / baseline_cooling - 1
        maximum = int(np.argmax(np.abs(relative)))
        differences.append({
            'density_cm3': float(density),
            'maximum_absolute_cooling_difference_percent':
                float(100 * np.max(np.abs(relative))),
            'maximum_cooling_difference_temperature_K': float(TEMPERATURE[maximum]),
            'maximum_absolute_population_difference':
                float(np.max(np.abs(measured_pop - baseline_pop))),
        })
    np.savez_compressed(profile_path, **payload)
    make_pdf(figure_path, original, measured)
    with np.load(profile_path) as saved:
        if (set(saved.files) != set(payload)
                or any(not np.array_equal(saved[key], value)
                       for key, value in payload.items())):
            raise ValueError('Incomplete saved Si I gauge profiles')
    pdf = subprocess.run(['pdfinfo', str(figure_path)], check=True,
                         capture_output=True, text=True, timeout=30).stdout
    if re.search(r'^Pages:\s+6\s*$', pdf, re.MULTILINE) is None or figure_path.stat().st_size < 10000:
        raise ValueError('Incomplete Si I gauge profile PDF')
    result_path.write_text(json.dumps({
        'passed': True,
        'interpretation': 'Conditional E1 gauge sensitivity, not an atomic accuracy interval',
        'candidate_sha256': summary['candidate_sha256'],
        'source_module_sha256': summary['source_sha256'][summary['source_module_path']],
        'original_nist_sha256': summary['source_sha256'][str(DEFAULT_NIST.resolve())],
        'temperature_points': 401, 'densities_cm3': DENSITY.tolist(),
        'changed_E1_pairs': len(identified),
        'unchanged_source_collision_rows': len(collisions),
        'numerical_floor_K': 1003.02,
        'differences': differences,
        'profiles_path': str(profile_path.resolve()),
        'profiles_sha256': digest(profile_path),
        'figure_path': str(figure_path.resolve()),
        'figure_sha256': digest(figure_path),
        'seconds': time.monotonic() - started,
    }, indent=2) + '\n')


if __name__ == '__main__':
    run()
