"""Prepare and verify Mg I with the adopted cold constant-Upsilon limit.

The supplied v2 atom contains one collision temperature, 10,000 K. Below it,
the user-approved constant-strength approximation defines a restricted atom.
No high-temperature continuation or published accuracy claim is supplied.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mocassin_v2_atom_inventory import read_atom, MOCASSIN
from validate_source_excitation_atom import compile_driver, grid_file, DENSITY
from validate_excitation_density_profiles import full_excitation_balance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    if work.exists():
        raise FileExistsError('Preserve earlier calculations')
    work.mkdir(parents=True)
    source = MOCASSIN/'data/mgi.dat'
    atom = read_atom(source, include_tables=True)
    if atom['temperatures_K'] != [10000.] or atom['nlev'] != 5:
        raise ValueError('Unexpected Mg I source interval or state count')
    if atom['nonpositive_photon_energy_pairs']:
        raise ValueError('Nonpositive radiative energy')
    levels = {r['index']: (r['statistical_weight'], r['energy_cm_inverse'])
              for r in atom['levels']}
    radiation = [(r['lower'], r['upper'], 0., 0., r['A_s_inverse'])
                 for r in atom['spontaneous_probabilities'] if r['A_s_inverse'] > 0]
    rows = sorted(atom['collision_tables'], key=lambda r: (r['lower'], r['upper']))
    if any(len(r['upsilon']) != 1 or r['upsilon'][0] <= 0 for r in rows):
        raise ValueError('Require positive single-node strengths')
    text = ['# Restricted Mg I cold constant-Upsilon approximation, not an accuracy certificate.',
            '# Source has exactly one collision temperature: 10000 K.',
            '# The 1000 K representation node is synthetic and holds the supplied strength constant.',
            '# Valid numerical preparation interval: 1000 <= T <= 10000 K.',
            '# Below 1000 K the native cold rule still holds Upsilon constant; not validated here.',
            '# Above 10000 K reject evaluation; no high-temperature model is specified.',
            '# Radiative probabilities are retained from the supplied source, with its accuracy unverified.',
            f'# Source SHA-256: {hashlib.sha256(source.read_bytes()).hexdigest()}',
            'NLEVEL_VERSION 2', 'NLEV 5']
    for i, (g, energy) in levels.items():
        text.append(f'{i} {g:.17g} {energy:.17g}')
    text.append(f'NRAD {len(radiation)}')
    text += [f'{lo} {up} {prob:.17g}' for lo, up, _, _, prob in radiation]
    text += [f'NODES {2*len(rows)}', f'NUPS {len(rows)}']
    collisions = []
    for r in rows:
        lo, up, strength = r['lower'], r['upper'], r['upsilon'][0]
        if levels[up][1] <= levels[lo][1]:
            raise ValueError('Collision endpoints are not energy ordered')
        text += [f'{lo} {up} TEMPERATURE_TABLE 0 0 0 LOG_PCHIP_T 2',
                 f'1000 {strength:.17g}', f'10000 {strength:.17g}']
        collisions.append(dict(ll=lo, ul=up, temperature_grid=np.array([1000., 10000.]),
                               upsilon_grid=np.array([strength, strength])))
    text += ['# Complete only for this supplied five-level approximation.', 'DENSITY_RATIO 1']
    native = work/'nlevel_mg_1.txt'
    native.write_text('\n'.join(text)+'\n')
    exe = compile_driver(work)
    temperature = np.geomspace(1000., 10000., 101)
    grid_file(work/'grid.in', 2, temperature)
    subprocess.run([str(exe), str(native), str(work/'grid.in'), str(work/'profile.txt')],
                   check=True, timeout=120)
    measured = np.loadtxt(work/'profile.txt')
    if measured.shape != (404, 8) or not np.isfinite(measured).all():
        raise ValueError('Incomplete or nonfinite Mg I profiles')
    independent = full_excitation_balance(levels, radiation, collisions, temperature, DENSITY)
    cool_error = population_error = 0.
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for j, ne in enumerate(DENSITY):
        power, pop = independent[ne]
        rows = measured[j::4]
        if np.any(power <= 0):
            raise ValueError('Nonpositive cooling')
        difference = rows[:, 2]/power-1
        cool_error = max(cool_error, float(np.max(np.abs(difference))))
        population_error = max(population_error, float(np.max(np.abs(rows[:, 3:]-pop))))
        label = r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$'
        axes.flat[0].loglog(temperature, rows[:, 2], label=label)
        axes.flat[0].loglog(temperature, power, 'k:', linewidth=.6)
        ax = axes.flat[j+1]
        for k in range(5):
            ax.loglog(temperature, np.maximum(rows[:, 3+k], 1e-30))
            ax.loglog(temperature, np.maximum(pop[:, k], 1e-30), 'k:', linewidth=.5)
        ax.set_title(label)
        ax.set_ylabel('All level populations')
        axes.flat[5].semilogx(temperature, difference, label=label)
    axes.flat[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
    axes.flat[0].legend(fontsize=7)
    axes.flat[5].set_ylabel('Relative cooling difference')
    for ax in axes.flat:
        ax.set_xlabel(r'$T$ [K]')
    fig.suptitle('Mg I cold approximation: native solid / independent dotted')
    fig.savefig(work/'cold_excitation_profiles.pdf')
    plt.close(fig)
    grid_file(work/'above_source.in', 1, np.array([10000.1]))
    rejected = subprocess.run([str(exe), str(native), str(work/'above_source.in'),
                               str(work/'above_source.txt')], capture_output=True, text=True, timeout=30)
    (work/'above_source.log').write_text(rejected.stdout+rejected.stderr)
    if rejected.returncode == 0 or 'temperature exceeds collision table' not in rejected.stdout+rejected.stderr:
        raise ValueError('Native evaluator did not reject unsupported hot temperature')
    if cool_error > 1e-8 or population_error > 1e-8:
        raise ValueError('Native and independent excitation solutions differ')
    result = dict(passed=True, samples=404, nlev=5, nrad=len(radiation), nups=len(collisions),
                  source_temperature_K=10000., numerical_interval_K=[1000., 10000.],
                  cooling_relative_error=cool_error, population_absolute_error=population_error,
                  hot_evaluation_rejected=True, production_installed=False, visual_review_completed=False,
                  atomic_accuracy_verified=False,
                  source_sha256={str(f): hashlib.sha256(f.read_bytes()).hexdigest()
                                 for f in (source, Path(__file__), native, exe)})
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
