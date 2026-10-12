"""Compare published O III optical A values in an isolated six-level atom."""

import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'tools/fitting'))
from verify_nlevel_pyneb import read_nlevel, solve_pops, emissivity, HC_ERG_CM

HERE = Path(__file__).resolve().parent
ORIGINAL = HERE / 'nlevel_o_3_before_correction.txt'
CANDIDATE = HERE / 'nlevel_o_3_storey_zeippen_2000.txt'
# Final journal Table 3, MNRAS 312, 813, in s^-1.
PUBLISHED = {(2, 4): 6.791e-3, (3, 4): 2.046e-2}


def write_candidate():
    lines = ORIGINAL.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('NRAD '))
    count = int(lines[start].split()[1])
    found = set()
    for row in range(start + 1, start + count + 1):
        fields = lines[row].split()
        pair = int(fields[0]), int(fields[1])
        if pair in PUBLISHED:
            lines[row] = f'{pair[0]:5d} {pair[1]:3d}  {PUBLISHED[pair]:.6e}'
            found.add(pair)
    if found != set(PUBLISHED):
        raise ValueError('Source optical pairs are missing')
    lines.insert(1, '# Isolated candidate: O III optical A values from Storey and Zeippen (2000), MNRAS 312, 813, Table 3.')
    lines.insert(2, '# Table 3 totals include relativistically corrected M1 and E2, in s^-1; all other rows match installed source.')
    CANDIDATE.write_text('\n'.join(lines) + '\n')


def cooling(levels, rates, pops, density):
    return sum(pops[u] * a * (levels[u][1] - levels[l][1]) * HC_ERG_CM
               for (l, u), a in rates.items()) / density


def main():
    write_candidate()
    levels, old, collisions = read_nlevel(ORIGINAL)
    new_levels, new, new_collisions = read_nlevel(CANDIDATE)
    assert levels == new_levels and len(collisions) == len(new_collisions)
    assert all(old[pair] == new[pair] for pair in old if pair not in PUBLISHED)
    assert old[(2, 4)] == 6.951e-3 and old[(3, 4)] == 2.029e-2
    temperatures = np.geomspace(1e3, 1e5, 801)
    densities = [1e-6, 1e2, 1e4, 1e6, 1e10]
    profiles = {'temperature_K': temperatures}
    samples = []
    figure, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    for density in densities:
        old_cooling = []
        new_cooling = []
        old_ratio = []
        new_ratio = []
        for temperature in temperatures:
            before = solve_pops(levels, old, collisions, temperature, density)
            after = solve_pops(levels, new, collisions, temperature, density)
            old_cooling.append(cooling(levels, old, before, density))
            new_cooling.append(cooling(levels, new, after, density))
            old_ratio.append(emissivity(levels, old, before, 4, 3) /
                             emissivity(levels, old, before, 4, 2))
            new_ratio.append(emissivity(levels, new, after, 4, 3) /
                             emissivity(levels, new, after, 4, 2))
        old_cooling, new_cooling = np.array(old_cooling), np.array(new_cooling)
        old_ratio, new_ratio = np.array(old_ratio), np.array(new_ratio)
        label = r'$n_e = ' + f'{density:.0e}' + r'\,\mathrm{cm}^{-3}$'
        axes[0].semilogx(temperatures, 100 * (new_cooling / old_cooling - 1), label=label)
        if density == densities[0]:
            axes[1].semilogx(temperatures, old_ratio, '--', label='Installed')
            axes[1].semilogx(temperatures, new_ratio, label='Published Table 3')
        prefix = f'ne_{density:.0e}_'
        profiles[prefix + 'old_cooling_erg_cm3_s'] = old_cooling
        profiles[prefix + 'new_cooling_erg_cm3_s'] = new_cooling
        profiles[prefix + 'old_5007_over_4959'] = old_ratio
        profiles[prefix + 'new_5007_over_4959'] = new_ratio
        samples.append({'electron_density_cm3': density,
                        'max_absolute_cooling_change_percent': float(np.max(np.abs(100 * (new_cooling / old_cooling - 1)))),
                        'T_10000_old_cooling': float(old_cooling[400]),
                        'T_10000_new_cooling': float(new_cooling[400]),
                        'T_10000_old_ratio': float(old_ratio[400]),
                        'T_10000_new_ratio': float(new_ratio[400])})
    axes[0].set_ylabel('Cooling change (percent)')
    axes[0].legend(fontsize=8)
    axes[1].set_ylabel(r'5007 / 4959 energy-flux ratio')
    axes[1].set_xlabel('Temperature (K)')
    axes[1].legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(HERE / 'o3_ground_radiation_profiles.pdf')
    plt.close(figure)
    np.savez_compressed(HERE / 'o3_ground_radiation_profiles.npz', **profiles)
    (HERE / 'o3_ground_radiation_result.json').write_text(json.dumps({
        'source': 'Storey and Zeippen 2000 MNRAS 312 813 Table 3',
        'table_units': 's^-1',
        'published_A_4959_s1': PUBLISHED[(2, 4)],
        'published_A_5007_s1': PUBLISHED[(3, 4)],
        'conditions': 'six-level installed atom, same Upsilon and energies, optically thin, electron collisions only',
        'temperature_grid_K': [1000, 100000, 801],
        'samples': samples}, indent=2) + '\n')
    for row in samples:
        print(row)


if __name__ == '__main__':
    main()
