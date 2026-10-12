"""Inspect the local OPEN-ADAS neutral-argon LS ADF04 without copying its data."""

import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path


SOURCE = Path('/tmp/arlike_cpb12_ar0.dat')
OUTPUT = Path(__file__).with_name('result.json')


def adas_number(value):
    return float(re.sub(r'(?<=\d)([+-]\d+)$', r'E\1', value))


def main():
    lines = SOURCE.read_text().splitlines()
    header = lines[0].split()
    assert header[:4] == ['Ar+', '0', '18', '1']
    ionization_boundary = float(header[4].split('(')[0])
    level_lines = lines[1:36]
    levels = []
    for line in level_lines:
        match = re.match(r'^\s*(\d+)\s+(.+?)\s+\((\d)\)(\w)\(\s*([\d.]+)\)\s+([\d.]+)\s*$', line)
        assert match, line
        index, configuration, multiplicity, orbital_term, xj, energy = match.groups()
        weight = 2 * float(xj) + 1
        levels.append({'id': int(index), 'configuration': configuration,
                       'multiplicity': int(multiplicity), 'term': orbital_term,
                       'weight': weight, 'energy_cm1': float(energy)})
    assert [x['id'] for x in levels] == list(range(1, 36))
    assert lines[36].strip() == '-1'
    grid = [adas_number(x) for x in lines[37].split()[2:]]
    assert len(grid) == 14 and all(b > a > 0 for a, b in zip(grid, grid[1:]))
    rows = []
    for line in lines[38:]:
        if line.strip() == '-1':
            break
        fields = line.split()
        assert len(fields) == 18, (line, len(fields))
        upper, lower = map(int, fields[:2])
        rate, *strengths = [adas_number(x) for x in fields[2:]]
        assert 1 <= lower < upper <= 35
        assert len(strengths) == 15
        rows.append({'upper': upper, 'lower': lower, 'A_s1': rate,
                     'upsilon': strengths[:14], 'bethe': strengths[14]})
    pairs = [(x['lower'], x['upper']) for x in rows]
    assert len(rows) == 595 and len(set(pairs)) == 595
    assert set(pairs) == {(l, u) for l in range(1, 36) for u in range(l + 1, 36)}
    assert all(math.isfinite(x['A_s1']) and x['A_s1'] >= 0 for x in rows)
    assert all(math.isfinite(y) and y >= 0 for x in rows for y in x['upsilon'])
    assert all(math.isfinite(x['bethe']) and x['bethe'] >= 0 for x in rows)
    outgoing = Counter(x['upper'] for x in rows if x['A_s1'] > 0)
    radiative_zero = [x['id'] for x in levels if outgoing[x['id']] == 0]
    energy = {x['id']: x['energy_cm1'] for x in levels}
    nonpositive_gap_A = [[x['lower'], x['upper'], x['A_s1']]
                         for x in rows if x['A_s1'] > 0
                         and energy[x['upper']] <= energy[x['lower']]]
    degenerate_collision = [[l, u] for l, u in pairs if energy[l] == energy[u]]
    edges = {x['id']: set() for x in levels}
    for x in rows:
        if x['A_s1'] > 0 and energy[x['upper']] > energy[x['lower']]:
            edges[x['upper']].add(x['lower'])
    can_reach_ground = {1}
    while True:
        enlarged = can_reach_ground | {u for u, ls in edges.items() if ls & can_reach_ground}
        if enlarged == can_reach_ground:
            break
        can_reach_ground = enlarged
    can_reach_either_sink = {1, 2}
    while True:
        enlarged = can_reach_either_sink | {
            u for u, ls in edges.items() if ls & can_reach_either_sink}
        if enlarged == can_reach_either_sink:
            break
        can_reach_either_sink = enlarged
    outlier = next(x for x in rows if (x['lower'], x['upper']) == (8, 15))
    tindex = grid.index(10000.0)
    upsilon_outlier = outlier['upsilon'][tindex]
    delta_cm1 = energy[15] - energy[8]
    rate_outlier = (8.629e-6 * upsilon_outlier /
                    (levels[7]['weight'] * math.sqrt(10000.0)) *
                    math.exp(-delta_cm1 / (0.6950348 * 10000.0)))
    result = {
        'source': str(SOURCE),
        'sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'header': header,
        'levels': len(levels),
        'temperature_grid_K': grid,
        'energy_range_cm1': [min(energy.values()), max(energy.values())],
        'ionization_boundary_cm1_file': ionization_boundary,
        'all_levels_below_file_boundary': all(
            x['energy_cm1'] < ionization_boundary for x in levels),
        'weight_range': [min(x['weight'] for x in levels), max(x['weight'] for x in levels)],
        'ls_weight_matches_multiplicity_term': all(
            x['weight'] == x['multiplicity'] * (2 * int(x['term']) + 1)
            for x in levels),
        'pairs': len(rows),
        'positive_A_rows': sum(x['A_s1'] > 0 for x in rows),
        'zero_A_rows': sum(x['A_s1'] == 0 for x in rows),
        'nonpositive_energy_gap_positive_A': nonpositive_gap_A,
        'degenerate_collision_pairs': degenerate_collision,
        'radiative_zero_outgoing_levels': radiative_zero,
        'levels_without_positive_A_path_to_ground': sorted(set(range(1, 36)) - can_reach_ground),
        'levels_without_positive_A_path_to_either_sink': sorted(
            set(range(1, 36)) - can_reach_either_sink),
        'min_positive_A_s1': min(x['A_s1'] for x in rows if x['A_s1'] > 0),
        'max_A_s1': max(x['A_s1'] for x in rows),
        'zero_upsilon_nodes': sum(y == 0 for x in rows for y in x['upsilon']),
        'pairs_with_zero_upsilon_node': sum(any(y == 0 for y in x['upsilon']) for x in rows),
        'min_positive_upsilon': min(y for x in rows for y in x['upsilon'] if y > 0),
        'max_upsilon': max(y for x in rows for y in x['upsilon']),
        'outlier_pair_8_15_at_10000_K': {
            'lower_weight': levels[7]['weight'],
            'energy_gap_cm1': delta_cm1,
            'upsilon': upsilon_outlier,
            'Maxwellian_excitation_rate_cm3_s': rate_outlier,
        },
        'pairs_upsilon_above_1e5_at_10000_K': sum(
            x['upsilon'][tindex] > 1e5 for x in rows),
        'zero_bethe_values': sum(x['bethe'] == 0 for x in rows),
        'source_file_explicit_ionization_or_recombination_rows': False,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
