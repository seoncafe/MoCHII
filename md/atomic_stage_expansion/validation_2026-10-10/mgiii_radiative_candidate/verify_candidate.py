"""Independently check Mg III dipole-strength physics and native source storage."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'tools/fitting'))
import chianti_cooling as C


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sections(path):
    lines = [line.strip() for line in path.read_text().splitlines()
             if line.strip() and not line.lstrip().startswith('#')]
    if lines.pop(0) != 'NLEVEL_VERSION 2':
        raise ValueError('Missing source version')
    count = int(lines.pop(0).split()[1])
    levels = lines[:count]
    del lines[:count]
    count = int(lines.pop(0).split()[1])
    radiative = lines[:count]
    del lines[:count]
    return levels, radiative, lines


def collision_rows(tail):
    lines = tail[:]
    if lines.pop(0).split()[0] != 'NODES':
        raise ValueError('Missing collision node count')
    count = int(lines.pop(0).split()[1])
    rows = []
    for _ in range(count):
        head = lines.pop(0)
        nodes = int(head.split()[-1])
        rows.append((head, lines[:nodes]))
        del lines[:nodes]
    return rows


def main():
    source = json.loads((HERE.parent/'mgiii_radiative_consistency/result.json').read_text())
    baseline = HERE/'original_A/nlevel_mg_3.txt'
    candidate = HERE/'nlevel_mg_3_line_strength_candidate.txt'
    small = HERE/'first27_original_A/nlevel_mg_3.txt'
    archived = HERE/'rejected_fixed_oscillator_strength/nlevel_mg_3_line_strength_candidate.txt'
    base_levels, base_radiative, base_tail = sections(baseline)
    new_levels, new_radiative, new_tail = sections(candidate)
    small_levels, small_radiative, small_tail = sections(small)
    if new_levels != base_levels or new_tail != base_tail:
        raise ValueError('Energy or collision source changed during radiative correction')
    if new_levels[:27] != small_levels:
        raise ValueError('First 27 level records changed')
    full_collisions = collision_rows(new_tail)
    small_collisions = collision_rows(small_tail)
    first27_collisions = [row for row in full_collisions
                          if max(map(int, row[0].split()[:2])) <= 27]
    if first27_collisions != small_collisions or len(small_collisions) != 351:
        raise ValueError('First 27 collision records changed')
    def rates(rows):
        return {tuple(map(int, line.split()[:2])): float(line.split()[2])
                for line in rows}
    base_rates, new_rates, small_rates = map(rates,
                                             (base_radiative, new_radiative,
                                              small_radiative))
    if set(base_rates) != set(new_rates) or len(new_rates) != 3364:
        raise ValueError('Radiative pair identity changed')
    if {pair: rate for pair, rate in new_rates.items() if max(pair) <= 27} != small_rates:
        raise ValueError('First 27 radiative source records changed')
    original, _, _ = C.read_excitation_atom('mg', 3)
    selected, radiation, collisions = C.bound_excitation_atom('mg', 3)
    source_ids = sorted(selected)
    index = {source_id: i for i, source_id in enumerate(source_ids, 1)}
    identities = {(r['source_lower'], r['source_upper']): r for r in source['rows']}
    dipole_error = []
    unchanged = 0
    reversed_count = 0
    for source_lower, source_upper, wavelength, gf, source_rate in radiation:
        if source_rate <= 0:
            continue
        physical_lower, physical_upper = sorted(
            (source_lower, source_upper), key=lambda i: selected[i][1])
        pair = index[physical_lower], index[physical_upper]
        actual = new_rates[pair]
        if max(source_lower, source_upper) <= 27 or gf == 0:
            if actual != source_rate:
                raise ValueError('Original Jönsson or unidentified-multipole A changed')
            unchanged += 1
            continue
        details = identities[source_lower, source_upper]
        if not details['dipole_selection']:
            raise ValueError('Positive gf without E1 selection')
        source_gap = original[source_upper][2]-original[source_lower][2]
        physical_gap = selected[physical_upper][1]-selected[physical_lower][1]
        if source_gap <= 0 or physical_gap <= 0:
            raise ValueError('Nonpositive E1 photon frequency')
        # Compare the source and candidate reduced electric-dipole strengths,
        # A * g_emitter / photon_frequency**3, independently of the builder.
        source_strength = source_rate*original[source_upper][0]/source_gap**3
        physical_strength = actual*original[physical_upper][0]/physical_gap**3
        dipole_error.append(abs(physical_strength/source_strength-1))
        reversed_count += physical_upper != source_upper
    if len(dipole_error) != 3135 or unchanged != 229 or reversed_count != 16:
        raise ValueError('Unexpected radiative class counts')
    maximum = max(dipole_error)
    if maximum > 2e-14:
        raise ValueError('Electric-dipole line strength was not preserved')
    if sha256(candidate) == sha256(archived):
        raise ValueError('Corrected atom is identical to rejected oscillator model')
    for path, expected in source['source_sha256'].items():
        if sha256(Path(path)) != expected:
            raise ValueError('Original CHIANTI source hash changed')
    # One temperature checks the linked native reader and all unchanged BT
    # source pairs without launching a full finite-density calculation.
    grid = HERE/'one_temperature_grid.txt'
    grid.write_text('1 1\n10000\n')
    output = HERE/'native_final_hook_one_temperature_strengths.txt'
    driver = HERE/'check_source_excitation_atom.x'
    if not output.exists():
        subprocess.run([str(driver), str(candidate), str(grid), str(output)],
                       check=True, timeout=90)
    measured = np.loadtxt(output)
    def pair(transition):
        return tuple(sorted((transition['ll'], transition['ul']),
                            key=lambda i: selected[i][1]))
    collisions.sort(key=pair)
    expected = np.array([C.upsilon(tr, np.array([10000.]))[0]
                         for tr in collisions])
    if measured.shape != (len(collisions)+1,) or measured[0] != 10000.:
        raise ValueError('Native source row count or temperature differs')
    positive = expected > 0
    source_error = float(np.max(np.abs(measured[1:][positive]/expected[positive]-1)))
    if source_error > 1e-9 or np.any(measured[1:][~positive] != 0):
        raise ValueError('Native source strengths differ')
    result = {
        'status': 'corrected preflight only; no finite-density cooling tested',
        'candidate_hash': sha256(candidate),
        'rejected_candidate_hash': sha256(archived),
        'baseline_hash': sha256(baseline),
        'first27_hash': sha256(small),
        'raw_source_hashes_preserved': True,
        'level_records_unchanged': len(new_levels),
        'collision_records_and_nodes_unchanged': True,
        'radiative_rows': len(new_rates),
        'first27_unchanged_rows': len(small_rates),
        'first27_collision_rows_unchanged': len(small_collisions),
        'higher_E1_dipole_strength_rows': len(dipole_error),
        'higher_unidentified_multipole_unchanged_rows': unchanged-len(small_rates),
        'reversed_E1_rows': reversed_count,
        'maximum_relative_dipole_strength_error': maximum,
        'native_source_rows_at_10000K': len(collisions),
        'maximum_relative_native_source_error': source_error,
        'native_object_hash': sha256(ROOT/'src/nlevel_mod.o'),
        'native_driver_hash': sha256(driver),
    }
    (HERE/'corrected_preflight.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
