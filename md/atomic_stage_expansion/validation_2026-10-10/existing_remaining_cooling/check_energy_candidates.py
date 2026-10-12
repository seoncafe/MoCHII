"""Read-only tests of identified CHIANTI energy-vector candidates."""

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'tools/fitting'))
import chianti_cooling as C
from atomic_excitation_energies import selected_energies
from write_excitation_atom import _source_files


def gaps(levels, radiation, collisions):
    bad_a = [(l, u) for l, u, _, _, a in radiation
             if a > 0 and levels[u][1] <= levels[l][1]]
    bad_c = [(row['ll'], row['ul']) for row in collisions
             if levels[row['ul']][1] <= levels[row['ll']][1]]
    return bad_a, bad_c


def candidate(element, stage, observed_through):
    raw, radiation_all, collisions_all = C.read_excitation_atom(element, stage)
    current_bound, radiation, collisions = C.bound_excitation_atom(element, stage)
    current = selected_energies(element, stage, raw, radiation_all,
                                collisions_all)
    hybrid = dict(current)
    for index in current_bound:
        weight, observed, theoretical = raw[index]
        if index <= observed_through and observed >= 0:
            hybrid[index] = (weight, observed, theoretical)
        else:
            hybrid[index] = (weight, theoretical, theoretical)
    bound = {i: hybrid[i] for i in current_bound}
    bad_a, bad_c = gaps(bound, radiation, collisions)
    paths = _source_files(element, stage)
    return {
        'element': element, 'stage': stage,
        'observed_through': observed_through,
        'bound_level_count': len(bound),
        'radiation_rows': len(radiation),
        'electron_rows': len(collisions),
        'nonpositive_radiative_gaps': bad_a,
        'nonpositive_electron_gaps': bad_c,
        'source_sha256': {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in paths},
        'status': 'private energy-identity diagnostic; no source or production change',
    }


def main():
    result = {
        'SiV_observed_1_to_23': candidate('si', 5, 23),
        'MgIII_observed_1_to_49': candidate('mg', 3, 49),
    }
    result['SiV_observed_1_to_23']['gap_screen_passed'] = (
        not result['SiV_observed_1_to_23']['nonpositive_radiative_gaps']
        and not result['SiV_observed_1_to_23']['nonpositive_electron_gaps'])
    result['MgIII_observed_1_to_49']['gap_screen_passed'] = (
        not result['MgIII_observed_1_to_49']['nonpositive_radiative_gaps']
        and not result['MgIII_observed_1_to_49']['nonpositive_electron_gaps'])
    destination = Path(sys.argv[1])
    if destination.exists():
        raise RuntimeError('Preserve the existing energy-candidate diagnostic')
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(destination)


if __name__ == '__main__':
    main()
