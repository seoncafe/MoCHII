"""Si I spontaneous probabilities from original NIST 2002 MCHF strengths.

The numerical application interval is 1,000--100,000 K. This interval is an
adopted test domain, not a published accuracy range for the atomic sources.
Original E1 strengths use observed photon frequencies. Existing CHIANTI
non-E1 rates remain approximate effective rates where a source multipole is
not identified. The numerical application interval is not an atomic accuracy
claim; the 2005 journal article's individual transition methods remain unverified.
"""

from collections import Counter, defaultdict
from hashlib import sha256
import math
from pathlib import Path

from silicon_source_records import (CHI_ELVLC_SHA256, CHI_SCUPS_SHA256,
    CHI_WGFA_SHA256, NIST_SHA256, NIST_SOURCE, chianti_levels,
    nist_records, parity)


RADIATIVE_SOURCE_COUNTS = {'E1': 196, 'E2': 13, 'M1': 10, 'M2': 16}
MULTIPOLE_CONSTANT = {'E1': (2.02613e18, 3), 'E2': (1.11995e18, 5),
                      'M1': (2.69735e13, 3), 'M2': (1.4910e13, 5)}
# The first positive BT abscissae imply 997.095--1003.020 K after source
# rounding. This adopted numerical floor is not an atomic accuracy limit.
INFERRED_BT_FLOOR_K = 1003.02


def _sha(path):
    return sha256(path.read_bytes()).hexdigest()


def _radiative_rate(record, emitting_weight, gap_cm_inverse):
    constant, exponent = MULTIPOLE_CONSTANT[record['multipole']]
    if gap_cm_inverse <= 0 or emitting_weight <= 0:
        raise ValueError('Nonpositive photon energy or emitting weight')
    wavelength = 1e8 / gap_cm_inverse
    rate = constant * record['length_S_atomic_units'] / (
        emitting_weight * wavelength**exponent)
    if not math.isfinite(rate) or rate <= 0:
        raise ValueError('Invalid transformed spontaneous probability')
    return rate


def _first_positive_node_temperature(transition):
    from chianti_cooling import KB_OVER_RY
    scaled = next(float(x) for x in transition['xs'] if 0 < x < 1)
    scale = float(transition['cups'])
    if transition['ttype'] in (1, 4):
        energy_ratio = math.exp(math.log(scale) / (1 - scaled)) - scale
    else:
        energy_ratio = scaled * scale / (1 - scaled)
    return energy_ratio * float(transition['de']) / KB_OVER_RY


def bound_sii_spontaneous_probabilities(raw, selected, original_radiation,
                                        collisions, stem):
    """Correct selected bound Si I radiation; retain all collision nodes.

    This function is called after the ordinary bound-state selection. It does
    not call that selector and therefore cannot recurse through its dispatch.
    """
    source_stem = Path(stem)
    source_dir = source_stem.parent
    inputs = {source_dir / 'si_1.elvlc': CHI_ELVLC_SHA256,
              source_dir / 'si_1.wgfa': CHI_WGFA_SHA256,
              source_dir / 'si_1.scups': CHI_SCUPS_SHA256,
              NIST_SOURCE: NIST_SHA256}
    if any(_sha(path) != expected for path, expected in inputs.items()):
        raise ValueError('Si I original source bytes changed')
    levels, identities = chianti_levels(source_dir / 'si_1.elvlc')
    if set(levels) != set(range(1, 37)) or set(raw) != set(levels):
        raise ValueError('Si I state identities changed')
    records, unmapped = nist_records(NIST_SOURCE, identities)
    if unmapped or Counter(row['multipole'] for row in records) != RADIATIVE_SOURCE_COUNTS:
        raise ValueError('Si I original transition identities changed')
    if (set(selected) != set(range(1, 36)) or len(original_radiation) != 229
            or len(collisions) != 595 or sum(len(row['xs']) for row in collisions) != 8925):
        raise ValueError('Si I selected source size changed')
    first_temperatures = [_first_positive_node_temperature(row) for row in collisions]
    if (not all(math.isfinite(value) and value > 0 for value in first_temperatures)
            or not (997.0 < min(first_temperatures) < 998.0)
            or not (1003.0 < max(first_temperatures) < INFERRED_BT_FLOOR_K)):
        raise ValueError('Si I rounded BT finite-node temperatures changed')
    for row in collisions:
        if row.get('bt_tmin_K') != INFERRED_BT_FLOOR_K:
            raise ValueError('Si I adopted BT floor changed')
    if any((selected[index][0] != levels[index]['identity'][2] + 1 or
            selected[index][1] != levels[index]['observed_cm_inverse'])
           for index in selected):
        raise ValueError('Si I selected state identity or energy changed')
    if (selected[1][1] != 0 or any(not math.isfinite(weight) or weight <= 0
                                  or not math.isfinite(energy) or energy <= 0
                                  for index, (weight, energy) in selected.items()
                                  if index != 1)):
        raise ValueError('Si I selected physical state values changed')
    source_pairs = defaultdict(list)
    for record in records:
        first, second = record['first_chianti_id'], record['second_chianti_id']
        if first not in selected or second not in selected:
            raise ValueError('Si I source transition outside bound states')
        first_state, second_state = levels[first]['identity'], levels[second]['identity']
        if (record['first_g'] != selected[first][0]
                or record['second_g'] != selected[second][0]
                or not math.isfinite(record['source_gap_cm_inverse'])
                or record['source_gap_cm_inverse'] <= 0):
            raise ValueError('Si I NIST transition energy or weight changed')
        odd = parity(first_state) != parity(second_state)
        rank = 1 if record['multipole'] in ('E1', 'M1') else 2
        if (odd != (record['multipole'] in ('E1', 'M2'))
                or abs(first_state[2] - second_state[2]) > 2 * rank
                or (first_state[2] == 0 and second_state[2] == 0)):
            raise ValueError('Si I original multipole state selection changed')
        source_pairs[tuple(sorted((first, second)))].append(record)
        source_rate = _radiative_rate(record, record['second_g'],
                                      1e8 / record['source_wavelength_angstrom'])
        if abs(source_rate / record['length_A_per_second'] - 1) > 4e-5:
            raise ValueError('Si I NIST multipole strength convention changed')
    if len(source_pairs) != 216:
        raise ValueError('Si I original source pair count changed')

    emitted = []
    seen = set()
    matched_e1 = 0
    retained_non_e1 = 0
    for first, second, _, _, original_a in original_radiation:
        pair = tuple(sorted((first, second)))
        if pair in seen:
            raise ValueError('Duplicate CHIANTI radiation pair')
        seen.add(pair)
        matching = source_pairs.get(pair, [])
        e1 = [row for row in matching if row['multipole'] == 'E1']
        if e1:
            if len(e1) != 1:
                raise ValueError('Ambiguous NIST E1 pair')
            matched_e1 += 1
            lower, upper = sorted(pair, key=lambda index: selected[index][1])
            if (parity(levels[first]['identity']) == parity(levels[second]['identity'])
                    or abs(levels[first]['identity'][2] - levels[second]['identity'][2]) > 2):
                raise ValueError('Unphysical E1 state identity')
            gap = selected[upper][1] - selected[lower][1]
            rate = sum(_radiative_rate(row, selected[upper][0], gap)
                       for row in matching)
            emitted.append((lower, upper, 1e8 / gap, 0., rate))
        else:
            retained_non_e1 += 1
            if original_a <= 0 or selected[second][1] <= selected[first][1]:
                raise ValueError('Non-E1 source rate has invalid energy direction')
            # Its CHIANTI gf is not used as an E1 absorption strength.
            emitted.append((first, second, 1e8 /
                            (selected[second][1] - selected[first][1]), 0., original_a))
    if matched_e1 != 185 or retained_non_e1 != 44:
        raise ValueError('Si I identified source radiation count changed')

    added = Counter()
    for pair, matching in sorted(source_pairs.items()):
        if pair in seen:
            continue
        lower, upper = sorted(pair, key=lambda index: selected[index][1])
        gap = selected[upper][1] - selected[lower][1]
        if gap <= 0:
            raise ValueError('Nonpositive added photon energy')
        rate = sum(_radiative_rate(row, selected[upper][0], gap)
                   for row in matching)
        emitted.append((lower, upper, 1e8 / gap, 0., rate))
        added.update(row['multipole'] for row in matching)
    if added != {'E1': 11, 'E2': 6, 'M1': 5, 'M2': 4} or len(emitted) != 250:
        raise ValueError('Si I completed radiation count changed')
    if len({tuple(sorted(row[:2])) for row in emitted}) != len(emitted):
        raise ValueError('Duplicate completed radiation pair')
    return emitted


