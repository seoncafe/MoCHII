"""Restore two identified Ne I spontaneous channels omitted by CHIANTI."""

from hashlib import sha256
import math
from pathlib import Path
import re

CHIANTI_SHA256 = (
    '7e813574e3157b9735a7d89660abcf463019e47f1df228ad41ad8c84b6167788',
    '92ae3ab4b0ce76f358617a0ba8d30c0758f0bcfbd7e0529d85ab4eca7509bd8f',
    'd51fbb34fb71c533a53668e3425a213672ddf66eaf0c7fc51c0b72b6363e4b06')
STOUT_SHA256 = (
    '6f0e609cb82255357cf870504a9a423c04f0da7f90b5dacebf8619a0748048e4',
    '1aa24d2f7125780e258baea02b9c7439f7b6d5251176a7fe2611be05ea9b88eb')
STOUT_STEM = Path.home() / 'CLOUDY/c25.00/data/stout/ne/ne_1/ne_1'
IDENTIFIED_RATES = {(3, 15): 9.00e5, (5, 15): 6.82e7}


def nei_identified_spontaneous_probabilities(levels, radiation, collisions,
                                              source_stem):
    """Return the CHIANTI Ne I rows plus two NIST-dated Stout decay channels.

    The first 31 Stout states match the CHIANTI observed energies, weights,
    orbital configurations, and J. CHIANTI omits both 3p 1S0 (level 15)
    decays to 3s J=1 states. Preserve all 142 original A values and collision
    nodes; Stout supplies only the two missing A values here. The effective
    collision strengths are adopted as an explicit 1,000--100,000 K numerical
    approximation; the source conversion and atomic accuracy are unverified.
    Below the inferred 1000.36 K finite-node floor, Upsilon is held constant
    while actual gas temperature enters Maxwellian rates. Neither floor nor
    numerical use interval is a published physical validity limit.
    """
    chianti_paths = [Path(f'{source_stem}.{ext}')
                     for ext in ('elvlc', 'wgfa', 'scups')]
    stout_paths = [Path(f'{STOUT_STEM}.{ext}') for ext in ('nrg', 'tp')]
    if tuple(sha256(path.read_bytes()).hexdigest() for path in chianti_paths) != CHIANTI_SHA256:
        raise ValueError('Ne I CHIANTI source identity changed')
    if tuple(sha256(path.read_bytes()).hexdigest() for path in stout_paths) != STOUT_SHA256:
        raise ValueError('Ne I Stout source identity changed')
    if (set(levels) != set(range(1, 32)) or len(radiation) != 142
            or len(collisions) != 465
            or sum(len(row['xs']) for row in collisions) != 6975):
        raise ValueError('Ne I CHIANTI source counts changed')

    chianti_lines = {}
    for line in chianti_paths[0].read_text().splitlines():
        if line.strip() == '-1':
            break
        chianti_lines[int(line[:7])] = line
    stout_lines = stout_paths[0].read_text().splitlines()[1:32]
    if len(stout_lines) != 31:
        raise ValueError('Ne I Stout bound identity list changed')
    for index, line in enumerate(stout_lines, 1):
        fields = line.split(maxsplit=3)
        if len(fields) != 4 or int(fields[0]) != index:
            raise ValueError('Ne I Stout level order changed')
        energy, weight = float(fields[1]), int(fields[2])
        source = chianti_lines[index]
        configuration = source[7:37].strip()
        label = fields[3].strip('"')
        orbital = re.search(r'\.(\d+[spd])\s', label)
        if (weight != levels[index][0]
                or abs(energy - levels[index][1]) > .002
                or (index == 1 and (configuration != '2p6' or '2p6' not in label))
                or (index > 1 and (orbital is None or orbital.group(1) != configuration))
                or int(source[52:57].strip().split('.')[0]) * 2 + 1 != weight):
            raise ValueError(f'Ne I source level {index} identity changed')
    for index, energy, weight, orbital in ((3, 134459.287, 3, '3s'),
                                           (5, 135888.717, 3, '3s'),
                                           (15, 152970.733, 1, '3p')):
        if (levels[index][:2] != (weight, energy)
                or chianti_lines[index][7:37].strip() != orbital):
            raise ValueError('Ne I identified radiative state changed')

    stout_rates = {}
    for line in stout_paths[1].read_text().splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[0] == 'A' and fields[1].isdigit() and fields[2].isdigit():
            pair = int(fields[1]), int(fields[2])
            if pair in IDENTIFIED_RATES:
                if pair in stout_rates:
                    raise ValueError('Ne I duplicate identified Stout rate')
                stout_rates[pair] = float(fields[3])
    if stout_rates != IDENTIFIED_RATES:
        raise ValueError('Ne I identified Stout rates changed')
    source_pairs = {(lower, upper) for lower, upper, _, _, rate in radiation
                    if rate > 0}
    if len(source_pairs) != 142 or source_pairs & IDENTIFIED_RATES.keys():
        raise ValueError('Ne I original spontaneous pairs changed')
    corrected = list(radiation)
    for (lower, upper), rate in IDENTIFIED_RATES.items():
        gap = levels[upper][1] - levels[lower][1]
        if not math.isfinite(gap) or gap <= 0:
            raise ValueError('Ne I identified photon energy invalid')
        corrected.append((lower, upper, 1e8 / gap, 0.0, rate))
    return corrected
