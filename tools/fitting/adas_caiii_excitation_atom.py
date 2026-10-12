"""Private Ca III excitation candidate from the OPEN-ADAS type-3 ADF04 file.

Source: adf04/cophps#ar/dw/ic#ca2.dat (Alessandra Giunta, 2012).
The source is an Autostructure distorted-wave calculation. It does not resolve
near-threshold close-coupling resonances, so the printed temperature interval
alone does not establish accuracy at low temperature. OPEN-ADAS limits use of
downloaded files and restricts redistribution with modeling code. This reader
keeps the original outside the repository and does not select a production atom.

ADF04 type 3 is defined in https://www.adas.ac.uk/man/appxa-04.pdf: the header
temperatures are kelvin; a transition row stores upper, lower, Einstein A in
s^-1, one Maxwellian effective collision strength for each temperature, and
an optional Bethe high-energy limit. For intermediate-coupling levels the
parenthesized quantum number is J, giving statistical weight 2*J+1.
An independent OPEN-ADAS Ne I ADF04 file labels A=1.0E-30 as an unphysical
placeholder. The Ca III source does not explain its identical values, so
omitting them from radiative rates is a provisional physical interpretation.
"""

from hashlib import sha256
from pathlib import Path
import math
import re

import numpy as np


SOURCE_SHA256 = '8d9d89220687b518899f98893ca4e53834e1990179cfc4816cf5bf3a9f02c5d7'
SOURCE_URL = ('https://open.adas.ac.uk/element/detail/adf04/'
              'cophps%5D%5Bar/dw/ic%5D%5Bca2.dat')
IONIZATION_CM1 = 410668.5947
BOUND_LEVEL_COUNT = 39
TEMPERATURE_COUNT = 14
ALL_LEVEL_COUNT = 180
ALL_COLLISION_COUNT = 15639
ABSENT_A_SENTINEL = 1e-30
BOUND_ABSENT_A_COUNT = 78
CM1_TO_RY = 1.0 / 109737.31568539

_LEVEL = re.compile(
    r'^\s*(\d+)\s+(\S+)\s+\((\d+)\)(\d+)\(\s*(\d+(?:\.\d+)?)\)'
    r'\s+(\d+(?:\.\d+)?)\s+\{X\}\s*$')
_ROW = re.compile(r'^\s*(\d+)\s+(\d+)\s+')
_COMPACT_REAL = re.compile(r'([+-]?(?:\d+\.\d*|\.\d+))([+-]\d+)')


def _compact_values(text, count):
    """Read adjacent Fortran E-format fields without an exponent letter."""
    matches = list(_COMPACT_REAL.finditer(text))
    if len(matches) != count:
        raise ValueError(f'Expected {count} compact real fields, found {len(matches)}')
    if _COMPACT_REAL.sub('', text).strip():
        raise ValueError('Unparsed ADF04 compact real field')
    values = [float(match.group(1) + 'e' + match.group(2)) for match in matches]
    if not all(math.isfinite(value) for value in values):
        raise ValueError('Nonfinite ADF04 value')
    return values


def _parse_caiii_adf04(raw):
    """Validate this source's complete file and return its bound-state triplet."""
    lines = raw.decode('ascii').splitlines()
    if not lines:
        raise ValueError('Empty Ca III ADF04 source')
    header = re.fullmatch(
        r'Ca\+\s*2\s+20\s+3\s+(\d+(?:\.\d+)?)\(1s\)',
        lines[0].strip(), flags=re.IGNORECASE)
    if header is None:
        raise ValueError('Unexpected Ca III ion or ionization header')
    threshold = float(header.group(1))
    if threshold != IONIZATION_CM1:
        raise ValueError('Ca III source ionization boundary changed')

    levels_all = {}
    cursor = 1
    while cursor < len(lines) and not lines[cursor].lstrip().startswith('-1'):
        match = _LEVEL.fullmatch(lines[cursor])
        if match is None:
            raise ValueError(f'Invalid ADF04 level record at line {cursor+1}')
        index = int(match.group(1))
        multiplicity, orbital_l = int(match.group(3)), int(match.group(4))
        j, energy = float(match.group(5)), float(match.group(6))
        spin = (multiplicity - 1) / 2
        if (index != len(levels_all) + 1 or multiplicity < 1 or orbital_l < 0
            or not math.isfinite(energy) or energy < 0
            or j < abs(orbital_l - spin) or j > orbital_l + spin
            or not math.isclose(j - abs(orbital_l-spin),
                                round(j - abs(orbital_l-spin)), abs_tol=1e-9)
            or (index == 1 and energy != 0)
            or (index > 1 and energy <= levels_all[index-1][1])):
            raise ValueError(f'Invalid Ca III level index, quantum number, or energy {index}')
        levels_all[index] = (2*j + 1, energy, energy)
        cursor += 1
    if len(levels_all) != ALL_LEVEL_COUNT or cursor >= len(lines):
        raise ValueError('Incomplete Ca III fine-structure level list')
    if _ROW.match(lines[cursor]):
        raise ValueError('Missing ADF04 level terminator')
    cursor += 1
    if cursor >= len(lines):
        raise ValueError('Missing ADF04 temperature header')
    temperature_header = lines[cursor].split(maxsplit=2)
    if len(temperature_header) != 3 or temperature_header[:2] != ['3.00', '3']:
        raise ValueError('Expected type-3 ADF04 Kelvin Upsilon grid')
    temperatures = np.asarray(_compact_values(temperature_header[2], TEMPERATURE_COUNT))
    if (np.any(temperatures <= 0) or np.any(np.diff(temperatures) <= 0)
        or temperatures[0] != 900 or temperatures[-1] != 18_000_000):
        raise ValueError('Invalid Ca III electron temperature grid')
    cursor += 1
    bound_indices = [i for i, (_, energy, _) in levels_all.items() if energy < threshold]
    if bound_indices != list(range(1, BOUND_LEVEL_COUNT+1)):
        raise ValueError('Ca III bound-state set changed')
    levels = {index: levels_all[index] for index in bound_indices}

    from chianti_cooling import KB_OVER_RY

    radiation = []
    collisions = []
    seen = set()
    absent_a_count = 0
    while cursor < len(lines) and lines[cursor].strip() != '-1':
        line = lines[cursor]
        match = _ROW.match(line)
        if match is None:
            raise ValueError(f'Invalid ADF04 electron transition at line {cursor+1}')
        upper, lower = int(match.group(1)), int(match.group(2))
        if (not 1 <= lower < upper <= ALL_LEVEL_COUNT
            or (upper, lower) in seen):
            raise ValueError(f'Duplicate or reversed Ca III electron pair {lower}->{upper}')
        seen.add((upper, lower))
        values = _compact_values(line[match.end():], TEMPERATURE_COUNT + 2)
        a_value, strengths, bethe_limit = values[0], values[1:-1], values[-1]
        if a_value < 0 or any(value <= 0 for value in strengths):
            raise ValueError(f'Nonphysical Ca III A or Upsilon {lower}->{upper}')
        if upper <= BOUND_LEVEL_COUNT:
            gap = levels[upper][1] - levels[lower][1]
            de = gap * CM1_TO_RY
            et = KB_OVER_RY * temperatures / de
            # A separate ADAS file identifies this value as unphysical.
            # Omission here is provisional because this Ca III file has no note.
            if a_value == ABSENT_A_SENTINEL:
                absent_a_count += 1
            else:
                radiation.append((lower, upper, 1e8/gap, 0., a_value))
            collisions.append(dict(
                ll=lower, ul=upper, de=de, gf=0., lim=bethe_limit,
                ttype=2, cups=1., xs=et/(1+et), ys=np.asarray(strengths),
                temperature_grid=temperatures, upsilon_grid=np.asarray(strengths)))
        cursor += 1
    if (len(seen) != ALL_COLLISION_COUNT or len(collisions) !=
        BOUND_LEVEL_COUNT*(BOUND_LEVEL_COUNT-1)//2 or cursor >= len(lines)
        or absent_a_count != BOUND_ABSENT_A_COUNT):
        raise ValueError('Incomplete Ca III electron-pair set')
    if (cursor+1 >= len(lines) or lines[cursor].strip() != '-1'
        or lines[cursor+1].split() != ['-1', '-1']):
        raise ValueError('Invalid ADF04 transition terminator')
    return levels, radiation, collisions


def adas_caiii_excitation_atom(path):
    """Return checked Ca III (levels, A rows, electron rows) from a private file.

    The caller supplies the locally downloaded original. No source rows are
    included in this repository. This routine is an investigation candidate;
    source licensing and distorted-wave accuracy require separate decisions.
    """
    raw = Path(path).read_bytes()
    if sha256(raw).hexdigest() != SOURCE_SHA256:
        raise ValueError('Changed OPEN-ADAS Ca III source')
    return _parse_caiii_adf04(raw)
