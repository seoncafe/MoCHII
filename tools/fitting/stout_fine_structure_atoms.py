"""Published Ar II and Ca IV ground-term electron excitation source atoms.

Pelan and Berrington (1995, A&AS 110, 209), Table 3 tabulates the sole
chlorine-like ground fine-structure collision pair. Their R-matrix calculation
uses LS coupling followed by algebraic recoupling. The table covers log10(T/K)
3.2--5.2 for Ar II and 4.0--6.0 for Ca IV. The paper excludes Ar II from its
better-than-10-percent accuracy estimate. Constant Upsilon below the first
node is an adopted approximation for 1,000--100,000 K numerical use, not part
of the published validity interval. No extrapolation above the final node is
allowed. Only the ground pair has identified electron collision data here.
"""

from hashlib import sha256
from pathlib import Path

import numpy as np


STOUT_ROOT = Path.home() / 'CLOUDY/c25.00/data/stout'
SOURCES = {
    ('ar', 2): {
        'sha256': (
            '2f8151e6d5b081488478f3252f37134a72cfafb534a846d205575c1d0571ff30',
            '0dc0d90620687c7f67e1b9369b38158d1d1778956796ae335812951d8cd5f655',
            'c6f860bc6d8e2f88c5a5e43e487fbf440c56431a593bc9feca22ba0828624e0e'),
        'energy_cm1': 1431.583,
        'a_s1': .0530,
        'source_log_t': np.arange(3.2, 5.21, .2),
        'published_upsilon': (2.48, 2.54, 2.63, 2.77, 2.93, 3.09,
                              3.19, 3.20, 3.13, 2.97, 2.71),
    },
    ('ca', 4): {
        'sha256': (
            '333ea9ead15655b240d573036201aea288c0b4f898d05b96599f816bb7aefaaf',
            'fc71397bc3cc9b9225c180a987143fc5ec82bdf1a77d630e0563032dfea77e32',
            'e68cce06232964d061e6ef945b353ad7aaaff68ea5522a10aee767b33828587f'),
        'energy_cm1': 3118.200,
        'a_s1': .543,
        'source_log_t': np.arange(4.0, 6.01, .2),
        'published_upsilon': (1.00, 1.12, 1.35, 1.77, 2.39, 2.80,
                              2.98, 2.82, 2.43, 1.97, 1.53),
    },
}


def _records(path, expected_hash):
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != expected_hash:
        raise ValueError(f'Changed Stout source: {path}')
    lines = raw.decode('ascii').splitlines()
    if not lines or lines[0].strip() != '17 09 05':
        raise ValueError(f'Unexpected Stout header: {path}')
    return [line.split() for line in lines[1:]
            if line.strip() and not line.startswith(('*', '#'))]


def stout_fine_structure_atom(element, stage, root=STOUT_ROOT):
    """Return checked two-level (levels, A rows, electron rows) source data."""
    element = element.lower()
    key = (element, stage)
    if key not in SOURCES:
        raise ValueError(f'No identified two-level source for {element}_{stage}')
    specification = SOURCES[key]
    stem = root / element / f'{element}_{stage}' / f'{element}_{stage}'
    nrg, tp, coll = (_records(Path(f'{stem}.{extension}'), digest)
                     for extension, digest in zip(('nrg', 'tp', 'coll'),
                                                  specification['sha256']))
    energies = {}
    for fields in nrg:
        if fields[0].isdigit() and int(fields[0]) in (1, 2):
            index = int(fields[0])
            if index in energies or len(fields) < 3:
                raise ValueError('Duplicate or incomplete fine-structure level')
            energy, weight = float(fields[1]), float(fields[2])
            energies[index] = (weight, energy, energy)
    if (set(energies) != {1, 2} or energies[1] != (4., 0., 0.)
        or energies[2] != (2., specification['energy_cm1'],
                            specification['energy_cm1'])):
        raise ValueError('Stout fine-structure energies changed')

    a_rows = [fields for fields in tp if fields[:3] == ['A', '1', '2']]
    if len(a_rows) != 1 or float(a_rows[0][3]) != specification['a_s1']:
        raise ValueError('Stout ground-term A value changed')
    rate = float(a_rows[0][3])

    temperature_rows = [fields for fields in coll if fields[0] == 'TEMP']
    collision_rows = [fields for fields in coll if fields[:2] == ['CS', 'ELECTRON']]
    if (len(temperature_rows) != 1 or len(collision_rows) != 1
        or collision_rows[0][2:4] != ['1', '2']):
        raise ValueError('Incomplete fine-structure electron source')
    temperatures = np.asarray([float(value) for value in temperature_rows[0][1:]])
    strengths = np.asarray([float(value) for value in collision_rows[0][4:]])
    if (len(temperatures) != 11 or strengths.shape != temperatures.shape
        or np.any(~np.isfinite(temperatures)) or np.any(~np.isfinite(strengths))
        or np.any(np.diff(temperatures) <= 0) or np.any(strengths <= 0)
        or not np.array_equal(strengths, specification['published_upsilon'])
        or not np.allclose(np.log10(temperatures), specification['source_log_t'],
                           atol=.0015, rtol=0)):
        raise ValueError('Stout Table 3 temperature or Upsilon nodes changed')

    # The physical bound-state cutoff comes from the identified NIST snapshot.
    from ionization_potentials import ELEMENTS, read_ionization_potentials, HC_EV_CM
    z = dict(ELEMENTS)[element]
    edge_cm1 = read_ionization_potentials()[z, z-stage+1] / HC_EV_CM
    if energies[2][1] >= edge_cm1:
        raise ValueError('Fine-structure upper state is unbound')

    from chianti_cooling import HC_ERG_CM, RY_ERG, KB_OVER_RY
    gap = energies[2][1] - energies[1][1]
    de = gap * HC_ERG_CM / RY_ERG
    et = KB_OVER_RY * temperatures / de
    transition = dict(ll=1, ul=2, de=de, gf=0., lim=-1., ttype=2,
                      cups=1., xs=et/(1.+et), ys=strengths,
                      temperature_grid=temperatures, upsilon_grid=strengths)
    radiation = [(1, 2, 1e8/gap, 0., rate)]
    return energies, radiation, [transition]
