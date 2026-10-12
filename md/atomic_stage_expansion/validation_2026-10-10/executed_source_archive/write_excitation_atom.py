"""Write a version-2 n-level atom with the identified source collision nodes.

Numerical use is 1,000--100,000 K. This is not a claim that the underlying
scattering calculation is physically accurate throughout that interval. Stout
Cl VII and Fe VII tables start at 9,800 K; below that, their adopted Upsilon
is held at its first node while Maxwellian rates use the actual temperature.
"""

import argparse
from hashlib import sha256
import io
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from atomic_output import atomic_write_text
from chianti_cooling import (BT_SOURCE_TMIN_K, bound_excitation_atom, ion_dir,
                             scaled_collision_strength_interpolator, upsilon)


USE_TMIN_K = 1.0e3
USE_TMAX_K = 1.0e5
BT_SOURCE_TMAX_K = {('c', 4): 3.2e7, ('mg', 3): 1.8e7, ('si', 5): 5.0e7}


def _source_files(element, stage):
    """Paths whose raw bytes determined the returned bound atom."""
    from ionization_potentials import SOURCE as ionization_source
    if (element, stage) in {('cl', 7), ('fe', 7), ('ar', 2), ('ca', 4)}:
        stem = Path.home() / f'CLOUDY/c25.00/data/stout/{element}/{element}_{stage}/{element}_{stage}'
        return [Path(f'{stem}.{ext}') for ext in ('nrg', 'tp', 'coll')] + [ionization_source]
    stem = Path(ion_dir(element, stage)) / f'{element}_{stage}'
    paths = [Path(f'{stem}.{ext}') for ext in ('elvlc', 'wgfa', 'scups')]
    if (element, stage) == ('ne', 1):
        from neon_spontaneous_probabilities import STOUT_STEM
        paths.extend(Path(f'{STOUT_STEM}.{ext}') for ext in ('nrg', 'tp'))
    if (element, stage) == ('n', 4):
        from niv_level_energies import STOUT_NRG
        paths.append(STOUT_NRG)
    elif (element, stage) == ('o', 5):
        from atomic_excitation_energies import STOUT_OV
        paths.append(STOUT_OV)
    elif (element, stage) in {('s', 7)}:
        paths.append(Path.home() / 'CLOUDY/c25.00/data/stout/s/s_7/s_7.nrg')
    paths.append(ionization_source)
    return paths


def _selection_code_files(element, stage):
    """Code that chooses bound identities, energies, and interpolation."""
    root = Path(__file__).resolve().parent
    names = ['write_excitation_atom.py', 'chianti_cooling.py', 'atomic_excitation_energies.py',
             'ionization_potentials.py']
    if (element, stage) in {('o', 4), ('ne', 5), ('s', 6), ('s', 7),
                            ('ar', 5), ('mg', 4)}:
        names.append('valence_excitation_energies.py')
    if (element, stage) == ('n', 4):
        names.append('niv_level_energies.py')
    if (element, stage) in {('o', 2), ('o', 3), ('mg', 3)}:
        names.append('electric_dipole_probabilities.py')
    if (element, stage) == ('ne', 1):
        names.append('neon_spontaneous_probabilities.py')
    if (element, stage) == ('cl', 7):
        names.append('stout_excitation_cooling.py')
    if (element, stage) == ('fe', 7):
        names.append('stout_excitation_atoms.py')
    if (element, stage) in {('ar', 2), ('ca', 4)}:
        names.append('stout_fine_structure_atoms.py')
    return [root / name for name in names]


def _number(value):
    """Round-trip an IEEE binary64 value in a portable decimal record."""
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Nonfinite atomic source value')
    return format(value, '.17g')


def _energy_pair(levels, first, second, *, collision=False):
    """Orient one transition by the selected physical level energies."""
    if first == second or first not in levels or second not in levels:
        raise ValueError('Invalid source transition levels')
    energy_first, energy_second = levels[first][1], levels[second][1]
    if energy_first == energy_second:
        # Degenerate states still exchange electrons. Detailed balance uses
        # their statistical weights and a unit Boltzmann factor; the positive
        # BT scaling energy is independent of this zero physical level gap.
        if collision:
            return min(first, second), max(first, second)
        raise ValueError(f'Zero-energy transition {first}-{second}')
    if not collision and energy_first > energy_second:
        raise ValueError(f'Radiative transition {first}->{second} reverses physical energy order; '
                         'identify its multipole and transform A before writing')
    return (first, second) if energy_first < energy_second else (second, first)


def _collision_row(levels, transition, index):
    first, second = int(transition['ll']), int(transition['ul'])
    lower, upper = _energy_pair(levels, first, second, collision=True)
    if 'temperature_grid' in transition:
        x = np.asarray(transition['temperature_grid'], dtype=float)
        y = np.asarray(transition['upsilon_grid'], dtype=float)
        if x.size < 2 or x.shape != y.shape or np.any(y <= 0):
            raise ValueError(f'Invalid temperature table for collision {index}')
        model, kind, ttype, cups, de = 'TEMPERATURE_TABLE', 'LOG_PCHIP_T', 0, 0., 0.
        if x[0] > USE_TMAX_K or x[-1] < USE_TMAX_K:
            raise ValueError(f'Temperature source does not cover numerical use: {index}')
        # The shared evaluator supplies the adopted constant cold continuation.
        samples = np.geomspace(USE_TMIN_K, USE_TMAX_K, 33)
        values = upsilon(transition, samples)
        check = np.exp(PchipInterpolator(np.log(x), np.log(y))(
            np.log(np.maximum(samples, x[0]))))
        if not np.allclose(values, check, rtol=1e-12, atol=0):
            raise ValueError(f'Temperature-table model differs from source: {index}')
    else:
        x = np.asarray(transition['xs'], dtype=float)
        y = np.asarray(transition['ys'], dtype=float)
        ttype = int(transition['ttype'])
        cups = float(transition['cups'])
        de = float(transition['de'])
        if ttype not in range(1, 7) or cups <= 0 or de <= 0:
            raise ValueError(f'Invalid Burgess-Tully parameters: {index}')
        if x.size < 2 or x.shape != y.shape:
            raise ValueError(f'Invalid scaled nodes: {index}')
        if ttype != 6 and np.any(y < 0):
            raise ValueError(f'Negative nonlogarithmic source nodes unsupported: {index}')
        interpolator = scaled_collision_strength_interpolator(transition)
        kind = 'PCHIP' if isinstance(interpolator, PchipInterpolator) else 'CUBIC'
        model = 'BT_SOURCE'
        # These source nodes, rather than a fitted polynomial, are serialized.
        samples = np.geomspace(USE_TMIN_K, USE_TMAX_K, 33)
        values = upsilon(transition, samples)
    if (not np.all(np.isfinite(x)) or not np.all(np.isfinite(y))
        or np.any(np.diff(x) <= 0) or not np.all(np.isfinite(values))
        or np.any(values < 0)):
        raise ValueError(f'Invalid source collision curve: {index}')
    return lower, upper, model, ttype, cups, de, kind, x, y


def write_excitation_atom(element, stage, levels=None, output_dir=None):
    """Write one deterministic source-node atom and return structural metadata.

    ``levels`` selects the first N collision-connected bound source indices;
    omitting it retains the full bound set. The returned ``passed`` flag means
    structural generation completed, not that compiled native rates passed.
    """
    element = element.lower()
    if not element.isalpha() or not isinstance(stage, int) or stage < 1:
        raise ValueError('Specify an element symbol and positive ion stage')
    if levels is not None and (not isinstance(levels, int) or levels < 1):
        raise ValueError('levels must be a positive integer')
    if output_dir is None:
        raise ValueError('An explicit output directory is required')
    bound, radiation, collisions = bound_excitation_atom(element, stage)
    source_ids = sorted(bound)
    if levels is not None and levels > len(source_ids):
        raise ValueError(f'Requested {levels} levels; only {len(source_ids)} bound levels exist')
    source_ids = source_ids[:levels]
    selected = {original: bound[original] for original in source_ids}
    native_index = {original: compact for compact, original in enumerate(source_ids, 1)}
    if not selected or any(g <= 0 or not math.isfinite(g) or not math.isfinite(e)
                           for g, e in selected.values()):
        raise ValueError('Invalid selected levels')
    if min(e for _, e in selected.values()) != 0:
        raise ValueError('Selected atom lacks a zero-energy ground state')

    radiative = {}
    for first, second, _, _, rate in radiation:
        if first not in selected or second not in selected:
            continue
        if not math.isfinite(rate) or rate < 0:
            raise ValueError(f'Invalid radiative rate {first}-{second}')
        if rate == 0:
            continue
        lower, upper = _energy_pair(selected, first, second)
        pair = (native_index[lower], native_index[upper])
        radiative[pair] = radiative.get(pair, 0.) + float(rate)
    if not all(math.isfinite(rate) and rate > 0 for rate in radiative.values()):
        raise ValueError('Invalid summed radiative rate')

    rows = []
    for index, tr in enumerate(collisions, 1):
        if tr['ll'] not in selected or tr['ul'] not in selected:
            continue
        lower, upper, *rest = _collision_row(selected, tr, index)
        rows.append((native_index[lower], native_index[upper], *rest))
    if not rows:
        raise ValueError('Selected atom has no collision transitions')
    rows.sort(key=lambda row: (row[0], row[1]))
    pairs = [(row[0], row[1]) for row in rows]
    if len(set(pairs)) != len(pairs):
        raise ValueError('Duplicate collision pairs in source atom')
    total_nodes = sum(len(row[-1]) for row in rows)
    bt_tmin = BT_SOURCE_TMIN_K.get((element, stage))
    if bt_tmin is not None and any(
            tr.get('bt_tmin_K') != bt_tmin or 'temperature_grid' in tr
            for tr in collisions if tr['ll'] in selected and tr['ul'] in selected):
        raise ValueError('Inconsistent BT source temperature minimum')

    output = io.StringIO()
    output.write(f'# MoCHII source-node atom: {element}_{stage}\n')
    output.write(f'# numerical use: {_number(USE_TMIN_K)}..{_number(USE_TMAX_K)} K; source validity may differ\n')
    if (element, stage) == ('mg', 3):
        output.write('# Mg III source ID 167 uses theoretical 634430 cm^-1; its J=0 label rejects the observed energy assigned to a J=1 7d state\n')
        output.write('# Higher E1 A values preserve dipole strength with physical frequency cubed and emitting statistical weight; first 27 A values are unchanged\n')
        output.write('# Higher gf=0 A values are approximate effective rates with unidentified multipoles; wavefunction and collision accuracy remain unverified\n')
    if (element, stage) in {('cl', 7), ('fe', 7)}:
        output.write('# Stout temperature nodes: 9800..9.8e7 K; below 9800 K hold first Upsilon constant; Maxwellian rates use actual T\n')
        output.write('# cold resonance accuracy below 9800 K is unverified\n')
    if (element, stage) in {('ar', 2), ('ca', 4)}:
        temperatures = collisions[0]['temperature_grid']
        output.write(f'# Pelan and Berrington (1995), Table 3: source temperature nodes {_number(temperatures[0])}..{_number(temperatures[-1])} K\n')
        output.write('# Adopted constant Upsilon below the first node; no continuation above the last node\n')
        output.write('# Ground fine-structure electron pair only; higher excited-state cooling is not established\n')
    if (element, stage) == ('ne', 1):
        output.write('# radiation: 142 original CHIANTI A rows plus NIST-dated Stout Ne I A(3,15)=900000 and A(5,15)=68200000 s^-1\n')
        output.write('# Ne I CHIANTI BT first finite nodes infer 999.410754..1000.355160 K from rounded parameters; not a published physical temperature range\n')
        output.write('# adopted numerical Upsilon floor: 1000.36 K; below it hold Upsilon constant while Maxwellian rates use actual T\n')
        output.write('# numerical atomic-data approximation: 1000..100000 K; source conversion and atomic accuracy unverified\n')
    elif bt_tmin is not None:
        output.write(f'# Published BT collision temperature interval: {_number(bt_tmin)}..'
                     f'{_number(BT_SOURCE_TMAX_K[(element, stage)])} K\n')
        output.write('# Adopt constant Upsilon below the published lower bound\n')
        output.write('# Maxwellian rates and detailed balance use actual gas temperature; '
                     'cold continuation is an approximation\n')
    output.write('# source level IDs: ' + ' '.join(map(str, source_ids)) + '\n')
    for path in _source_files(element, stage):
        output.write(f'# source SHA-256 {path}: {sha256(path.read_bytes()).hexdigest()}\n')
    for path in _selection_code_files(element, stage):
        output.write(f'# selection code SHA-256 {path}: {sha256(path.read_bytes()).hexdigest()}\n')
    output.write('# NODES is the total number of (x, y) records across all collision rows\n')
    output.write('NLEVEL_VERSION 2\n')
    output.write(f'NLEV {len(selected)}\n')
    for original, (g, energy) in selected.items():
        output.write(f'{native_index[original]} {_number(g)} {_number(energy)}\n')
    output.write(f'NRAD {len(radiative)}\n')
    for (lower, upper), rate in sorted(radiative.items()):
        output.write(f'{lower} {upper} {_number(rate)}\n')
    output.write(f'NODES {total_nodes}\n')
    output.write(f'NUPS {len(rows)}\n')
    for lower, upper, model, ttype, cups, de, kind, x, y in rows:
        output.write(f'{lower} {upper} {model} {ttype} {_number(cups)} {_number(de)} {kind} {len(x)}\n')
        for abscissa, ordinate in zip(x, y):
            output.write(f'{_number(abscissa)} {_number(ordinate)}\n')
    if bt_tmin is not None:
        output.write(f'BT_TMIN_K {_number(bt_tmin)}\n')
    # A complete bound atom has no omitted excitation cooling. Partial atoms
    # may use the same density ratio only after independent convergence checks.
    complete = len(selected) == len(bound)
    if complete:
        output.write('# Full collision-connected bound source: use its density ratio for fitted cooling\n')
        output.write('DENSITY_RATIO 1\n')
    path = Path(output_dir) / f'nlevel_{element}_{stage}.txt'
    atomic_write_text(path, output.getvalue())
    return {'passed': True, 'path': str(path.resolve()), 'element': element,
            'stage': stage, 'nlev': len(selected), 'nrad': len(radiative),
            'nups': len(rows), 'nodes': total_nodes,
            'selected_ids': source_ids, 'density_ratio_cooling': complete,
            'bt_tmin_K': bt_tmin}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--element', required=True)
    parser.add_argument('--stage', type=int, required=True)
    parser.add_argument('--levels', type=int)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args(argv)
    result = write_excitation_atom(args.element, args.stage, args.levels,
                                   args.output_dir)
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
