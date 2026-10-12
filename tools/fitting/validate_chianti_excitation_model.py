"""Checks of the 'chianti' excitation model atoms against their sources.

For each nlevel_<el>_<stage>.txt in the model directory:

1. Radiative rows: every photon energy E_upper - E_lower from the file levels
   is positive and every A is positive.
2. Detailed balance of the collision rates: the Fortran loader's collision
   matrix (nlevel_cooling_prepT) solved with the radiative rates set to zero
   (nlevel_populations) must give the Boltzmann populations of the file's
   levels; reported at 3e3, 1e4, 1e5 K as the largest relative deviation over
   levels with Boltzmann fraction > 1e-6 and the largest absolute deviation.
3. Serialization: Upsilon from the file nodes (evaluated as nlevel_mod does:
   Burgess-Tully descaling of the stored nodes, not-a-knot cubic or PCHIP,
   BT_TMIN_K hold, log-PCHIP temperature tables) against
   chianti_cooling.upsilon on the source transitions (mapped through the
   header's source level IDs).
4. Low-density cooling: the Fortran loader and solver (nlevel_mod, compiled
   from a snapshot of src by the caller, driver excitation_atom_cooling.x) at
   n_e = 1e-6 cm^-3 against chianti_cooling.cooling_low_density_limit on the
   CHIANTI (or Stout) source directly, 81 temperatures over 1e3-1e5 K;
   temperatures where the source cooling underflows (< 1e-250) are counted,
   not compared.
5. Cooling fit: cooling_<el>_<stage>.txt against the same source cooling on
   801 temperatures over 1e3-1e5 K, and its term count.

Run: python3 validate_chianti_excitation_model.py --model-dir DIR --tag o_3
     --fortran-cooling FILE --output DIR
"""
import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

import chianti_cooling as C

KB_OVER_RY = C.K_B_ERG / C.RY_ERG
HC = C.HC_ERG_CM


def parse_atom(path):
    """Read a NLEVEL_VERSION 2 source-node atom as nlevel_mod does."""
    header, records = [], []
    for line in Path(path).read_text().splitlines():
        if line.lstrip().startswith('#'):
            header.append(line)
        elif line.strip():
            records.append(line.split())
    atom = dict(header=header, bt_tmin=None, density_ratio=None)
    i = 0
    while i < len(records):
        key = records[i][0]
        if key == 'NLEVEL_VERSION':
            atom['version'] = int(records[i][1]); i += 1
        elif key == 'NLEV':
            n = int(records[i][1])
            atom['g'] = np.array([float(r[1]) for r in records[i+1:i+1+n]])
            atom['E'] = np.array([float(r[2]) for r in records[i+1:i+1+n]])
            if [int(r[0]) for r in records[i+1:i+1+n]] != list(range(1, n+1)):
                raise ValueError('Level indices are not 1..NLEV')
            i += 1 + n
        elif key == 'NRAD':
            n = int(records[i][1])
            atom['rad'] = [(int(r[0]), int(r[1]), float(r[2])) for r in records[i+1:i+1+n]]
            i += 1 + n
        elif key == 'NODES':
            atom['nodes'] = int(records[i][1]); i += 1
        elif key == 'NUPS':
            n = int(records[i][1]); i += 1
            rows = []
            for _ in range(n):
                r = records[i]
                count = int(r[7])
                xy = np.array([[float(a), float(b)] for a, b in records[i+1:i+1+count]])
                rows.append(dict(l=int(r[0]), u=int(r[1]), model=r[2], ttype=int(r[3]),
                                 cups=float(r[4]), de=float(r[5]), kind=r[6],
                                 x=xy[:, 0], y=xy[:, 1]))
                i += 1 + count
            atom['coll'] = rows
        elif key == 'BT_TMIN_K':
            atom['bt_tmin'] = float(records[i][1]); i += 1
        elif key == 'DENSITY_RATIO':
            atom['density_ratio'] = int(records[i][1]); i += 1
        else:
            raise ValueError(f'Unknown section {key}')
    ids = next((h for h in header if h.startswith('# source level IDs:')), None)
    atom['source_ids'] = [int(v) for v in ids.split(':', 1)[1].split()] if ids else None
    return atom


def file_upsilon(atom, row, T):
    """Effective collision strength from the stored nodes, as nlevel_mod."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    if row['model'] == 'TEMPERATURE_TABLE':
        x, y = np.log(row['x']), np.log(row['y'])
        s = np.minimum(np.log(T), x[-1])
        value = np.exp(PchipInterpolator(x, y)(np.maximum(s, x[0])))
        return np.where(s <= x[0], row['y'][0], value)
    Ts = np.maximum(T, atom['bt_tmin']) if atom['bt_tmin'] else T
    et = KB_OVER_RY * Ts / row['de']
    c = row['cups']
    st = 1 - np.log(c) / np.log(et + c) if row['ttype'] in (1, 4) else et / (et + c)
    st = np.clip(st, row['x'][0], row['x'][-1])
    spline = (CubicSpline(row['x'], row['y']) if row['kind'] == 'CUBIC'
              else PchipInterpolator(row['x'], row['y']))
    y = spline(st)
    t = row['ttype']
    return {1: lambda: y * np.log(et + np.e), 2: lambda: y, 3: lambda: y / (et + 1),
            4: lambda: y * np.log(et + c), 5: lambda: y / et, 6: lambda: 10.0 ** y}[t]()


def collision_balance_deviation(atom, path):
    """Deviation of collision-only Fortran populations from Boltzmann."""
    table = np.atleast_2d(np.loadtxt(path))
    g, E = atom['g'], atom['E']
    result = {}
    for row in table:
        T, population = row[0], row[1:]
        log_weight = np.log(g) - (E - E.min())*HC/(C.K_B_ERG*T)
        boltzmann = np.exp(log_weight - np.logaddexp.reduce(log_weight))
        large = boltzmann > 1e-6
        result[f'{T:.0e}'] = {
            'max_relative_where_above_1e-6': float(np.max(np.abs(population[large]/boltzmann[large]-1))),
            'max_absolute': float(np.max(np.abs(population-boltzmann)))}
    return result


def serialized_upsilon_difference(atom, element, stage):
    """Max relative difference of file and source Upsilon at five temperatures."""
    if atom['source_ids'] is None:
        return None
    levels, _, collisions = C.bound_excitation_atom(element, stage)
    native = {s: k for k, s in enumerate(atom['source_ids'], 1)}
    source = {}
    for tr in collisions:
        if tr['ll'] in native and tr['ul'] in native:
            pair = tuple(sorted((native[tr['ll']], native[tr['ul']])))
            source[pair] = tr
    temperatures = np.array([1e3, 3e3, 1e4, 3e4, 1e5])
    worst = 0.0
    if len(source) != len(atom['coll']):
        raise ValueError('Collision rows of file and source differ in number')
    for row in atom['coll']:
        tr = source[tuple(sorted((row['l'], row['u'])))]
        a = file_upsilon(atom, row, temperatures)
        b = C.upsilon(tr, temperatures)
        scale = np.maximum(np.abs(b), 1e-300)
        worst = max(worst, float(np.max(np.abs(a - b) / scale)))
    return worst


def read_cooling_fit(path):
    records = [l.split() for l in Path(path).read_text().splitlines()
               if l.strip() and not l.lstrip().startswith('#')]
    n = int(records[0][0])
    data = np.array([[float(a), float(b)] for a, b in records[1:1+n]])
    return data[:, 0], data[:, 1]


def evaluate_fit(A, Ti, T):
    return np.sum(A[:, None] * np.exp(-Ti[:, None] / T[None, :]), axis=0) / np.sqrt(T)


def relative_error(fit, source):
    positive = source > 0
    error = np.abs(fit[positive] / source[positive] - 1)
    return float(np.max(error)) if error.size else None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-dir', type=Path, required=True)
    parser.add_argument('--tag', required=True, help='element_stage, e.g. o_3')
    parser.add_argument('--fortran-cooling', type=Path,
                        help='excitation_atom_cooling.x output for this atom')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    start = time.time()
    element, stage = args.tag.split('_')
    stage = int(stage)
    path = args.model_dir / f'nlevel_{args.tag}.txt'
    atom = parse_atom(path)
    result = dict(tag=args.tag, file=str(path), file_bytes=path.stat().st_size,
                  levels=len(atom['g']), radiative_rows=len(atom['rad']),
                  collision_rows=len(atom['coll']), nodes=atom['nodes'],
                  density_ratio=atom['density_ratio'], bt_tmin_K=atom['bt_tmin'])
    gaps = [(atom['E'][u-1] - atom['E'][l-1], a) for l, u, a in atom['rad']]
    result['radiative_nonpositive_gaps'] = sum(1 for gap, _ in gaps if gap <= 0)
    result['radiative_nonpositive_A'] = sum(1 for _, a in gaps if a <= 0)
    result['radiative_min_gap_cm'] = min(gap for gap, _ in gaps) if gaps else None
    if args.fortran_cooling is not None:
        balance = Path(str(args.fortran_cooling)+'.collision_balance')
        if balance.exists():
            result['collision_only_populations_vs_boltzmann'] = collision_balance_deviation(atom, balance)
    result['upsilon_file_vs_source_max_relative'] = serialized_upsilon_difference(
        atom, element, stage)
    if args.fortran_cooling is not None and args.fortran_cooling.exists():
        table = np.loadtxt(args.fortran_cooling)
        T81 = table[:, 0]
        source81 = C.cooling_low_density_limit(element, stage, T81)
        normal = source81 > 1e-250
        result['source_cooling_underflow_temperatures'] = int(np.sum(~normal))
        result['fortran_low_density_vs_source_max_relative'] = relative_error(
            table[normal, 1], source81[normal])
        positive = source81 > 1e-3 * source81.max()
        result['fortran_low_density_vs_source_max_relative_where_above_1e-3_max'] = (
            float(np.max(np.abs(table[positive, 1] / source81[positive] - 1))))
        result['fortran_cooling_at_1e4K'] = dict(zip(
            ['ne_1e-6', 'ne_1e2', 'ne_1e4', 'ne_1e6'],
            [float(np.interp(4.0, np.log10(T81), table[:, k])) for k in range(1, 5)]))
    fit_path = args.model_dir / f'cooling_{args.tag}.txt'
    if fit_path.exists():
        A, Ti = read_cooling_fit(fit_path)
        T801 = np.logspace(3, 5, 801)
        source801 = C.cooling_low_density_limit(element, stage, T801)
        fit801 = evaluate_fit(A, Ti, T801)
        result['fit_terms'] = int(len(A))
        result['fit_signed_terms'] = bool(np.any(A < 0))
        normal = source801 > 1e-250
        result['fit_vs_source_max_relative'] = relative_error(fit801[normal], source801[normal])
        mask = source801 > 1e-3 * source801.max()
        result['fit_vs_source_max_relative_where_above_1e-3_max'] = float(
            np.max(np.abs(fit801[mask] / source801[mask] - 1)))
        result['fit_negative_points'] = int(np.sum(fit801 < 0))
        result['fit_zero_points_where_source_below_1e-250'] = int(np.sum((fit801 == 0) & (source801 < 1e-250)))
        result['fit_zero_points_where_source_above_1e-250'] = int(np.sum((fit801 == 0) & (source801 >= 1e-250)))
    if args.tag == 'fe_2':
        decays = {k: sum(a for l, u, a in atom['rad'] if u == k) for k in range(14, 36)}
        result['feii_levels_14_35'] = dict(
            rows=sum(1 for l, u, a in atom['rad'] if 14 <= u <= 35),
            levels_without_decay=[k for k, a in decays.items() if a <= 0],
            min_total_A=min(decays.values()), max_total_A=max(decays.values()),
            total_A=decays)
    result['seconds'] = time.time() - start
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / f'{args.tag}.json').write_text(json.dumps(result, indent=1) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'feii_levels_14_35'}))


if __name__ == '__main__':
    main()
