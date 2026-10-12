"""Low-density Fe II line cooling fit for the 'mocassin_v2' excitation model.

Reads levels, spontaneous probabilities and the log-PCHIP effective-strength
tables from the atom file itself (data/atomic/excitation/mocassin_v2/
nlevel_fe_2.txt), evaluates the strict n_e -> 0 statistical-equilibrium
cooling (chianti_cooling.cooling_low_density_from_levels, the population
model of every 'mocassin_v2' fit) and fits it with the signed-exponential
linear program of mocassin_cooling_fits.exponential_fit (acceptance:
relative error < 0.5% on 1601 temperatures over 1e3-1e5 K, cancellation
< 1e5), the method that produced the other 'mocassin_v2' cooling files,
with nonnegativity on 1-1e9 K added as fit_cooling.py requires (the plain
exponential_fit solution for this atom turns negative above 5.8e5 K).
fit_cooling.py builds its sources from the CHIANTI or Stout databases and
cannot read a finished atom file, so it is not used here; its checks are
repeated: the error on 801 temperatures over 1e3-1e5 K and positivity on
1-1e9 K.

The fit file and the manifest SHA-256 entries of the atom and the fit are
rewritten together.

Run:  python3 tools/fitting/feii_low_density_cooling.py [--report REPORT.json]
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import chianti_cooling as C
from atomic_output import atomic_write_text
from mocassin_cooling_fits import log_cooling

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / 'data/atomic/excitation/mocassin_v2'
ATOM = MODEL / 'nlevel_fe_2.txt'
FIT = MODEL / 'cooling_fe_2.txt'
MANIFEST = MODEL / 'manifest.json'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_temperature_table_atom(path):
    """Levels, radiation and TEMPERATURE_TABLE collision rows of a MoCHII atom."""
    records = [s.split() for s in Path(path).read_text().splitlines()
               if s.strip() and not s.lstrip().startswith('#')]
    k = next(i for i, r in enumerate(records) if r[0] == 'NLEV')
    n = int(records[k][1])
    levels = {int(r[0]): (float(r[1]), float(r[2])) for r in records[k+1:k+1+n]}
    k = next(i for i, r in enumerate(records) if r[0] == 'NRAD')
    n = int(records[k][1])
    radiation = []
    for r in records[k+1:k+1+n]:
        l, u, a = int(r[0]), int(r[1]), float(r[2])
        if levels[u][1] < levels[l][1]:
            l, u = u, l
        radiation.append((l, u, 0., 0., a))
    k = next(i for i, r in enumerate(records) if r[0] == 'NUPS')
    n = int(records[k][1])
    collisions = []
    i = k + 1
    for _ in range(n):
        r = records[i]
        if r[2] != 'TEMPERATURE_TABLE' or r[6] != 'LOG_PCHIP_T':
            raise ValueError(f'Unsupported collision row {r}')
        count = int(r[7])
        nodes = np.array([[float(x) for x in s[:2]] for s in records[i+1:i+1+count]])
        collisions.append(dict(ll=int(r[0]), ul=int(r[1]), temperature_grid=nodes[:, 0],
                               upsilon_grid=nodes[:, 1]))
        i += 1 + count
    return levels, radiation, collisions


def read_fit(path):
    rows = [s.split() for s in Path(path).read_text().splitlines()
            if s.strip() and not s.lstrip().startswith('#')]
    n = int(rows[0][0])
    data = np.array([[float(x) for x in r[:2]] for r in rows[1:1+n]])
    return data[:, 0], data[:, 1]


def signed_exponential_fit(levels, radiation, collisions, maximum=1e5):
    """mocassin_cooling_fits.exponential_fit with positivity on 1-1e9 K.

    The same linear program (relative error bound on 401 temperatures,
    cancellation < 1e5, acceptance < 0.5% on 1601 temperatures), with the
    cancellation rows of fit_cooling.constrained_signed_exponentials added on
    401 temperatures over 1-1e9 K so the expansion stays nonnegative outside
    the fit interval; positivity is then checked on 10001 temperatures.
    """
    temperature = np.geomspace(1000., maximum, 401)
    dense = np.geomspace(1000., maximum, 1601)
    extended = np.logspace(0, 9, 401)
    check = np.logspace(0, 9, 10001)
    log_source = log_cooling(levels, radiation, collisions, temperature)
    dense_source = log_cooling(levels, radiation, collisions, dense)
    closed, probability = C.radiative_cascade_probabilities(levels, collisions, radiation)
    end_energy = np.array([levels[i][1] for i in closed])
    thresholds = []
    for tr in collisions:
        lo, up = C.excitation_levels(levels, tr)
        if lo in closed and levels[up][1]-np.dot(probability[up], end_energy) > 0:
            energy = levels[up][1]-min(e for g, e in levels.values())
            if energy > 0:
                thresholds.append(energy*C.HC_ERG_CM/C.K_B_ERG)
    minimum = min(thresholds)
    target = log_source+.5*np.log(temperature)+minimum/temperature
    scale = float(target.max())
    normalized = np.exp(target-scale)
    for count in (12, 18, 24, 32, 48, 64, 96, 128):
        for spacing in (.00125, .0025, .005, .01, .02, .04, .08):
            offsets = np.arange(count)*spacing*maximum
            basis = np.exp(-offsets[None, :]/temperature[:, None])
            outer = np.exp(-offsets[None, :]/extended[:, None])
            relative = basis/normalized[:, None]
            ones, zeros = np.ones(len(temperature)), np.zeros(len(temperature))
            constraints = np.vstack((
                np.column_stack((relative, -relative, -ones)),
                np.column_stack((-relative, relative, -ones)),
                np.column_stack((basis*(1e-5-1), basis*(1e-5+1), zeros)),
                np.column_stack((outer*(1e-5-1), outer*(1e-5+1), np.zeros(len(extended))))))
            rhs = np.r_[ones, -ones, zeros, np.zeros(len(extended))]
            solution = linprog(np.r_[np.zeros(2*count), 1.], A_ub=constraints, b_ub=rhs,
                               bounds=[(0, None)]*(2*count+1), method='highs')
            if not solution.success:
                continue
            coefficient = (solution.x[:count]-solution.x[count:2*count])*np.exp(scale)
            excitation = minimum+offsets
            parts = coefficient[:, None]*np.exp(-(excitation[:, None]-minimum)/dense)
            total = parts.sum(axis=0)
            if np.any(total <= 0):
                continue
            condition = float(np.max(np.abs(parts).sum(axis=0)/total))
            log_fit = np.log(total)-minimum/dense-.5*np.log(dense)
            error = float(np.max(np.abs(np.expm1(log_fit-dense_source))))
            keep = coefficient != 0
            wide = evaluate(coefficient[keep], excitation[keep], check)
            if (condition < 1e5 and error < .005 and np.all(np.isfinite(wide))
                    and np.all(wide >= 0)):
                return coefficient[keep], excitation[keep], error, condition
    raise ValueError('No accepted nonnegative exponential fit')


def evaluate(amplitude, excitation, temperature):
    temperature = np.atleast_1d(temperature)
    base = excitation.min()
    terms = amplitude[:, None]*np.exp(-(excitation[:, None]-base)/temperature)
    return terms.sum(axis=0)*np.exp(-base/temperature)/np.sqrt(temperature)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--report', type=Path)
    parser.add_argument('--fortran-cooling', type=Path,
                        help='excitation_atom_cooling.f90 output for this atom (n_e = 1e-6 column)')
    args = parser.parse_args()
    levels, radiation, collisions = read_temperature_table_atom(ATOM)
    amplitude, excitation, error_1601, condition = signed_exponential_fit(
        levels, radiation, collisions)

    t801 = np.logspace(3, 5, 801)
    source_801 = log_cooling(levels, radiation, collisions, t801)
    fitted_801 = evaluate(amplitude, excitation, t801)
    error_801 = float(np.max(np.abs(np.expm1(np.log(fitted_801)-source_801))))
    wide = evaluate(amplitude, excitation, np.logspace(0, 9, 10001))
    if not (np.all(np.isfinite(wide)) and np.all(wide >= 0)):
        raise ValueError('Fit is negative or nonfinite on 1-1e9 K')
    if error_801 >= 0.005:
        raise ValueError(f'801-temperature error {error_801} exceeds 0.5%')

    previous = read_fit(FIT)
    lines = [
        '# Fe II, MOCASSIN v2 142-level atom: 1144 supplied radiative rows, 40 earlier forbidden A values',
        '# and 141 CHIANTI 11.0.2 forbidden A values from levels 14-35 (feii_spontaneous_probabilities.py).',
        '# This fit uses the same levels, radiation and collision tables as the atom file.',
        '# Model interval: 1000 <= T <= 100000 K; physical accuracy remains unverified.',
        '# Strict low-density statistical equilibrium, including radiatively closed states.',
        '# Signed terms represent total positive cooling, not individual emission rates.',
        f'# Maximum relative error: {error_1601:.6g} on 1601 temperatures, {error_801:.6g} on 801; '
        f'cancellation {condition:.3g}; nonnegative on 1-1e9 K.',
        '# generated by tools/fitting/feii_low_density_cooling.py (exponential_fit with nonnegativity on 1-1e9 K)',
        '# Source atom SHA-256: ' + sha256(ATOM),
        str(len(amplitude))]
    lines += [f'{a:.17g} {e:.17g}' for a, e in zip(amplitude, excitation)]
    atomic_write_text(FIT, '\n'.join(lines) + '\n')

    manifest = json.loads(MANIFEST.read_text())
    for record in manifest['files']:
        if record['file'] in ('nlevel_fe_2.txt', 'cooling_fe_2.txt'):
            record['sha256'] = sha256(MODEL / record['file'])
    atomic_write_text(MANIFEST, json.dumps(manifest, indent=2) + '\n')
    for record in json.loads(MANIFEST.read_text())['files']:
        if sha256(MODEL / record['file']) != record['sha256']:
            raise ValueError(f'Manifest disagrees with {record["file"]}')

    chianti = read_fit(ROOT / 'data/atomic/excitation/chianti/cooling_fe_2.txt')
    probe = np.array([5e3, 1e4, 2e4])
    new = evaluate(amplitude, excitation, probe)
    report = dict(terms=len(amplitude), error_1601=error_1601, error_801=error_801,
                  cancellation=condition, atom_sha256=sha256(ATOM), fit_sha256=sha256(FIT),
                  temperature_K=probe.tolist(), new_fit=new.tolist(),
                  new_over_previous_fit=(new/evaluate(*previous, probe)).tolist(),
                  new_over_chianti_fit=(new/evaluate(*chianti, probe)).tolist())
    if args.fortran_cooling:
        data = np.loadtxt(args.fortran_cooling)
        source = np.exp(log_cooling(levels, radiation, collisions, data[:, 0]))
        report['fortran_low_density_over_python_max_deviation'] = float(
            np.max(np.abs(data[:, 1]/source-1)))
        report['fortran_low_density_over_fit_max_deviation'] = float(
            np.max(np.abs(data[:, 1]/evaluate(amplitude, excitation, data[:, 0])-1)))
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
