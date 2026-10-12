"""Check every level atom and cooling file of both excitation models.

For each model in data/atomic/excitation/ (chianti, mocassin_v2):
  populations  every atom: statistical equilibrium of the level populations
               (each level's gains equal its losses to 1e-10 of that level's
               own throughput, plus the underflow floor (n+1) tiny times the
               level's rates, tiny = 2.2e-308; see
               check_source_excitation_atom.f90), unit normalization, and
               the emitted power
               sum n_u A_ul dE_ul against nlevel_cooling_solve and
               nlevel_emissivities, at 41 temperatures 1e2-1e6 K and n_e =
               1e-6, 1e2, 1e6, 1e10 cm^-3 (check_source_excitation_atom.f90);
  registered   every atom: the (T, n_e) cooling table against the direct
               n-level solve, 1e-10 at its nodes and 0.5% between them, and
               the fit-combination rule (check_registered_excitation_atom.f90,
               check_small_excitation_cooling.f90);
  fit          every cooling file: cooling_eval against the formula
               Lambda = max(sum A_i exp(-T_i/T), 0)/sqrt(T) evaluated here
               (check_cooling_fit.f90), and against the low-density
               (n_e = 1e-6 cm^-3) cooling of the same stage's atom on 801
               temperatures, 1e3-1e5 K, the interval the files state they
               were fitted on, within the maximum fit error each file states
               in its header (CHIANTI files print it in percent to two
               decimals, so 5e-5 is added for that rounding), where Lambda
               exceeds 1e-3 of its maximum there, and within that error times
               1e-3 of the maximum elsewhere; files whose stated error exceeds
               0.5% (the table interpolation criterion; the acceptance target
               of tools/fitting/fit_cooling.py is 1%) are listed
               (check_bound_excitation_atom.f90);
  no fit       every atom without a cooling file: the species cooling of a
               one-element gas against the direct cooling of its atoms,
               0.5%, in both cooling models (check_atomic_cooling_without_fit.f90).
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import re
from pathlib import Path
import shutil

import numpy as np

from run_checks import ROOT, link, run

MODELS = ('chianti', 'mocassin_v2')
WORKERS = 4
FIT_TOLERANCE = 5e-3


def atoms(model):
    directory = ROOT/'data/atomic/excitation'/model
    names = sorted(p.stem.removeprefix('nlevel_') for p in directory.glob('nlevel_*.txt'))
    return [(name.rsplit('_', 1)[0], int(name.rsplit('_', 1)[1])) for name in names]


def stated_fit_error(path):
    """The maximum fit error a cooling file states in its header."""
    text = path.read_text()
    found = re.search(r'max fit error on these temperatures: ([0-9.]+)%', text)
    if found:
        return float(found.group(1))/100+5e-5
    found = re.search(r'Maximum relative error(?: on \d+ temperatures)?: (\d+\.\d+(?:[eE][+-]?\d+)?)', text)
    if found:
        return float(found.group(1))*(1+1e-6)
    raise RuntimeError(f'{path.name} states no fit error')


def fit_formula(path, temperature):
    rows = [s.split() for s in path.read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
    n = int(rows[0][0])
    coefficients = np.array(rows[1:1+n], dtype=float)
    terms = coefficients[:, 0, None]*np.exp(-coefficients[:, 1, None]/temperature)
    return np.maximum(terms.sum(axis=0), 0)/np.sqrt(temperature), np.abs(terms).sum(axis=0)/np.sqrt(temperature)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = {name: link(work/'build', f'tests/ion_stages/{name}.f90') for name in (
        'check_source_excitation_atom', 'check_registered_excitation_atom',
        'check_small_excitation_cooling', 'check_cooling_fit', 'check_bound_excitation_atom',
        'check_atomic_cooling_without_fit')}
    grid = work/'temperatures.in'
    grid.write_text('0 41\n'+' '.join(f'{10**(2+0.1*i):.17g}' for i in range(41))+'\n')
    failures = []
    summary = {}

    def atom_checks(model, element, stage):
        tag = f'{model}_{element}_{stage}'
        directory = work/tag
        directory.mkdir()
        model_dir = ROOT/'data/atomic/excitation'/model
        atom = model_dir/f'nlevel_{element}_{stage}.txt'
        record = {}
        try:
            # A balance residual above 1e-10 of the level throughput is recorded
            # as a failure; the other checks of the atom still run.
            try:
                log = run([exe['check_source_excitation_atom'], atom, grid, directory/'populations.txt'],
                          log=directory/'populations.log')
            except RuntimeError as error:
                log = str(error)
                if 'Populations not stationary' not in log:
                    raise
            balance = log.split('BALANCE worst relative residual')[1].split()
            record['balance_residual'] = float(balance[0])
            record['balance_worst_at'] = dict(T=float(balance[3]), n_e=float(balance[5]), level=int(balance[7]),
                                             level_throughput=float(balance[11]))
            if record['balance_residual'] > 1e-10:
                failures.append(f'{tag}: level balance residual {balance[0]} of the level throughput '
                                f'at T {balance[3]} K, n_e {balance[5]}, level {balance[7]}, throughput {balance[11]} s^-1')
            (directory/'populations.txt').unlink(missing_ok=True)
            run([exe['check_registered_excitation_atom'], ROOT/'data/atomic', model, element, stage,
                 directory/'registered.txt'], log=directory/'registered.log')
            text = (directory/'registered.txt').read_text()
            for key in ('worst_grid_relative_error', 'worst_middle_relative_error',
                        'worst_apply_relative_error'):
                record[key] = float(text.split('# '+key)[1].split()[0])
            run([exe['check_small_excitation_cooling'], ROOT/'data/atomic', model, element, stage,
                 directory/'small.txt'], log=directory/'small.log')
            fit = model_dir/f'cooling_{element}_{stage}.txt'
            if fit.exists():
                run([exe['check_cooling_fit'], fit, directory/'fit.txt'], log=directory/'fit.log')
                data = np.loadtxt(directory/'fit.txt')
                expected, scale = fit_formula(fit, data[:, 0])
                # Relative to the sum of the term magnitudes, with tiny/eps
                # (1.0e-292 erg cm^3 s^-1) as the floor: the production build
                # flushes intermediates below tiny to zero, so within a factor
                # 1/eps of tiny cooling_eval returns 0 where the formula is
                # still positive (measured: C V, C VI, N VI, N VII, O VII,
                # S VII, Si V between 1800 and 9700 K, formula 5e-324 ..
                # 1e-308, Fortran 0); above the floor the two agree to
                # 1.1e-13. Without a floor, 0/0 where every term underflows
                # gives NaN, which no comparison rejects.
                difference = np.abs(data[:, 1]-expected)
                floor = np.finfo(float).tiny/np.finfo(float).eps
                record['fit_formula_error'] = float(np.max(difference/np.maximum(scale, floor)))
                if record['fit_formula_error'] > 1e-12:
                    raise RuntimeError(f'cooling_eval differs from the fit formula: {record}')
                run([exe['check_bound_excitation_atom'], fit, atom, directory/'bound.txt', 'cooling'],
                    log=directory/'bound.log')
                data = np.loadtxt(directory/'bound.txt')
                low = data[data[:, 1] == 1e-6]
                direct, fitted = low[:, 3], low[:, 2]
                significant = direct > 1e-3*direct.max()
                record['fit_atom_error'] = float(np.max(np.abs(fitted[significant]/direct[significant]-1)))
                record['fit_atom_error_below_significance'] = float(
                    np.max(np.abs(fitted[~significant]-direct[~significant]))/direct.max()) \
                    if (~significant).any() else 0.
                record['stated_fit_error'] = stated = stated_fit_error(fit)
                if record['fit_atom_error'] > stated or \
                        record['fit_atom_error_below_significance'] > stated*1e-3:
                    raise RuntimeError(f'cooling file differs from its atom beyond its stated error: {record}')
            if record['balance_residual'] <= 1e-10:
                shutil.rmtree(directory)
        except Exception as error:
            failures.append(f'{tag}: {str(error)[-1500:]}')
        return tag, record

    jobs = [(model, element, stage) for model in MODELS for element, stage in atoms(model)]
    with ThreadPoolExecutor(WORKERS) as pool:
        for tag, record in pool.map(lambda job: atom_checks(*job), jobs):
            summary[tag] = record

    # Stages with a level atom and no cooling file, in a fixture holding only
    # level atoms: the element file, the stage's atom, and, for the form-9
    # sulfur charge exchange, the S II atom it reads.
    no_fit = {}
    for model in MODELS:
        model_dir = ROOT/'data/atomic/excitation'/model
        for element, stage in atoms(model):
            if (model_dir/f'cooling_{element}_{stage}.txt').exists():
                continue
            tag = f'{model}_{element}_{stage}'
            fixture = work/'no_fit'/tag
            target = fixture/'excitation'/model
            target.mkdir(parents=True)
            (target/'manifest.json').write_text('{"files": [], "scope": "level atoms without cooling files"}\n')
            shutil.copy2(ROOT/f'data/atomic/element_{element}.txt', fixture)
            for name in {f'nlevel_{element}_{stage}.txt', *(['nlevel_s_2.txt'] if element == 's' else [])}:
                shutil.copy2(model_dir/name, target)
            for mode in ('local_ne', 'low_density'):
                inp = fixture/f'{mode}.in'
                inp.write_text(f"&parameters\n par%metal_excitation_model='{model}'\n"
                               f" par%cooling_model='{mode}'\n/\n")
                try:
                    log = run([exe['check_atomic_cooling_without_fit'], fixture, element, inp,
                               fixture/f'{mode}.txt'], log=fixture/f'{mode}.log')
                    if 'uses absolute n-level cooling without an exponential fit' not in log:
                        raise RuntimeError('the stage was not registered for absolute cooling')
                    data = np.loadtxt(fixture/f'{mode}.txt')
                    error = float(np.max(np.abs(data[:, 2]/data[:, 3]-1)))
                    no_fit[f'{tag}_{mode}'] = error
                    if data.shape != (132, 4) or error > FIT_TOLERANCE:
                        raise RuntimeError(f'species cooling differs from the atoms by {error}')
                except Exception as error:
                    failures.append(f'{tag} {mode}: {str(error)[-1500:]}')
    worst = {key: max((r.get(key, 0.) for r in summary.values()), default=0.)
             for key in ('balance_residual', 'worst_grid_relative_error',
                         'worst_middle_relative_error', 'fit_formula_error', 'fit_atom_error',
                         'fit_atom_error_below_significance')}
    above_target = sorted(k for k, r in summary.items() if r.get('stated_fit_error', 0.) > FIT_TOLERANCE)
    result = dict(passed=not failures, stated_fit_error_above_target=above_target, atoms={m: len(atoms(m)) for m in MODELS},
                  cooling_files={m: sum((ROOT/'data/atomic/excitation'/m/f'cooling_{e}_{s}.txt').exists()
                                        for e, s in atoms(m)) for m in MODELS},
                  worst=worst, no_fit_relative_error=no_fit, failures=failures, atoms_detail=summary)
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'atoms_detail'}, indent=1))


if __name__ == '__main__':
    main()
