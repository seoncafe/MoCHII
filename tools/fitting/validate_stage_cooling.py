"""Generate and validate one stationary low-density cooling fit in staging."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import chianti_cooling as C
import fit_cooling as F
from atomic_output import atomic_write_text
from write_excitation_atom import _source_files


def coefficients(path):
    rows = [line.split() for line in path.read_text().splitlines()
            if line.strip() and not line.startswith('#')]
    number = int(rows[0][0])
    values = np.array(rows[1:], dtype=float)
    if values.shape != (number, 2):
        raise ValueError(f'Incomplete coefficient file: {path}')
    return values[:, 0], values[:, 1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ion', required=True)
    parser.add_argument('--element', required=True)
    parser.add_argument('--stage', type=int, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    F.generate({args.ion: (args.element, args.stage)}, directory)
    product = directory/f'cooling_{args.element}_{args.stage}.txt'
    amplitude, temperature = coefficients(product)
    grid = np.logspace(3, 5, 801)
    source = C.cooling_low_density_limit(args.element, args.stage, grid)
    fitted = F.evaluate(amplitude, temperature, grid)
    underflow=bool(np.any(source==0))
    if underflow:
        log_source=C.log_cooling_ground_state_limit(args.element,args.stage,grid)
        log_fit=F.log_cooling_fit(amplitude,temperature,grid)
        if not np.all(np.isfinite(log_fit)):raise ValueError('Nonpositive logarithmic cooling fit')
        differences=np.expm1(log_fit-log_source)
    else:
        differences=fitted/source-1
    error = np.max(np.abs(differences))
    if not np.all(np.isfinite(fitted)) or np.any(fitted < 0) or error > F.ERR_TARGET:
        raise ValueError(f'Rounded fit fails dense-grid acceptance: {error:.8g}')
    levels, radiation, collisions = C.bound_excitation_atom(args.element,args.stage)
    source_paths = _source_files(args.element.lower(), args.stage)
    closed, _ = C.radiative_cascade_probabilities(
        levels, collisions, radiation)
    previous = Path(F.OUTDIR)/product.name
    previous_values = None
    if previous.exists():
        a, t = coefficients(previous)
        previous_values = F.evaluate(a, t, grid)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True,
                                gridspec_kw={'height_ratios': [3, 1]})
    axes[0].loglog(grid, source, label='Stationary source')
    axes[0].loglog(grid, fitted, '--', label='Rounded fit')
    if previous_values is not None:
        axes[0].loglog(grid, previous_values, ':', label='Previous production fit')
    axes[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
    axes[0].set_title(args.ion)
    axes[0].legend(fontsize=9)
    axes[1].semilogx(grid, 100*differences)
    axes[1].set_ylabel(r'Fit error [\%]')
    axes[1].set_xlabel(r'$T$ [K]')
    figure.tight_layout()
    figure.savefig(directory/f'{args.ion}.pdf')
    plt.close(figure)
    result = dict(passed=True, ion=args.ion, chianti_version=F.CHIANTI_VERSION,
                  temperatures_K=[1000, 100000], dense_grid_points=801,
                  max_relative_fit_error=float(error), radiatively_closed_states=closed,
                  source_underflow_points=int(np.sum(source==0)),
                  underflow_treatment='Logarithmic one-end-state excitation sum; no rate floor' if underflow else None,
                  population='stationary electron collisions and spontaneous radiation, n_e -> 0',
                  excluded_physics=['photoexcitation', 'excited-state ionization', 'finite-density suppression'],
                  sources={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths},
                  product_sha256=hashlib.sha256(product.read_bytes()).hexdigest())
    if (args.element.lower(),args.stage) == ('cl',7):
        result.update(chianti_version=None,excitation_source=f'Cloudy c25.00 Stout, {len(levels)} bound states',
            cold_approximation='Constant Upsilon below 9800 K; actual T in rates; unverified resonance accuracy')
    if (args.element.lower(),args.stage) in [('ar',2),('ca',4)]:
        first = float(collisions[0]['temperature_grid'][0])
        last = float(collisions[0]['temperature_grid'][-1])
        result.update(chianti_version=None,
            excitation_source='Pelan and Berrington (1995), A&AS 110, 209, Table 3; Cloudy c25.00 Stout ground pair',
            source_temperature_interval_K=[first, last],
            published_log10_temperature_interval=(
                [3.2, 5.2] if args.element.lower() == 'ar' else [4.0, 6.0]),
            cold_approximation=f'Constant first-node Upsilon below {first:g} K; actual T in Maxwellian rates',
            above_source_temperature='Rejected above the last source node; no high-temperature continuation',
            source_accuracy=('The published better-than-10-percent estimate explicitly excludes Ar II' if args.element.lower() == 'ar'
                             else 'The published better-than-10-percent estimate includes Ca IV within Table 3 range'),
            selected_levels=2, electron_collision_rows=1,
            omitted_higher_levels='No identified electron collision strengths for the higher radiative levels')
    if args.element.lower() == 'n' and args.stage == 4:
        from niv_level_energies import STOUT_NRG, STOUT_SHA256
        result['sources'][str(STOUT_NRG)] = STOUT_SHA256
        result['energy_accuracy'] = ('Observed energies at levels 1--29; BSR above. '
                                     'Low-temperature resonance uncertainty and '
                                     'convergence above level 87 remain unquantified. '
                                     'The fit interval does not guarantee physical accuracy.')
    if (args.element.lower(),args.stage) in [('o',5),('cl',6)]:
        result['energy_accuracy']='Observed low levels and ICFT above; scattering resonance positions and higher-level convergence are not certified by the fit interval.'
        if args.element.lower()=='o':
            from atomic_excitation_energies import STOUT_OV,STOUT_OV_SHA256
            result['sources'][str(STOUT_OV)]=STOUT_OV_SHA256
    if previous_values is not None:
        difference = source/previous_values-1
        result.update(previous_sha256=hashlib.sha256(previous.read_bytes()).hexdigest(),
                      source_change_from_previous_min=float(difference.min()),
                      source_change_from_previous_max=float(difference.max()),
                      source_change_at_10000K=float(difference[400]))
    atomic_write_text(directory/f'{args.ion}.json', json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
