"""Fit Ne-like Mg III or Si V only within published temperature coverage."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'tools/fitting'))
import chianti_cooling as C
import fit_cooling as F
from validate_stage_cooling import coefficients

IONS = {'MgIII': ('mg', 3, 1800.0), 'SiV': ('si', 5, 5000.0)}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ion', choices=IONS, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--driver', type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise RuntimeError('Preserve the existing candidate products')
    output.mkdir(parents=True, exist_ok=True)
    element, stage, minimum = IONS[args.ion]
    F.TMIN = minimum
    F.TMAX = 100000.0
    F.T = np.geomspace(minimum, F.TMAX, 201)
    levels, radiation, collisions = C.bound_excitation_atom(element, stage)
    closed, _ = C.radiative_cascade_probabilities(levels, collisions, radiation)
    if closed != [1]:
        raise ValueError('This diagnostic requires a single radiatively closed ground state')
    if any(C.excitation_levels(levels, row)[0] == 1
           and C.transition_energy_erg(levels, row) <= 0 for row in collisions):
        raise ValueError('Nonpositive ground excitation gap')
    F.generate({args.ion: (element, stage)}, output)
    product = output / f'cooling_{element}_{stage}.txt'
    amplitude, excitation = coefficients(product)
    grid = np.geomspace(minimum, F.TMAX, 801)
    source = C.cooling_low_density_limit(element, stage, grid)
    fitted = F.evaluate(amplitude, excitation, grid)
    if not np.all(np.isfinite(source)) or np.any(source <= 0):
        raise ValueError('Invalid source cooling on published temperature overlap')
    difference = fitted / source - 1
    error = float(np.max(np.abs(difference)))
    if not np.all(np.isfinite(fitted)) or np.any(fitted < 0) or error > F.ERR_TARGET:
        raise ValueError(f'{args.ion} source-range fit error: {error}')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True,
                                gridspec_kw={'height_ratios': [3, 1]})
    axes[0].loglog(grid, source, label='Source')
    axes[0].loglog(grid, fitted, '--', label='Fit')
    axes[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
    axes[0].set_title(args.ion + ': published temperature overlap')
    axes[0].legend()
    axes[1].semilogx(grid, 100*difference)
    axes[1].set_xlabel(r'$T$ [K]')
    axes[1].set_ylabel(r'Fit error [\%]')
    figure.tight_layout()
    figure.savefig(output / f'{args.ion}.pdf')
    plt.close(figure)
    native = output / 'native'
    subprocess.run([sys.executable,
                    str(ROOT / 'tools/fitting/validate_native_cooling_coefficients.py'),
                    '--coefficients', str(product), '--driver', str(args.driver.resolve()),
                    '--work-dir', str(native)], cwd=ROOT, check=True, timeout=180)
    native_result = json.loads((native / 'result.json').read_text())
    figure_path = output / f'{args.ion}.pdf'
    if not native_result['passed'] or figure_path.stat().st_size < 1000:
        raise ValueError('Incomplete native or figure diagnostic')
    result = {
        'passed': True, 'ion': args.ion,
        'temperature_interval_K': [minimum, F.TMAX],
        'max_relative_fit_error': error,
        'native_reader': native_result,
        'physical_acceptance': False,
        'reason': 'Published lower-temperature boundary is above 1,000 K; source radiation/energy identity requires separate native excitation-atom review.',
        'outputs_sha256': {str(path): digest(path) for path in
                           (product, figure_path, native / 'result.json',
                            native / 'native_cooling.txt')},
    }
    (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(output / 'result.json')


if __name__ == '__main__':
    main()
