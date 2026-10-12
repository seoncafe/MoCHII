"""Fit C IV only over its published 3,200--100,000 K overlap interval."""

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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output = Path(sys.argv[1]).resolve()
    driver = Path(sys.argv[2]).resolve()
    if output.exists() and any(output.iterdir()):
        raise RuntimeError('Preserve the existing C IV candidate products')
    output.mkdir(parents=True, exist_ok=True)
    F.TMIN = 3200.0
    F.TMAX = 100000.0
    F.T = np.geomspace(F.TMIN, F.TMAX, 201)
    F.generate({'CIV': ('c', 4)}, output)
    product = output / 'cooling_c_4.txt'
    amplitude, excitation = coefficients(product)
    grid = np.geomspace(F.TMIN, F.TMAX, 801)
    source = C.cooling_low_density_limit('c', 4, grid)
    fitted = F.evaluate(amplitude, excitation, grid)
    if not np.all(np.isfinite(source)) or np.any(source <= 0):
        raise ValueError('Invalid source cooling in published C IV interval')
    difference = fitted / source - 1
    error = float(np.max(np.abs(difference)))
    if not np.all(np.isfinite(fitted)) or np.any(fitted < 0) or error > F.ERR_TARGET:
        raise ValueError(f'C IV source-range fit error: {error}')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True,
                                gridspec_kw={'height_ratios': [3, 1]})
    axes[0].loglog(grid, source, label='Source')
    axes[0].loglog(grid, fitted, '--', label='Fit')
    axes[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
    axes[0].set_title('C IV: published temperature overlap')
    axes[0].legend()
    axes[1].semilogx(grid, 100*difference)
    axes[1].set_xlabel(r'$T$ [K]')
    axes[1].set_ylabel(r'Fit error [\%]')
    figure.tight_layout()
    figure.savefig(output / 'CIV.pdf')
    plt.close(figure)
    native = output / 'native'
    subprocess.run([sys.executable,
                    str(ROOT / 'tools/fitting/validate_native_cooling_coefficients.py'),
                    '--coefficients', str(product), '--driver', str(driver),
                    '--work-dir', str(native)], cwd=ROOT, check=True, timeout=180)
    native_result = json.loads((native / 'result.json').read_text())
    if not native_result['passed'] or (output / 'CIV.pdf').stat().st_size < 1000:
        raise ValueError('Incomplete C IV numerical diagnostics')
    result = {
        'passed': True, 'ion': 'CIV',
        'temperature_interval_K': [3200.0, 100000.0],
        'source_interval_K': [3200.0, 32000000.0],
        'max_relative_fit_error': error,
        'native_reader': native_result,
        'physical_acceptance': False,
        'reason': 'The 1,000--3,200 K part of the intended interval lacks source Upsilon; no cold continuation is implemented.',
        'outputs_sha256': {str(path): digest(path) for path in
                           (product, output / 'CIV.pdf', native / 'result.json',
                            native / 'native_cooling.txt')},
    }
    (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(output / 'result.json')


if __name__ == '__main__':
    main()
