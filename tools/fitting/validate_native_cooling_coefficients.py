"""Check the actual compiled cooling reader against serialized coefficients."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from atomic_output import atomic_write_text
from validate_stage_cooling import coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coefficients', type=Path, required=True)
    parser.add_argument('--driver', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    directory = args.work_dir.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    output = directory / 'native_cooling.txt'
    if output.exists() or (directory / 'result.json').exists():
        raise RuntimeError('Preserve prior products; use a new validation directory')
    subprocess.run([str(args.driver.resolve()), str(args.coefficients.resolve()),
                    str(output)], check=True, timeout=120)
    rows = np.loadtxt(output)
    if rows.shape != (801, 2) or not np.all(np.isfinite(rows)) or np.any(rows[:, 1] < 0):
        raise ValueError('Incomplete or invalid native cooling coefficients')
    if not np.allclose(rows[:, 0], np.geomspace(1e3, 1e5, 801), rtol=1e-13, atol=0):
        raise ValueError('Unexpected native temperature grid')
    amplitude, excitation = coefficients(args.coefficients)
    expected = np.sum(amplitude[:, None] * np.exp(-excitation[:, None]/rows[:, 0]),
                      axis=0) / np.sqrt(rows[:, 0])
    if not np.all(np.isfinite(expected)) or np.any(expected < 0):
        raise ValueError('Nonphysical serialized cooling fit')
    # Optimized builds may flush subnormal intermediate exponential terms.
    # The guard allows five orders for accepted signed-sum conditioning below
    # 1e5; it is a numerical comparison range, not a physical cooling floor.
    comparison_limit = 1e5 * np.finfo(float).tiny
    resolved = expected >= comparison_limit
    error = float(np.max(np.abs(rows[resolved, 1]/expected[resolved]-1))) if resolved.any() else 0.
    if error > 1e-9 or np.any(rows[expected == 0, 1] != 0):
        raise ValueError(f'Native reader differs from serialized coefficients: {error}')
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    result = dict(passed=True, points=len(rows), resolved_points=int(resolved.sum()),
                  numerical_comparison_limit=comparison_limit,
                  maximum_resolved_relative_error=error,
                  native_zero_vs_positive_unresolved=int(np.sum((rows[:, 1] == 0) & (expected > 0) & ~resolved)),
                  coefficients_path=str(args.coefficients.resolve()),
                  coefficients_sha256=digest(args.coefficients),
                  driver_path=str(args.driver.resolve()), driver_sha256=digest(args.driver),
                  scope='Serialized coefficients and actual cooling_load/cooling_eval; no source-fit accuracy or spatial-model claim')
    atomic_write_text(directory / 'result.json', json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
