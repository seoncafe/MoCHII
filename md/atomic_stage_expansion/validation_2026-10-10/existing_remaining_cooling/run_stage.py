"""Stage one source cooling fit and check it with an immutable native reader.

These are numerical diagnostics. Published temperature coverage and atomic
accuracy require a separate physical review before installing any product.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


STAGES = {
    'NeI': ('ne', 1),
    'MgI': ('mg', 1),
    'MgIII': ('mg', 3),
    'SiI': ('si', 1),
    'SiV': ('si', 5),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ion', choices=STAGES, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--driver', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[4]
    output = args.output_dir.resolve()
    if any(output.iterdir()) if output.exists() else False:
        raise RuntimeError('Preserve existing candidate products; use a new directory')
    output.mkdir(parents=True, exist_ok=True)
    element, stage = STAGES[args.ion]
    fit_script = root / 'tools/fitting/validate_stage_cooling.py'
    native_script = root / 'tools/fitting/validate_native_cooling_coefficients.py'
    subprocess.run([sys.executable, str(fit_script), '--ion', args.ion,
                    '--element', element, '--stage', str(stage),
                    '--output-dir', str(output)], cwd=root, check=True,
                   timeout=1200)
    coefficients = output / f'cooling_{element}_{stage}.txt'
    fit_result_path = output / f'{args.ion}.json'
    figure_path = output / f'{args.ion}.pdf'
    fit_result = json.loads(fit_result_path.read_text())
    if (not fit_result.get('passed') or fit_result['ion'] != args.ion
        or not coefficients.is_file() or figure_path.stat().st_size < 1000):
        raise ValueError('Incomplete source-fit diagnostic')
    native = output / 'native'
    subprocess.run([sys.executable, str(native_script), '--coefficients',
                    str(coefficients), '--driver', str(args.driver.resolve()),
                    '--work-dir', str(native)], cwd=root, check=True,
                   timeout=180)
    native_result = json.loads((native / 'result.json').read_text())
    if (not native_result.get('passed') or native_result.get('points') != 801
        or not (native / 'native_cooling.txt').is_file()):
        raise ValueError('Incomplete native cooling-reader diagnostic')
    result = {
        'passed': True,
        'meaning': 'Numerical fit and native reader diagnostics; physical source acceptance pending',
        'ion': args.ion,
        'source_fit': fit_result,
        'native_reader': native_result,
        'outputs_sha256': {str(path): digest(path) for path in
                           (coefficients, fit_result_path, figure_path,
                            native / 'result.json', native / 'native_cooling.txt')},
    }
    (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(output / 'result.json')


if __name__ == '__main__':
    main()
