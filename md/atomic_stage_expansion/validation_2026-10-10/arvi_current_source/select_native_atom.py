"""Select the first Ar VI bound atom meeting the measured source cutoff."""

import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT/'tools/fitting'))
import chianti_cooling as C
from validate_excitation_density_profiles import full_excitation_balance
from validate_source_excitation_atom import DENSITY, source_subset


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--driver', type=Path, required=True)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
    base = args.work_dir.resolve()
    full_result = json.loads((base/'fullsource/result.json').read_text())
    profiles = base/'fullsource/full_source_profiles.npz'
    full = np.load(profiles)
    temperature = full['temperature_K']
    levels, radiation, collisions = C.bound_excitation_atom('ar', 6)
    candidates = [dict(levels=40,
                       maximum_relative_cutoff_errors=full_result['maximum_relative_cutoff_errors'])]
    selected_count = 40 if full_result['selected_40_meets_cutoff'] else None
    if selected_count is None:
        for count in list(range(45, len(levels), 5)) + [len(levels)]:
            chosen, chosen_radiation, chosen_collisions = source_subset(
                'ar', 6, sorted(levels)[:count])
            balance = full_excitation_balance(chosen, chosen_radiation,
                                              chosen_collisions, temperature, DENSITY)
            errors = []
            for index, density in enumerate(DENSITY):
                reference = full[f'cooling_density_{index}']
                actual = balance[density][0]
                normal = reference > 64*np.finfo(float).tiny
                errors.append(float(np.max(np.abs(actual[normal]/reference[normal]-1))))
                if np.any((reference > 0) & (actual == 0)):
                    raise ValueError('Expanded atom lost positive full-source cooling')
            candidates.append(dict(levels=count,
                                   maximum_relative_cutoff_errors=errors))
            if max(errors) <= .003:
                selected_count = count
                break
    if selected_count is None:
        raise ValueError('No source cutoff candidate passed')
    directory = base/'native'/f'level{selected_count}'
    directory.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, str(ROOT/'tools/fitting/validate_source_excitation_atom.py'),
               '--element', 'ar', '--stage', '6', '--levels', str(selected_count),
               '--work-dir', str(directory), '--points', '801',
               '--full-source-profiles', str(profiles), '--driver', str(args.driver.resolve())]
    with (directory/'validation.log').open('w') as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=1500)
    native = json.loads((directory/'result.json').read_text())
    if not native['passed'] or not native['density_ratio_cooling']:
        raise ValueError('Native atom lacks measured density-ratio acceptance')
    result = dict(passed=True, legacy_levels=40, selected_levels=selected_count,
                  candidates=candidates, cutoff_target=.003,
                  source_profiles_sha256=digest(profiles),
                  native_result=str(directory/'result.json'),
                  native_result_sha256=digest(directory/'result.json'),
                  native_atom_sha256=native['product_sha256'],
                  driver_sha256=digest(args.driver))
    target = base/'native/result.json'
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
