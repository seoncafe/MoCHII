"""Measure active ion stages and selected excitation coverage against v2.

This read-only inventory does not certify atomic accuracy or numerical
consumer validation. Source file availability does not imply invocation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

from mocassin_v2_atom_inventory import ELEMENTS, inventory

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(path):
    return [s.split() for s in path.read_text().splitlines()
            if s.strip() and not s.lstrip().startswith('#')]


def atom_counts(path):
    rows = records(path)
    result = {}
    levels = {}
    invalid = []
    for i, row in enumerate(rows):
        if row[0] not in ('NLEV', 'NRAD', 'NUPS'):
            continue
        n = int(row[1])
        result[row[0]] = n
        if row[0] == 'NLEV':
            levels = {int(r[0]): float(r[2]) for r in rows[i+1:i+1+n]}
            if sorted(levels) != list(range(1, n+1)):
                raise ValueError(f'Incomplete level index set: {path}')
        elif row[0] == 'NRAD':
            for r in rows[i+1:i+1+n]:
                lower, upper, rate = int(r[0]), int(r[1]), float(r[2])
                gap = levels[upper]-levels[lower]
                if not math.isfinite(rate) or rate < 0 or (rate > 0 and gap <= 0):
                    invalid.append([lower, upper, rate, gap])
    result['invalid_radiation'] = invalid
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=['chianti', 'mocassin_v2'], default='chianti')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve earlier coverage measurements')
    selected = ROOT/'data/atomic/excitation'/args.model
    manifest_path = selected/'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    hashes = {str(manifest_path): sha(manifest_path), str(Path(__file__)): sha(Path(__file__))}
    for row in manifest['files']:
        path = selected/row['file']
        if sha(path) != row['sha256']:
            raise ValueError(f'Selected atom/fit differs from its manifest: {path}')
        hashes[str(path)] = row['sha256']
    source = inventory()
    hashes.update(source['source_sha256'])
    stages = {}
    for element, z in ELEMENTS:
        path = ROOT/'data/atomic'/f'element_{element}.txt'
        rows = records(path)
        declaration = next(r for r in rows if r[0] == 'ELEMENT')
        count = int(declaration[3])
        assert declaration[1] == element and int(declaration[2]) == z
        assert [int(r[1]) for r in rows if r[0] == 'TRANSITION'] == list(range(1, count))
        stages[element] = count
        hashes[str(path)] = sha(path)
    ions = []
    for ion in source['ions']:
        element, stage = ion['element'], ion['stage']
        atom = selected/f'nlevel_{element}_{stage}.txt'
        fit = selected/f'cooling_{element}_{stage}.txt'
        installed = atom_counts(atom) if atom.exists() else {}
        ref = ion.get('atom', {})
        ions.append(dict(element=element, stage=stage, active=stage <= stages[element],
            reference_atom_available=ion['line_file_available'],
            reference_levels=ref.get('nlev'), reference_positive_radiation=ref.get('positive_radiation_pair_count'),
            reference_positive_collision_pairs=ref.get('positive_collision_pair_count'),
            selected_atom_available=atom.exists(), selected_fit_available=fit.exists(),
            selected_counts=installed, source_sha256=ion.get('source_sha256')))
    missing = [f"{r['element']}_{r['stage']}" for r in ions
               if r['reference_atom_available'] and not r['selected_atom_available']]
    insufficient = [f"{r['element']}_{r['stage']}" for r in ions
                    if r['reference_atom_available'] and r['selected_atom_available'] and
                    (r['selected_counts']['NLEV'] < r['reference_levels'] or
                     r['selected_counts']['NRAD'] < r['reference_positive_radiation'] or
                     r['selected_counts']['NUPS'] < r['reference_positive_collision_pairs'])]
    result = dict(model=args.model, production_modified=False, production_stage_counts=stages,
        active_metal_stages=sum(stages.values()), active_transitions=sum(n-1 for n in stages.values()),
        target_metal_stages=77, target_transitions=66, missing_reference_atoms=missing,
        insufficient_reference_row_counts=insufficient,
        invalid_radiation={f"{r['element']}_{r['stage']}": r['selected_counts']['invalid_radiation']
                           for r in ions if r['selected_counts'].get('invalid_radiation')},
        atomic_accuracy_verified=False, numerical_consumers_verified=False,
        scope='Current source mapping, selected manifest hashes, row counts, positive photon gaps, and active stage declarations only',
        ions=ions, source_sha256=hashes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('model','active_metal_stages','active_transitions',
                                          'missing_reference_atoms','insufficient_reference_row_counts','invalid_radiation')}, indent=2))


if __name__ == '__main__':
    main()
