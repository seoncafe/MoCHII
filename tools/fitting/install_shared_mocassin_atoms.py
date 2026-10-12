"""Expand six existing atoms with cooling fits from the same source model.

Preserve original production files and require compatibility of every retained
level's weight, energy, and positive radiative transition. Fe II is excluded
because its forbidden radiation requires a separate physical decision.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

from mocassin_v2_atom_inventory import MOCASSIN, read_atom

ROOT = Path(__file__).resolve().parents[2]
TARGETS = {('c', 3): 5, ('n', 3): 5, ('o', 1): 5,
           ('ne', 3): 5, ('s', 4): 5, ('si', 2): 5}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--prepare-only', action='store_true',
                        help='Preserve and inspect current files; select candidates only in isolated diagnostics')
    parser.add_argument('--matching-cooling-review', type=Path,
                        help='Verified six matching atom/fit pairs; required for production replacement')
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if work.exists():
        raise FileExistsError('Preserve prior installation records')
    if not args.prepare_only and args.matching_cooling_review is None:
        raise ValueError('Atom-only replacement would mix thermal and radiative source data')
    paired = {}
    if args.matching_cooling_review:
        review = json.loads(args.matching_cooling_review.read_text())
        if not review['passed'] or review['selected_cooling_fits'] != 6:
            raise ValueError('Incomplete matching-source review')
        paired = {(r['element'], r['stage']): r for r in review['atoms']}
        if set(paired) != set(TARGETS):
            raise ValueError('Matching cooling scope differs')
    manifest = ROOT/'md/atomic_stage_expansion/phase1_77_mocassin_v2/candidate_atoms/manifest.json'
    native_checks = ROOT/'md/atomic_stage_expansion/validation_2026-10-11/mocassin_source_atoms/native/checks_result.json'
    checks = json.loads(native_checks.read_text())
    if not checks['passed']:
        raise ValueError('Native candidate validation failed')
    accepted = {(r['element'], r['stage']): r for r in checks['atoms']}
    changes = []
    for record in json.loads(manifest.read_text())['atoms']:
        key = (record['element'], record['stage'])
        if key not in TARGETS:
            continue
        candidate = Path(record['path'])
        destination = ROOT/'data/atomic'/candidate.name
        digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
        if digest != record['native_sha256'] or digest != accepted[key]['native_sha256']:
            raise ValueError('Candidate changed after numerical validation')
        if paired:
            pair = paired[key]
            fit = Path(pair['selected_cooling_path'])
            old_fit = ROOT/'data/atomic'/fit.name
            if digest != pair['source_atom_sha256']:
                raise ValueError('Matching cooling uses a different atom')
            if hashlib.sha256(fit.read_bytes()).hexdigest() != pair['cooling_sha256']:
                raise ValueError('Matching cooling bytes changed')
            if hashlib.sha256(old_fit.read_bytes()).hexdigest() != pair['previous_cooling_sha256']:
                raise ValueError('Production cooling changed after review')
        external = MOCASSIN/record['relative_file']
        if hashlib.sha256(external.read_bytes()).hexdigest() != record['source_sha256']:
            raise ValueError('External source changed')
        source = read_atom(external, include_tables=True)
        rows = [s.split() for s in destination.read_text().splitlines()
                if s.strip() and not s.lstrip().startswith('#')]
        first = next(i for i, s in enumerate(rows) if s[0] == 'NLEV')
        old_count = int(rows[first][1])
        if old_count != TARGETS[key]:
            raise ValueError('Installed model changed; inspect before replacing')
        maximum_energy_difference = 0.
        for old, new in zip(rows[first+1:first+1+old_count], source['levels']):
            if int(old[0]) != new['index'] or float(old[1]) != new['statistical_weight']:
                raise ValueError('Retained statistical weights or ordering differ')
            maximum_energy_difference = max(maximum_energy_difference,
                abs(float(old[2])-new['energy_cm_inverse']))
        if maximum_energy_difference > 0.05:
            raise ValueError('Retained energy mismatch exceeds source precision')
        start = next(i for i, s in enumerate(rows) if s[0] == 'NRAD')
        old_radiation = {(int(s[0]), int(s[1])): float(s[2])
                        for s in rows[start+1:start+1+int(rows[start][1])] if float(s[2]) > 0}
        new_radiation = {(s['lower'], s['upper']): s['A_s_inverse']
                        for s in source['spontaneous_probabilities']
                        if s['upper'] <= old_count and s['A_s_inverse'] > 0}
        if old_radiation != new_radiation:
            raise ValueError('Retained positive radiative transitions differ')
        changes.append((record, candidate, destination, dict(
            previous_nlev=old_count, maximum_retained_energy_difference_cm_inverse=maximum_energy_difference,
            retained_positive_radiative_pairs=len(old_radiation),
            previous_sha256=hashlib.sha256(destination.read_bytes()).hexdigest())))
    if len(changes) != 6:
        raise ValueError('Incomplete six-ion scope')
    work.mkdir(parents=True)
    previous = work/'previous_atoms'
    previous.mkdir()
    previous_cooling = work/'previous_cooling'
    if paired:
        previous_cooling.mkdir()
    result = dict(passed=False, atoms=[], expanded_atoms=0,
                  atomic_accuracy_verified=False, production_stages_changed=False,
                  production_installed=not args.prepare_only)
    # Preserve every affected original before replacing any production file.
    for record, candidate, destination, compatibility in changes:
        preserved = previous/destination.name
        shutil.copyfile(destination, preserved)
        if hashlib.sha256(preserved.read_bytes()).hexdigest() != compatibility['previous_sha256']:
            raise ValueError('Original atom preservation failed')
        if paired:
            pair = paired[record['element'], record['stage']]
            name = Path(pair['selected_cooling_path']).name
            preserved_fit = previous_cooling/name
            shutil.copyfile(ROOT/'data/atomic'/name, preserved_fit)
            if hashlib.sha256(preserved_fit.read_bytes()).hexdigest() != pair['previous_cooling_sha256']:
                raise ValueError('Original cooling preservation failed')
    for record, candidate, destination, compatibility in changes:
        preserved = previous/destination.name
        if not args.prepare_only:
            shutil.copyfile(candidate, destination)
        selected = candidate if args.prepare_only else destination
        if hashlib.sha256(selected.read_bytes()).hexdigest() != record['native_sha256']:
            raise ValueError('Installed bytes differ')
        cooling_record = {}
        if paired:
            pair = paired[record['element'], record['stage']]
            fit = Path(pair['selected_cooling_path'])
            fit_destination = ROOT/'data/atomic'/fit.name
            if not args.prepare_only:
                shutil.copyfile(fit, fit_destination)
            selected_fit = fit if args.prepare_only else fit_destination
            if hashlib.sha256(selected_fit.read_bytes()).hexdigest() != pair['cooling_sha256']:
                raise ValueError('Installed cooling bytes differ')
            cooling_record = dict(selected_cooling_path=str(selected_fit),
                cooling_sha256=pair['cooling_sha256'],
                cooling_maximum_K=pair['cooling_maximum_K'],
                preserved_previous_cooling_path=str(previous_cooling/fit.name),
                previous_cooling_sha256=pair['previous_cooling_sha256'])
        result['atoms'].append({**record, **compatibility, **cooling_record,
            'selected_path': str(selected), 'production_path': str(destination),
            'preserved_previous_path': str(preserved)})
        (work/'progress.json').write_text(json.dumps(result, indent=2)+'\n')
    result.update(passed=True, expanded_atoms=6)
    result['selected_cooling_fits'] = 6 if paired else 0
    result['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in (Path(__file__), manifest, native_checks)}
    if args.matching_cooling_review:
        result['source_sha256'][str(args.matching_cooling_review.resolve())] = hashlib.sha256(
            args.matching_cooling_review.read_bytes()).hexdigest()
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(('Prepared six candidate atom/fit pairs without production changes.' if paired
           else 'Prepared six candidate atoms without production changes.') if args.prepare_only
          else 'Expanded six shared-ion atoms and matching cooling fits; twelve originals preserved.')


if __name__ == '__main__':
    main()
