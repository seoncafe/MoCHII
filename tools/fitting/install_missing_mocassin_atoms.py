"""Install previously validated, absent finite-level atoms without overwrites.

This does not activate expanded rate chains, change shared-ion level selection,
install exponential fits, or certify published atomic accuracy.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work-dir', type=Path, required=True)
    work = p.parse_args().work_dir.resolve()
    if work.exists():
        raise FileExistsError('Preserve previous installation records')
    source = ROOT/'md/atomic_stage_expansion/phase1_77_mocassin_v2/candidate_atoms'
    records = json.loads((source/'manifest.json').read_text())['atoms']
    review = ROOT/'md/atomic_stage_expansion/validation_2026-10-11/mocassin_source_atoms'
    checks = json.loads((review/'native/checks_result.json').read_text())
    if not checks['passed']:
        raise ValueError('Earlier native source checks failed')
    accepted = {(r['element'], r['stage']): r for r in checks['atoms']}
    missing = []
    for r in records:
        src = Path(r['path'])
        destination = ROOT/'data/atomic'/src.name
        if destination.exists():
            continue
        if hashlib.sha256(src.read_bytes()).hexdigest() != r['native_sha256']:
            raise ValueError('Candidate changed since preparation')
        if accepted[(r['element'], r['stage'])]['native_sha256'] != r['native_sha256']:
            raise ValueError('Candidate does not match accepted native checks')
        missing.append((r, src, destination))
    if len(missing) != 23:
        raise ValueError('Installation scope changed; inspect before writing')
    work.mkdir(parents=True)
    result = dict(passed=False, atoms=[], overwrites=0, production_stages_changed=False,
                  atomic_accuracy_verified=False)
    for r, src, destination in missing:
        # Exclusive creation also rejects a concurrent installation.
        with destination.open('xb') as output, src.open('rb') as input_file:
            shutil.copyfileobj(input_file, output)
        h = hashlib.sha256(destination.read_bytes()).hexdigest()
        if h != r['native_sha256']:
            raise ValueError('Installed bytes differ')
        result['atoms'].append({**r, 'installed_path': str(destination), 'installed_sha256': h})
        (work/'progress.json').write_text(json.dumps(result, indent=2)+'\n')
    result.update(passed=True, installed_atoms=len(result['atoms']))
    result['source_sha256'] = {str(f): hashlib.sha256(f.read_bytes()).hexdigest()
                              for f in (Path(__file__), source/'manifest.json', review/'native/checks_result.json')}
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: installed 23 absent finite-level atoms without overwrites', flush=True)


if __name__ == '__main__':
    main()
