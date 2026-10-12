"""Install the approved 142-level Fe II atom and matching cooling together.

Historical record (2026-10-11 installation of the 142-level Fe II atom).
Do not rerun: it writes the root data/atomic/nlevel_fe_2.txt and
cooling_fe_2.txt, a layout the code no longer reads (metal atoms are under
data/atomic/excitation/<model>/, nlevel_mod::metal_excitation_data_file).
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'md/atomic_stage_expansion/validation_2026-10-11'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', required=True, type=Path)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if work.exists():
        raise FileExistsError('Preserve prior installation records')
    preparation = BASE / 'feii_radiative_cooling/prepare/result.json'
    native = BASE / 'feii_radiative_cooling/native/result.json'
    density = BASE / 'feii_density_cooling_checked/native/result.json'
    root_review = BASE / 'feii_density_cooling_checked/root_review.json'
    ready = json.loads(preparation.read_text())
    checks = json.loads(native.read_text())
    table = json.loads(density.read_text())
    review = json.loads(root_review.read_text())
    assert ready['passed'] and checks['passed'] and table['passed']
    assert review['source_and_artifact_hashes_verified']
    assert review['figure_pages_visually_inspected'] == 4
    atom = Path(ready['source_atom_path'])
    fit = Path(ready['fit_path'])
    assert digest(atom) == ready['source_atom_sha256'] == table['source_sha256'][str(atom)]
    assert digest(fit) == ready['fit_sha256'] == table['source_sha256'][str(fit)]
    destination_atom = ROOT / 'data/atomic/nlevel_fe_2.txt'
    destination_fit = ROOT / 'data/atomic/cooling_fe_2.txt'
    assert digest(destination_atom) == ready['source_sha256'][str(destination_atom)]
    # The candidate explicitly preserves these forty positive transitions.
    def rows(path):
        records = [s.split() for s in path.read_text().splitlines()
                   if s.strip() and not s.lstrip().startswith('#')]
        i = next(i for i, r in enumerate(records) if r[0] == 'NRAD')
        return {(int(r[0]), int(r[1])): float(r[2])
                for r in records[i+1:i+1+int(records[i][1])]}
    previous = rows(destination_atom)
    current = rows(atom)
    assert len(previous) == 40 and len(current) == 1184
    # The candidate construction maps state identities; its source report
    # verifies the index mapping before adding the forty forbidden rates.
    source_report = BASE / 'feii_forbidden_effect/products/result.json'
    assert json.loads(source_report.read_text())['passed']
    work.mkdir(parents=True)
    originals = work / 'previous'
    originals.mkdir()
    replacements = [(atom, destination_atom), (fit, destination_fit)]
    records = []
    for source, destination in replacements:
        saved = originals / destination.name
        shutil.copyfile(destination, saved)
        assert digest(saved) == digest(destination)
        records.append(dict(source_path=str(source), source_sha256=digest(source),
                            production_path=str(destination), previous_path=str(saved),
                            previous_sha256=digest(saved)))
    for source, destination in replacements:
        shutil.copyfile(source, destination)
        assert digest(destination) == digest(source)
    result = dict(passed=True, production_installed=True, user_approval_date='2026-10-11',
                  nlev=142, nrad=1184, retained_forbidden_probabilities=40,
                  nups=10011, matching_fit_terms=29, atomic_accuracy_verified=False,
                  records=records,
                  source_sha256={str(p): digest(p) for p in
                                 [Path(__file__), preparation, native, density, root_review, source_report]})
    (work / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Installed the approved 142-level Fe II atom and matching fit; both originals preserved.')


if __name__ == '__main__':
    main()
