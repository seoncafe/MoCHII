"""Measure installed and supplied v2 excitation models without changing data."""
import argparse
import hashlib
import json
from pathlib import Path
from mocassin_v2_atom_inventory import inventory

ROOT = Path(__file__).resolve().parents[2]


def counts(path):
    result = {}
    if not path.exists():
        return result
    with path.open() as stream:
        for line in stream:
            fields = line.split()
            if fields and fields[0] in ('NLEV', 'NRAD', 'NUPS'):
                result[fields[0]] = int(fields[1])
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work-dir', type=Path, required=True)
    a = p.parse_args()
    work = a.work_dir.resolve()
    if work.exists():
        raise FileExistsError('Preserve previous inventory')
    work.mkdir(parents=True)
    source = inventory()
    stages = {}
    hashes = dict(source['source_sha256'])
    for element in ('c','n','o','ne','s','ar','mg','fe','si','cl','ca'):
        path = ROOT/'data/atomic'/('element_'+element+'.txt')
        row = next(s.split() for s in path.read_text().splitlines() if s.startswith('ELEMENT '))
        stages[element] = int(row[3])
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    rows = []
    for ion in source['ions']:
        if not ion['line_file_available']:
            continue
        element, stage = ion['element'], ion['stage']
        atom = ROOT/'data/atomic'/f'nlevel_{element}_{stage}.txt'
        fit = ROOT/'data/atomic'/f'cooling_{element}_{stage}.txt'
        current = counts(atom)
        for file in (atom, fit):
            if file.exists():
                hashes[str(file)] = hashlib.sha256(file.read_bytes()).hexdigest()
        rows.append(dict(element=element, stage=stage,
            reference_levels=ion['atom']['nlev'],
            installed_levels=current.get('NLEV'),
            installed_radiative_rows=current.get('NRAD'),
            installed_collision_rows=current.get('NUPS'),
            production_stage_active=stage <= stages[element],
            source_temperature_maximum_K=max(ion['atom']['temperatures_K']),
            source_temperature_nodes=len(ion['atom']['temperatures_K']),
            reference_atomic_accuracy_verified=False,
            source_sha256=ion['source_sha256']))
    larger = [r for r in rows if (r['installed_levels'] or 0) > r['reference_levels']]
    result = dict(production_modified=False, source_models=len(rows),
        larger_installed_models=len(larger),
        larger_active_models=sum(r['production_stage_active'] for r in larger),
        production_metal_stages=sum(stages.values()),
        production_adjacent_transitions=sum(n-1 for n in stages.values()),
        larger_models=larger, rows=rows, source_sha256=hashes,
        scope='Direct source mapping and installed file counts; no state-identity, atomic-accuracy, or negligible-cooling claim')
    (work/'inventory.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f"Measured {len(larger)} larger installed shared-ion models; {result['larger_active_models']} are in active rate chains.")


if __name__ == '__main__':
    main()
