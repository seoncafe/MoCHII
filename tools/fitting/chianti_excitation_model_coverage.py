"""Stage coverage, manifest and source comparisons of the 'chianti' excitation model.

For the 77 stages I-VII of C N O Ne S Ar Mg Fe Si Cl Ca this program records
which excitation source supplies each stage (CHIANTI 11.0.2, Cloudy c25.00
Stout where CHIANTI lacks the ion, or none), the sizes of the written atoms,
the checks of each atom by validate_chianti_excitation_model.py, and writes
manifest.json (SHA-256 of every file) into the model directory. It also
compares the low-density cooling of the new fits with the fits of the
'mocassin_v2' model for the ions named in the October 11 audit (E3-E5;
the 'extended' model it was also compared with was removed on 2026-10-11,
and that comparison is recorded in md/atomic_stage_expansion/
chianti_model_2026-10-11/), and the CHIANTI Fe VII atom (Tayal & Zatsarinny 2014) with the
Stout Fe VII atom (Witthoeft & Badnell 2008) it replaces as the source.

Run: python3 chianti_excitation_model_coverage.py --model-dir DIR --checks DIR --output DIR
"""
import argparse
from datetime import datetime
from hashlib import sha256
import json
import os
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np

import chianti_cooling as C

ELEMENTS = ['c', 'n', 'o', 'ne', 's', 'ar', 'mg', 'fe', 'si', 'cl', 'ca']
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII']
NO_LINE_COOLING = {
    'c_7': 'bare nucleus; no bound electron',
    'ar_1': 'no electron-impact excitation data in CHIANTI 11.0.2 or Cloudy c25.00 Stout',
    'cl_1': 'no electron-impact excitation data in CHIANTI 11.0.2 or Cloudy c25.00 Stout',
    'ca_1': 'no electron-impact excitation data in CHIANTI 11.0.2 or Cloudy c25.00 Stout',
    'fe_1': 'no electron-impact excitation data in CHIANTI 11.0.2 or Cloudy c25.00 Stout',
    'cl_5': 'CHIANTI 11.0.2 has only .easplups (inner-shell excitation-autoionization); no valence collision data',
    'ca_3': 'CHIANTI 11.0.2 has only .easplups (inner-shell excitation-autoionization); no valence collision data',
}
STOUT_SOURCE = {
    'cl_7': 'Cloudy c25.00 Stout cl_7 (161 source levels, 32 bound outer-shell levels), Upsilon nodes 9800-9.8e7 K',
    'ar_2': 'Cloudy c25.00 Stout ar_2 ground fine-structure pair, Pelan & Berrington (1995) Table 3',
    'ca_4': 'Cloudy c25.00 Stout ca_4 ground fine-structure pair, Pelan & Berrington (1995) Table 3, from 1e4 K',
}
COMPARED_IONS = ['n_1', 'ca_2', 'c_1', 's_3', 'si_3', 'fe_2', 'fe_3', 'ca_7']


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def read_fit(path):
    records = [l.split() for l in Path(path).read_text().splitlines()
               if l.strip() and not l.lstrip().startswith('#')]
    n = int(records[0][0])
    data = np.array([[float(a), float(b)] for a, b in records[1:1+n]])
    return data[:, 0], data[:, 1]


def evaluate_fit(path, T):
    A, Ti = read_fit(path)
    return np.sum(A[:, None]*np.exp(-Ti[:, None]/T[None, :]), axis=0)/np.sqrt(T)


def chianti_files_present(element, stage):
    stem = Path(C.ion_dir(element, stage))/f'{element}_{stage}'
    return [e for e in ('elvlc', 'wgfa', 'scups', 'easplups')
            if Path(f'{stem}.{e}').exists()]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-dir', type=Path, required=True)
    parser.add_argument('--checks', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    model = args.model_dir
    stages = []
    for element in ELEMENTS:
        for stage in range(1, 8):
            tag = f'{element}_{stage}'
            atom = model/f'nlevel_{tag}.txt'
            fit = model/f'cooling_{tag}.txt'
            record = dict(stage=tag, ion=f'{element.capitalize()} {ROMAN[stage-1]}',
                          chianti_files=chianti_files_present(element, stage))
            if atom.exists():
                record['source'] = ('Stout' if tag in STOUT_SOURCE else 'CHIANTI 11.0.2')
                record['source_detail'] = STOUT_SOURCE.get(tag)
                record['nlevel_bytes'] = atom.stat().st_size
                record['cooling_bytes'] = fit.stat().st_size if fit.exists() else None
                check = args.checks/f'{tag}.json'
                if check.exists():
                    record['check'] = json.loads(check.read_text())
            else:
                record['source'] = 'none'
                record['reason'] = NO_LINE_COOLING[tag]
            stages.append(record)

    # Low-density cooling relative to the earlier models' fits.
    temperatures = np.array([5e3, 1e4, 2e4])
    root = model.parent
    comparison = {}
    for tag in COMPARED_IONS:
        new = evaluate_fit(model/f'cooling_{tag}.txt', temperatures)
        entry = dict(chianti=dict(zip(['5e3', '1e4', '2e4'], new.tolist())))
        for other in ('mocassin_v2',):
            path = root/other/f'cooling_{tag}.txt'
            if path.exists():
                entry[f'{other}_over_chianti'] = dict(
                    zip(['5e3', '1e4', '2e4'], (evaluate_fit(path, temperatures)/new).tolist()))
            else:
                entry[f'{other}_over_chianti'] = None
        comparison[tag] = entry

    # Fe VII: CHIANTI (installed) against the Stout atom.
    from stout_excitation_atoms import fevii_excitation_atom
    raw, radiation, transitions = fevii_excitation_atom()
    levels = {i: (g, e) for i, (g, e, _) in raw.items()}
    T = np.logspace(3, 5, 41)
    stout = C.cooling_low_density_from_levels(levels, radiation, transitions, T)
    chianti = C.cooling_low_density_limit('fe', 7, T)
    fevii = dict(temperature_K=T.tolist(), chianti=chianti.tolist(), stout=stout.tolist(),
                 stout_over_chianti_at={f'{t:.0e}': float(np.interp(np.log10(t), np.log10(T), stout/chianti))
                                        for t in (5e3, 1e4, 2e4, 5e4, 1e5)})

    files = sorted(p for p in model.iterdir() if p.suffix == '.txt')
    manifest = dict(
        model='chianti', generated_date=datetime.now(ZoneInfo('Asia/Seoul')).date().isoformat(),
        chianti_version=(Path(C.DBASE)/'VERSION').read_text().strip(),
        files=[dict(file=p.name, sha256=sha(p)) for p in files],
        stages_with_atom=[r['stage'] for r in stages if r['source'] != 'none'],
        stages_from_stout=sorted(STOUT_SOURCE),
        stages_without_line_cooling=NO_LINE_COOLING,
        atomic_accuracy_verified=False,
        generator_sha256={p.name: sha(p) for p in sorted(Path(__file__).parent.glob('*.py'))
                          if p.name in {'write_excitation_atom.py', 'chianti_cooling.py',
                                        'electric_dipole_probabilities.py', 'fit_cooling.py',
                                        'make_hcoll_fine_structure.py',
                                        'atomic_excitation_energies.py'}})
    (model/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'coverage.json').write_text(json.dumps(
        dict(stages=stages, cooling_comparison=comparison, fevii_chianti_vs_stout=fevii),
        indent=1)+'\n')
    lines = ['| Stage | Source | Levels | Radiative rows | Collision rows | nlevel size [kB] | '
             'Fit terms | Fit error [%] | Fortran vs source [%] |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in stages:
        c = r.get('check', {})
        if r['source'] == 'none':
            lines.append(f"| {r['ion']} | none: {r['reason']} | | | | | | | |")
            continue
        fit_error = c.get('fit_vs_source_max_relative')
        fortran = c.get('fortran_low_density_vs_source_max_relative')
        lines.append(
            f"| {r['ion']} | {r['source']} | {c.get('levels', '')} | {c.get('radiative_rows', '')} | "
            f"{c.get('collision_rows', '')} | {r['nlevel_bytes']/1024:.0f} | {c.get('fit_terms', '')} | "
            f"{'' if fit_error is None else f'{100*fit_error:.2f}'} | "
            f"{'' if fortran is None else f'{100*fortran:.1e}'} |")
    (args.output/'coverage_table.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
    print(json.dumps(comparison, indent=1))
    print(json.dumps(fevii['stout_over_chianti_at'], indent=1))


if __name__ == '__main__':
    main()
