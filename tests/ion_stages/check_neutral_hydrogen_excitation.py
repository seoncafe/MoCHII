#!/usr/bin/env python3
"""Check the atomic-hydrogen excitation consumer in both excitation models.

For carbon (C I, C II) and oxygen (O I), in each model directory,
check_neutral_hydrogen_excitation.f90 holds metal_cooling_H to n_ion times
the radiated-power increment of the joint electron + hydrogen level
equilibrium.  Here: the switch defaults on and then adds positive cooling in
neutral gas; off, it adds none; and the electron-impact cooling is the same
either way.
"""
import argparse
import json
from pathlib import Path

from run_checks import ROOT, link, run

HYDROGEN_STAGES = {'c': 2, 'o': 1}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = link(work, 'tests/ion_stages/check_neutral_hydrogen_excitation.f90')
    results = {}
    for model in ('chianti', 'mocassin_v2'):
        for element, stages in HYDROGEN_STAGES.items():
            rows = {}
            for name, setting in (('default', ''), ('off', ' par%neutral_hydrogen_excitation=.false.\n')):
                inp = work/f'{model}_{element}_{name}.in'
                inp.write_text(f"&parameters\n par%metal_excitation_model='{model}'\n"+setting+'/\n')
                log = run([exe, ROOT/'data/atomic', element, inp], log=inp.with_suffix('.log'))
                # stages with hydrogen rates, T, n_e, electron cooling, hydrogen increment
                rows[name] = [list(map(float, line.split()[1:])) for line in log.splitlines()
                              if line.startswith('COOLING')]
                if len(rows[name]) != 9 or any(int(r[0]) != stages for r in rows[name]):
                    raise RuntimeError(f'{model} {element} {name}: {rows[name]}')
            if not all(r[4] > 0 for r in rows['default']):
                raise RuntimeError(f'{model} {element}: no hydrogen-impact cooling in neutral gas')
            if any(r[4] != 0 for r in rows['off']):
                raise RuntimeError(f'{model} {element}: cooling with the switch off')
            if any(a[3] != b[3] for a, b in zip(rows['default'], rows['off'])):
                raise RuntimeError(f'{model} {element}: electron-impact cooling moved')
            results[f'{model}_{element}'] = rows
    result = {'passed': True, 'samples': 9*2*len(results), 'rows': results}
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'passed': True, 'samples': result['samples']}))


if __name__ == '__main__':
    main()
