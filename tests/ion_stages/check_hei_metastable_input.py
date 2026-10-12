"""Check the He I 2^3S input defaults and run the coupled-physics gate.

setup_mod::read_parameters_with_implicit_defaults: par%hei_metastable_coupling
defaults on; add_fuv, left out of the input, is on when something absorbs
below 13.6 eV -- here the 2^3S population (coupling or the hei_metastable
diagnostic, with helium) -- and efuv_min then defaults to the 4.78 eV edge of
the 2^3S photoionization cross section (whether or not add_fuv is then
switched off in the input), else 6.0 eV; a value the input gives
is kept.  With the coupling off the ion solve must not depend on the leaf
(check_hei_metastable_input.f90).  The physics of the coupling is the gate
tests/hei_metastable_coupling/check_hei_metastable_coupling.f90.
"""
import argparse
import json
from pathlib import Path

from run_checks import ROOT, link, run

# name, namelist lines, expected (coupling, hei_metastable, add_fuv), efuv_min
CASES = [
    ('default', '', (True, False, True), 4.78),
    ('enabled', ' par%hei_metastable_coupling=.true.\n', (True, False, True), 4.78),
    ('disabled', ' par%hei_metastable_coupling=.false.\n', (False, False, False), 6.0),
    ('diagnostic', ' par%hei_metastable=.true.\n', (True, True, True), 4.78),
    ('explicit_fuv', ' par%add_fuv=.true.\n', (True, False, True), 4.78),
    ('explicit_no_fuv', ' par%add_fuv=.false.\n', (True, False, False), 4.78),
    ('explicit_edge', ' par%efuv_min=5.5\n', (True, False, True), 5.5),
    ('no_helium', ' par%abund_He=0\n', (True, False, False), 6.0),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = link(work, 'tests/ion_stages/check_hei_metastable_input.f90')
    physics = ROOT/'tests/hei_metastable_coupling/check_hei_metastable_coupling.f90'
    gate = link(work, physics.relative_to(ROOT))
    results = []
    for name, settings, expected, efuv in CASES:
        inp = work/(name+'.in')
        inp.write_text('&parameters\n'+settings+'/\n')
        log = run([exe, inp], log=work/(name+'.log'))
        rows = [row.split()[1:] for row in log.splitlines() if row.startswith('INPUT ')]
        if len(rows) != 1:
            raise RuntimeError(name+': missing input summary')
        flags = tuple(v == 'T' for v in rows[0][:6])
        # helium_charge_exchange, metal_freefree, charge_exchange_heating default on
        if flags != (*expected, True, True, True):
            raise RuntimeError(f'{name}: read {flags}, expected {expected} + (T, T, T)')
        if float(rows[0][6]) != efuv:
            raise RuntimeError(f'{name}: efuv_min {rows[0][6]}, expected {efuv}')
        if not expected[0] and 'PASS: disabled metastable population feedback' not in log:
            raise RuntimeError(name+': disabled-coupling solve check missing')
        results.append(dict(name=name, flags=flags, efuv_min_eV=efuv))
    log = run([gate], cwd=physics.parent, log=work/'coupled_physics.log')
    if 'PASS' not in log or 'FAIL' in log:
        raise RuntimeError('coupled physics gate failed')
    result = dict(passed=True, input_cases=len(CASES), cases=results, coupled_physics_passed=True)
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
