#!/usr/bin/env python3
"""Check the thermal consumers of the optional kinetic-energy terms.

thermal_mod::net_rate (made public in a private copy of thermal_mod) must
equal the heating minus cooling assembled term by term from the public
consumers (check_optional_thermal_terms.f90), for C/O metals, helium, and
both, at 5000-20000 K and three radiation levels.  Across the inputs, the
ionization state and every consumer stay the same and net_rate moves by
exactly + charge-exchange heat (par%charge_exchange_heating) - metal
free-free (par%metal_freefree); both default on.
"""
import argparse
import json
from pathlib import Path

from run_checks import ROOT, compile_private_module, link, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    private = compile_private_module(ROOT/'src/thermal_mod.f90', work, [(
        '  public :: gas_thermal_update\n',
        '  public :: gas_thermal_update\n  public :: net_rate\n')])
    exe = link(work, 'tests/ion_stages/check_optional_thermal_terms.f90', [private])
    cases = {}
    for name, settings in (
            ('default', ''),
            ('both_off', ' par%charge_exchange_heating=.false.\n par%metal_freefree=.false.\n'),
            ('heat_only', ' par%metal_freefree=.false.\n'),
            ('freefree_only', ' par%charge_exchange_heating=.false.\n'),
            ('no_hydrogen_impact', ' par%neutral_hydrogen_excitation=.false.\n')):
        inp = work/(name+'.in')
        inp.write_text('&parameters\n'+settings+'/\n')
        log = run([exe, inp], log=work/(name+'.log'))
        if 'PASS: optional thermal consumers cases=27' not in log:
            raise RuntimeError(name+': incomplete')
        # k, T, radiation, n_e, net, base, cx, ff
        cases[name] = [list(map(float, line.split()[1:])) for line in log.splitlines()
                       if line.startswith('THERMAL ')]
    worst = 0.
    for i, off in enumerate(cases['both_off']):
        for name in ('default', 'heat_only', 'freefree_only'):
            row = cases[name][i]
            if row[:4] != off[:4] or row[5:] != off[5:]:
                raise RuntimeError(f'{name}: the state or a consumer moved with a thermal switch')
            delta = (row[6] if name in ('default', 'heat_only') else 0.) \
                - (row[7] if name in ('default', 'freefree_only') else 0.)
            error = abs((row[4]-off[4])-delta)/max(abs(off[4]), abs(delta), 1e-35)
            worst = max(worst, error)
    if worst > 1e-12:
        raise RuntimeError(f'net-rate differences differ from the switched terms by {worst}')
    result = {'passed': True, 'thermal_cases': 27*len(cases), 'choices': list(cases),
              'maximum_scaled_energy_difference_error': worst}
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
