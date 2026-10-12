#!/usr/bin/env python3
"""Check the q >= 5 Langevin charge exchange option on the production roster.

par%langevin_charge_exchange (default off) adds X^q + H -> X^(q-1) + H+ for
q >= 5 at the Ferland et al. (1997) eq. 2 polarization-capture rate
2 pi sqrt(2 C4/mu) with unit transfer probability, and the KF99 9.5-bohr
avoided-crossing kinetic release 2.86 (q-1) eV; both are evaluated here by
tools/fitting/charge_exchange_capture.py from the physical constants.  With
the option off those reactions carry no rate and no heat, and every q <= 4
reaction is the same either way.  The KF96 rates added at q = 3-4 must follow
the published fits (tools/fitting/charge_exchange_expansion.py) with the
KF99 mean heat.  Both runs must conserve charge reaction by reaction and close
the heat (check_added_charge_transfer_heat.f90).  The former name
par%high_charge_exchange must be refused with the name of its replacement.
"""
import argparse
import json
import math
from pathlib import Path
import sys

from run_checks import ROOT, link, run

sys.path.insert(0, str(ROOT/'tools/fitting'))
from charge_exchange_capture import avoided_crossing_mean_heat, langevin_recombination_rates  # noqa: E402
from charge_exchange_expansion import kf96_added_rates  # noqa: E402

STAGES = 77


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = link(work, 'tests/ion_stages/check_added_charge_transfer_heat.f90')
    results = {}
    for name, setting in (('default', ''), ('langevin', ' par%langevin_charge_exchange=.true.\n')):
        inp = work/(name+'.in')
        inp.write_text('&parameters\n'+setting+'/\n')
        output = work/(name+'.txt')
        log = run([exe, ROOT/'data/atomic', output, STAGES, inp], log=work/(name+'.log'))
        if 'PASS: staged ions=' not in log:
            raise RuntimeError(name+': no PASS line')
        rows = {}
        for line in output.read_text().splitlines():
            f = line.split()
            # element, transition, state; T, k_ion, k_rec, dE_ion, dE_rec, split, f(it+1)
            rows[(f[0], int(f[1]), int(f[2]))] = list(map(float, f[3:]))
        if len(rows) != 3*(STAGES-11):
            raise RuntimeError(f'{name}: {len(rows)} reaction samples')
        results[name] = rows
    capture = langevin_recombination_rates()
    worst_capture = 0.
    langevin_reactions = 0
    for key, default in results['default'].items():
        enabled = results['langevin'][key]
        element, charge, state = key
        if not all(math.isfinite(v) for v in default+enabled):
            raise RuntimeError(f'nonfinite sample {key}')
        if charge >= 5:
            if default[1:6] != [0.]*5:
                raise RuntimeError(f'q >= 5 reaction active with the option off: {key} {default}')
            expected = capture[element, charge]
            worst_capture = max(worst_capture, abs(enabled[2]/expected-1))
            if enabled[1] != 0 or abs(enabled[4]-avoided_crossing_mean_heat(charge)) > 1e-12:
                raise RuntimeError(f'Langevin reaction {key}: {enabled}')
            langevin_reactions += 1
        elif any(abs(x-y) > 1e-12*max(abs(x), abs(y), 1e-30)
                 for x, y in zip(default[:6], enabled[:6])):
            raise RuntimeError(f'q <= 4 reaction changed with the option: {key}')
    if langevin_reactions != 3*len(capture) or worst_capture > 5e-7:
        raise RuntimeError(f'{langevin_reactions} Langevin samples, worst {worst_capture}')
    worst_added = 0.
    for (element, charge), reaction in kf96_added_rates().items():
        for state in range(1, 4):
            row = results['default'][element, charge, state]
            worst_added = max(worst_added, abs(row[2]/reaction.rate(row[0])-1))
            if abs(row[4]-reaction.cloudy_energy_defect_ev) > 1e-12:
                raise RuntimeError(f'KF99 heat of {element} q={charge}: {row[4]}')
    if worst_added > 1e-5:
        raise RuntimeError(f'KF96 added rates differ by {worst_added}')
    renamed = work/'renamed.in'
    renamed.write_text('&parameters\n par%high_charge_exchange=.true.\n/\n')
    run([exe, ROOT/'data/atomic', work/'renamed.txt', STAGES, renamed], log=work/'renamed.log',
        expect_failure='use par%langevin_charge_exchange')
    result = dict(passed=True, samples_each_option=3*(STAGES-11),
                  langevin_samples=langevin_reactions,
                  maximum_capture_rate_relative_error=worst_capture,
                  kf96_added_reactions=len(kf96_added_rates()),
                  maximum_kf96_added_rate_relative_error=worst_added,
                  former_name_refused=True)
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
