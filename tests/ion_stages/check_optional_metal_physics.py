#!/usr/bin/env python3
"""Check metal electron-impact ionization in the stage balance and energy budget.

A synthetic seven-stage carbon (C VII bare) with analytic rates: photo-
ionization gamma = 2e-10 s^-1, RR alpha = A/[x (1+x)^2], x = sqrt(T/1e4 K)
(the Badnell form with B = 0.5, T0 = T1 = 1e4 K, C = 0), and Voronov
collisional ionization C = A e^-U/(1+U), U = E_th/kT.  At T = 1e5 K,
n_e = n_HII = 100 cm^-3, n_C = 1e-2 cm^-3 the printed fractions must satisfy
the stationary adjacent balance f(it+1)/f(it) = (gamma + n_e C)/(n_e alpha)
with and without electron-impact ionization, and the metal ionization-energy
loss must equal n_e n_C sum_it [f(it+1) beta_it + f(it) C_it E_th,it] with
beta = kT[(3/2) alpha + d alpha/d ln T], evaluated here analytically.  The
switch defaults on and is read from the namelist; a C VI terminal stage is
refused under a band reaching 1000 eV, bare C VII is not.
"""
import argparse
import json
import math
from pathlib import Path

from run_checks import link, run

ENERGIES = [11.2603, 24.3833, 47.8878, 64.4935, 392.0905, 489.9932]   # C I-VI [eV]
GAMMA, NE, NC, T = 2e-10, 100., 1e-2, 1e5
RR_A, CI_A = 1e-11, 1e-8
KB, EV = 1.380649e-16, 1.602176634e-12


def element_text(stages):
    rows = [f'ELEMENT c 6 {stages}']
    for i, e in enumerate(ENERGIES[:stages-1], 1):
        rows += [f'TRANSITION {i}', f'ETH {e}', f'PHOTO {e} 10 1 1 2 1 0 1', 'PHINNER 1000',
                 f'PHSHELL 0 {e} 10 1 1 2 1', f'CI {e} 0 {CI_A} 1 0',
                 f'RR {RR_A} 0.5 1e4 1e4 0 1', 'DR 0', 'CXI 0', 'CXR 0', 'CXHEAT 0 0 0']
    return '\n'.join(rows)+'\n'


def fixture(directory, stages):
    model = directory/'excitation/chianti'
    model.mkdir(parents=True, exist_ok=True)
    (model/'manifest.json').write_text('{"model": "chianti", "files": []}\n')
    (directory/'element_c.txt').write_text(element_text(stages))
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = link(work, 'tests/ion_stages/check_optional_metal_physics.f90')
    inp = work/'options.in'
    inp.write_text('&parameters\n par%metal_collisional_ionization=.false.\n/\n')
    log = run([exe, fixture(work/'bare', 7), inp], log=work/'native.log')
    states = {}
    for line in log.splitlines():
        values = line.split()
        if values and values[0] == 'METAL_STATE':
            states[int(values[1])] = [float(v) for v in values[2:]]
    if sorted(states) != [0, 1]:
        raise RuntimeError('missing printed states')
    x = math.sqrt(T/1e4)
    alpha = RR_A/(x*(1+x)**2)
    dalpha = alpha*(-0.5-x/(1+x))
    beta = KB*T*(1.5*alpha+dalpha)
    worst_balance = worst_loss = 0.
    for selection, row in states.items():
        frac, cool = row[:7], row[8]
        collisional = [CI_A*math.exp(-u)/(1+u) if selection else 0.
                       for u in (e*EV/(KB*T) for e in ENERGIES)]
        expected = 0.
        for it in range(6):
            ratio = (GAMMA+NE*collisional[it])/(NE*alpha)
            if frac[it] > 0:
                worst_balance = max(worst_balance, abs(frac[it+1]/(frac[it]*ratio)-1))
            expected += NE*NC*(frac[it+1]*beta+frac[it]*collisional[it]*ENERGIES[it]*EV)
        worst_loss = max(worst_loss, abs(cool/expected-1))
    if worst_balance > 1e-12 or worst_loss > 1e-7:
        raise RuntimeError(f'balance {worst_balance}, ionization-energy loss {worst_loss}')
    run([exe, fixture(work/'bound', 6), inp, 'bound'], log=work/'bound.log',
        expect_failure='the band reaches above the terminal stage of c (stage 6')
    result = dict(passed=True, states=states, maximum_stationary_ratio_error=worst_balance,
                  maximum_ionization_energy_loss_error=worst_loss,
                  bound_carbon_refused=True, bare_carbon_accepted=True)
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'states'}))


if __name__ == '__main__':
    main()
