"""Check the metal photoabsorption cross sections and their consumers.

Verner's phfit2 prescription, from the PHOTO / PHINNER / PHSHELL records of
each production element file: below the first inner-shell edge E_inn the
Verner et al. (1996) outer-shell fit,
  sigma = s0 [(x-1)^2 + yw^2] y^(P/2-5.5) (1+sqrt(y/ya))^-P,
  x = E/E0 - y0, y = sqrt(x^2 + y1^2);
at and above E_inn the sum over the occupied subshells with E >= E_th,s of
the Verner & Yakovlev (1995) fits,
  sigma_s = s0 [(y-1)^2 + yw^2] y^-(5.5+l-P/2) (1+sqrt(y/ya))^-P,  y = E/E0;
heating coefficient sum sigma_s (1 - E_th,s/E).  A transition with a PHTABLE
record (Ne III since 2026-10-11, Nahar 2019) takes that table, linear in
eps = E - E_th, below its last node, with heating sigma (1 - E_th/E); phfit2
at and above it.  species_cross_section must
match this evaluation for every transition of the 77-stage roster at 18
energies, and its packet, estimator and opacity consumers must carry the same
values (check_metal_photoabsorption.f90, with a private copy of species_mod
that exposes the transition order, stage fractions and estimators).
"""
import argparse
import json
from pathlib import Path

import numpy as np

from run_checks import ROOT, compile_private_module, link, run

ELEMENTS = ['c', 'n', 'o', 'ne', 's', 'ar', 'mg', 'fe', 'si', 'cl', 'ca']


def records(path):
    result = {}
    stage = None
    lines = iter(path.read_text().splitlines())
    for line in lines:
        x = line.split()
        if not x or x[0].startswith('#'):
            continue
        if x[0] == 'TRANSITION':
            stage = int(x[1])
            result[stage] = {'shells': []}
        if x[0] in ('ETH', 'PHINNER'):
            result[stage][x[0]] = float(x[1])
        if x[0] == 'PHOTO':
            result[stage]['PHOTO'] = list(map(float, x[1:]))
        if x[0] == 'PHSHELL':
            result[stage]['shells'].append(list(map(float, x[1:])))
        if x[0] == 'PHTABLE':
            rows = np.array([list(map(float, next(lines).split())) for _ in range(int(x[1]))])
            result[stage]['PHTABLE'] = (rows[:, 0], rows[:, 1]*1e-18)
    return result


def phfit2(r, energy):
    if energy < r['ETH']:
        return 0., 0.
    if 'PHTABLE' in r and energy - r['ETH'] < r['PHTABLE'][0][-1]:
        sigma = float(np.interp(energy - r['ETH'], *r['PHTABLE']))
        return sigma, sigma*(1-r['ETH']/energy)
    if energy < r['PHINNER']:
        eth, e0, s0, ya, p, yw, y0, y1 = r['PHOTO']
        x = energy/e0-y0
        y = np.sqrt(x*x+y1*y1)
        sigma = s0*((x-1)**2+yw*yw)*y**(p/2-5.5)*(1+np.sqrt(y/ya))**(-p)*1e-18
        return sigma, sigma*max(1-r['ETH']/energy, 0.)
    sigma = heat = 0.
    for l, eth, e0, s0, ya, p, yw in r['shells']:
        if energy < eth:
            continue
        y = energy/e0
        s = s0*((y-1)**2+yw*yw)*y**(-(5.5+l-p/2))*(1+np.sqrt(y/ya))**(-p)*1e-18
        sigma += s
        heat += s*max(1-eth/energy, 0.)
    return sigma, heat


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    private = compile_private_module(ROOT/'src/species_mod.f90', work, [(
        '  private\n',
        '  private\n  public :: flat_ie, flat_it, stage_frac, metal_rate_path\n')])
    exe = link(work, 'tests/ion_stages/check_metal_photoabsorption.f90', [private])
    output = work/'cross_sections.txt'
    run([exe, ROOT/'data/atomic', output], log=work/'native.log')
    measured = np.loadtxt(output)
    data = {i: records(ROOT/f'data/atomic/element_{el}.txt') for i, el in enumerate(ELEMENTS, 1)}
    transitions = sum(len(d) for d in data.values())
    if measured.shape != (transitions*18, 5) or not np.isfinite(measured).all():
        raise RuntimeError(f'incomplete output {measured.shape} for {transitions} transitions')
    worst = 0.
    inner = 0
    for ie, it, energy, sigma, heat in measured:
        r = data[int(ie)][int(it)]
        expected = np.array(phfit2(r, energy))
        native = np.array([sigma, heat])
        positive = expected > 0
        if np.any(native[~positive] != 0):
            raise RuntimeError(f'absorption where phfit2 has none: {ie} {it} {energy}')
        if positive.any():
            worst = max(worst, float(np.max(np.abs(native[positive]/expected[positive]-1))))
        inner += energy >= r['PHINNER']
    if worst > 1e-12:
        raise RuntimeError(f'species_cross_section differs from phfit2 by {worst}')
    result = dict(passed=True, transitions=transitions, samples=len(measured),
                  samples_at_or_above_inner_edge=int(inner), maximum_relative_error=worst)
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
