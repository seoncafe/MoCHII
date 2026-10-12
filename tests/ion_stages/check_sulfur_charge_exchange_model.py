"""Check the term-resolved sulfur charge exchange against detailed balance.

element_s.txt carries S + H+ -> S+(term j) + H by product term (Zhao et al.
2005, Table III; CXSTATE / CXSTATET / CXSTATEK).  The reverse from each term
is its detailed balance, k_r,j = k_f,j g_react/(2 g_j) exp(-(dIP - E_j)/kT),
weighted by the population of that term, from the S II atom of the selected
excitation model.  For both models the species code's reverse rate and mean
kinetic energies must equal that evaluation, made here from the element
file and the driver's direct term populations; the CXSTATEK level ranges
must name levels whose weights and mean energies are the terms'.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from run_checks import ROOT, link, run

KB, EV, CM_TO_EV = 1.380649e-16, 1.602176634e-12, 1.239841984e-4


def records(path):
    return [s.split() for s in path.read_text().splitlines() if s.strip() and not s.startswith('#')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = link(work, 'tests/ion_stages/check_sulfur_charge_exchange_model.f90')
    element = records(ROOT/'data/atomic/element_s.txt')
    state = next(r for r in element if r[0] == 'CXSTATE')
    greact, dip = map(float, state[3:5])
    tt = np.array(list(map(float, next(r for r in element if r[0] == 'CXSTATET')[1:])))
    terms = [r for r in element if r[0] == 'CXSTATEK']
    g = np.array([float(r[4]) for r in terms])
    energy = np.array([float(r[5]) for r in terms])
    forward = np.array([list(map(float, r[6:])) for r in terms])
    models, checks = {}, {}
    for model in ('chianti', 'mocassin_v2'):
        output = work/(model+'.txt')
        run([exe, ROOT/'data/atomic', model, output], log=work/(model+'.log'), timeout=1800)
        data = np.loadtxt(output)
        if data.shape != (404, 10) or not np.isfinite(data).all():
            raise RuntimeError(f'{model}: incomplete output {data.shape}')
        atom = records(ROOT/'data/atomic/excitation'/model/'nlevel_s_2.txt')
        start = next(i for i, r in enumerate(atom) if r[0] == 'NLEV')
        nlev = int(atom[start][1])
        levels = np.array([[float(r[1]), float(r[2])] for r in atom[start+1:start+1+nlev]])
        mismatch = []
        for r in terms:
            lo, hi = int(r[2]), int(r[3])
            part = levels[lo-1:hi]
            weight = part[:, 0].sum()
            e = (part[:, 0]*part[:, 1]).sum()/weight*CM_TO_EV
            if abs(weight-float(r[4])) > 1e-6 or abs(e-float(r[5])) > .002:
                raise RuntimeError(f'{model}: term {r[1]} is not levels {lo}-{hi}')
            mismatch.append(abs(e-float(r[5])))
        t = data[:, 0]
        kf = np.array([np.exp(np.interp(np.log(t), np.log(tt), np.log(row))) for row in forward]).T
        de = dip-energy
        kr = kf*greact/(2*g)*np.exp(-de*EV/(KB*t[:, None]))
        population = data[:, 7:10]
        expected_reverse = (population*kr).sum(axis=1)
        expected_ion = (kf*de).sum(axis=1)/kf.sum(axis=1)
        expected_rec = -(population*kr*de).sum(axis=1)/expected_reverse
        rate_error = float(max(abs(data[:, 3]/expected_reverse-1)))
        ion_error = float(max(abs(data[:, 4]-expected_ion)))
        rec_error = float(max(abs(data[:, 5]-expected_rec)))
        forward_error = float(max(abs(data[:, 2]/kf.sum(axis=1)-1)))
        checks[model] = dict(samples=404, nlev=nlev, maximum_term_energy_mismatch_eV=max(mismatch),
                             maximum_forward_rate_relative_error=forward_error,
                             maximum_reverse_rate_relative_error=rate_error,
                             maximum_forward_mean_energy_error_eV=ion_error,
                             maximum_reverse_mean_energy_error_eV=rec_error,
                             minimum_three_term_population_sum=float(population.sum(axis=1).min()))
        if forward_error > 1e-8 or rate_error > 1e-8 or ion_error > 1e-10 or rec_error > 1e-8:
            raise RuntimeError(f'{model}: {checks[model]}')
        models[model] = data
    ratio = models['chianti'][:, 3]/models['mocassin_v2'][:, 3]
    result = dict(passed=True, samples=808, models=checks,
                  reverse_rate_ratio_chianti_to_mocassin_v2=[float(ratio.min()), float(ratio.max())])
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
