"""Check the Fe recombination coefficients the stage balance actually uses.

check_fe_recombination.f90 reads alpha_q back from the stationary balance
of the production Fe file (photoionization and recombination only) over
q^2 (20 K - 1e7 K), the domain of the Zhang et al. (2026, ApJS 285, 50)
fits.  For q = 3-6 the expectation is the published ground-level RR (eqs.
5-6) + S1 DR fit, coefficients typed from the paper's Tables 6/7 in
tools/fitting/zhang_fe_recombination.py; for q = 1-2 it is the RR2 power
law A (T/1e4)^-eta plus the Badnell DR sum T^-1.5 sum c_i exp(-E_i/T) of
the element file's own records, evaluated here.
"""
import argparse
import json
import math
from pathlib import Path
import sys

from run_checks import ROOT, link, run

sys.path.insert(0, str(ROOT/'tools/fitting'))
from zhang_fe_recombination import dr_fit, rr_fit  # noqa: E402


def file_coefficients(path):
    records = {}
    stage = None
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith('#'):
            continue
        if fields[0] == 'TRANSITION':
            stage = int(fields[1])
            records[stage] = {}
        elif fields[0] in ('RR', 'RR2', 'DR', 'DR2'):
            records[stage][fields[0]] = list(map(float, fields[1:]))
    return records


def file_alpha(record, t):
    if 'RR2' in record:
        a, eta = record['RR2']
        alpha = a*(t/1e4)**(-eta)
    else:
        a, b, t0, t1, c, t2 = record['RR']
        bt = b+c*math.exp(-t2/t)
        x = math.sqrt(t/t0)
        alpha = a/(x*(1+x)**(1-bt)*(1+math.sqrt(t/t1))**(1+bt))
    n = int(record['DR'][0])
    coefficients, energies = record['DR'][1:1+n], record['DR'][1+n:1+2*n]
    alpha += math.fsum(c*math.exp(-e/t) for c, e in zip(coefficients, energies))/t**1.5
    if 'DR2' in record:
        a, b, t0, t1 = record['DR2']
        alpha += a*t**-1.5*math.exp(-t0/t)*(1+b*math.exp(-t1/t))
    return alpha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = link(work, 'tests/ion_stages/check_fe_recombination.f90')
    run([exe, ROOT/'data/atomic', work/'native_rates.txt'], log=work/'native.log')
    records = file_coefficients(ROOT/'data/atomic/element_fe.txt')
    worst = {}
    for line in (work/'native_rates.txt').read_text().splitlines():
        q, t, a = line.split()
        q, t, a = int(q), float(t), float(a)
        expected = rr_fit(q, t)+dr_fit(q, t) if q >= 3 else file_alpha(records[q], t)
        worst[q] = max(worst.get(q, 0.), abs(a/expected-1))
    if sorted(worst) != [1, 2, 3, 4, 5, 6] or max(worst.values()) > 1e-12:
        raise RuntimeError(f'Fe recombination differs: {worst}')
    result = dict(passed=True, samples=6*301, maximum_relative_error_by_charge=worst)
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
