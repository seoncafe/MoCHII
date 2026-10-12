"""Check the tabulated KF96 charge-transfer rates through the linked species code.

Kingdon & Ferland (1996) Table 1 gives X^q + H -> X^(q-1) + H+ as
k = a 1e-9 T4^b [1 + c exp(d T4)] on a stated temperature interval, outside
which the fit is held at its end value.  The production element files carry
those rates as CXRTAB tables (log-log in T); species_cx_coefficients must
reproduce the published fit, coefficients from
tools/fitting/charge_exchange_expansion.py, to the table's interpolation
accuracy at and beyond both ends.  N + H+ and N+ + H use KF96 Table 3 / Table
1 with the 100-50000 K interval and the 10860 K endothermic factor.  Three
malformed CXRTAB rows must be refused.
"""
import argparse
import json
import math
from pathlib import Path
import shutil
import sys

from run_checks import ROOT, link, run

sys.path.insert(0, str(ROOT/'tools/fitting'))
from charge_exchange_expansion import kf96_existing_rates  # noqa: E402

REGISTRY = ['c', 'n', 'o', 'ne', 's', 'ar', 'mg', 'fe', 'si', 'cl', 'ca']


def atomic_copy(directory):
    """Element files copied, the excitation models linked."""
    directory.mkdir(parents=True)
    for source in (ROOT/'data/atomic').glob('element_*.txt'):
        shutil.copy2(source, directory/source.name)
    (directory/'excitation').symlink_to(ROOT/'data/atomic/excitation')
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    index = {el: i+1 for i, el in enumerate(REGISTRY)}
    samples = []
    for (el, it), reaction in kf96_existing_rates().items():
        low, high = reaction.min_temperature, reaction.max_temperature
        temperatures = [low/10]+[low*(high/low)**(j/2000) for j in range(2001)]+[high*10]
        for t in temperatures:
            samples.append((index[el], it, t, 0., reaction.rate(t), el+str(it)))
    for t in [30., 100., 300., 1000., 2000., 3000., 5000., 10000., 50000., 100000.]:
        t4 = min(max(t, 100.), 50000.)/1e4
        base = t4**(-.29)*(1-.92*math.exp(-8.38*t4))
        samples.append((index['n'], 1, t, 4.55e-12*base*math.exp(-10860./t), 1.01e-12*base, 'n1'))
    (work/'sample_temperatures.txt').write_text(''.join(
        f'{ie} {it} {t:.17e}\n' for ie, it, t, ki, kr, key in samples))
    exe = link(work, 'tests/ion_stages/check_charge_transfer_ranges.f90')
    run([exe, ROOT/'data/atomic', work/'sample_temperatures.txt', work/'compiled_rates.txt'],
        log=work/'rates.log')
    output = (work/'compiled_rates.txt').read_text().splitlines()
    if len(output) != len(samples):
        raise RuntimeError('incomplete compiled samples')
    errors = {}
    for line, sample in zip(output, samples):
        ie, it, t, ki, kr, key = sample
        fields = line.split()
        values = list(map(float, fields[2:]))
        if int(fields[0]) != ie or int(fields[1]) != it or values[0] != t:
            raise RuntimeError('changed sample identity')
        for measured, expected in zip(values[1:], (ki, kr)):
            if not math.isfinite(measured):
                raise RuntimeError('nonfinite compiled rate')
            error = abs(measured/expected-1) if expected else abs(measured)
            errors[key] = max(errors.get(key, 0.), error)
    if max(errors.values()) > 1e-5 or errors['n1'] > 1e-13:
        raise RuntimeError('compiled charge-transfer fit mismatch: '+str(errors))
    malformed = {'negative': 'CXRTAB 2 1000 2000 -1e-9 1e-9',
                 'reversed': 'CXRTAB 2 2000 1000 1e-9 1e-9',
                 'nan': 'CXRTAB 2 1000 2000 NaN 1e-9'}
    for key, row in malformed.items():
        directory = atomic_copy(work/key)
        path = directory/'element_c.txt'
        lines = path.read_text().splitlines()
        pos = next(i for i, line in enumerate(lines) if line.startswith('CXRTAB'))
        lines[pos] = row
        path.write_text('\n'.join(lines)+'\n')
        run([exe, directory, work/'sample_temperatures.txt', directory/'unused.txt'],
            log=directory/'rejection.log', expect_failure='CXRTAB')
    result = dict(passed=True, samples=len(samples), reactions=len(errors),
                  maximum_relative_errors=errors, malformed_tables_rejected=len(malformed))
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
