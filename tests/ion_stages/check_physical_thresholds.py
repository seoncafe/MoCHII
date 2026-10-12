"""Check the loaded metal ionization thresholds against Cloudy's table.

The compiled potentials (metal_ionization_potential_mod) of stages I-VII of
the eleven elements are compared with Cloudy c25.00
data/ionization_potentials.dat, converted from cm^-1 with the exact SI
constants; the production element files must carry the same thresholds
(check_ionization_thresholds.f90).  A band reaching above the lowest
terminal potential must be refused, naming the first element in registry
order whose terminal stage lies below the band edge and the luminosity share
above it.
"""
import argparse
import json
from pathlib import Path
import re

from run_checks import ROOT, link, run

CLOUDY = Path.home()/'CLOUDY/c25.00/data/ionization_potentials.dat'
REGISTRY = [('c', 6), ('n', 7), ('o', 8), ('ne', 10), ('s', 16), ('ar', 18),
            ('mg', 12), ('fe', 26), ('si', 14), ('cl', 17), ('ca', 20)]
CM_TO_EV = 6.62607015e-34*299792458*100/1.602176634e-19
EEDGE = [13.598, 50.0, 200.0]    # the band of check_ionization_thresholds.f90
LUM = [1.0, 1.0]
REJECT_ENERGY = 130.0


def cloudy_potentials():
    expected = {}
    for line in CLOUDY.read_text().splitlines()[1:]:
        if line.lstrip().startswith('*'):
            break
        fields = line.split('#', 1)[0].split()
        if not fields:
            continue
        z, q = map(int, fields[:2])
        expected[z, q+1] = float(fields[2])*CM_TO_EV
    return expected


def luminosity_above(energy):
    total = 0.
    for k, l in enumerate(LUM):
        lo, hi = EEDGE[k], EEDGE[k+1]
        if lo >= energy:
            total += l
        elif hi > energy:
            total += l*(hi-energy)/(hi-lo)
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = link(work, 'tests/ion_stages/check_ionization_thresholds.f90')
    output = run([exe, ROOT/'data/atomic', work/'compiled_energies.txt'], log=work/'accept.log')
    transitions = int(re.search(r'PASS: (\d+) loaded transitions', output).group(1))
    expected = cloudy_potentials()
    error = 0.
    entries = 0
    for line in (work/'compiled_energies.txt').read_text().splitlines():
        z, s, value = line.split()
        key = int(z), int(s)
        value = float(value)
        if key[1] > key[0]:
            # A fully stripped ion has no ionization potential.
            if value < 1e300:
                raise RuntimeError(f'bare ion {key} permits further ionization')
        else:
            error = max(error, abs(value/expected[key]-1))
        entries += 1
    if entries != 77 or error > 1e-14:
        raise RuntimeError(f'compiled potentials differ from Cloudy: {entries} entries, {error}')
    if transitions != sum(min(7, z+1)-1 for _, z in REGISTRY):
        raise RuntimeError(f'{transitions} loaded transitions; the I-VII roster has 66')
    # The first element in registry order with a terminal potential below
    # the band edge is the one refused.
    first = next((name, z) for name, z in REGISTRY
                 if z+1 > 7 and expected[z, 7] < REJECT_ENERGY)
    terminal = expected[first[1], 7]
    share = luminosity_above(terminal)/luminosity_above(13.598)
    log = run([exe, ROOT/'data/atomic', work/'unused.txt', REJECT_ENERGY], log=work/'reject.log',
              expect_failure=f'the band reaches above the terminal stage of {first[0]} (stage 7')
    printed = float(re.search(r'eV; *([0-9.Ee+-]+) of the ionizing', log).group(1))
    if abs(printed/share-1) > 1e-3:
        raise RuntimeError(f'refusal reports share {printed}, expected {share}')
    result = dict(passed=True, compiled_potentials=entries, maximum_relative_difference=error,
                  loaded_transitions=transitions, refused_element=first[0],
                  refused_terminal_potential_eV=terminal, refused_share=share,
                  refused_share_printed=printed)
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
