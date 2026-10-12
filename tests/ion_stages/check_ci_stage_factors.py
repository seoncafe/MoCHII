#!/usr/bin/env python3
"""Compare the selected HYBRID stage factors with the actual Cloudy table."""
import argparse
from pathlib import Path
import re
import subprocess

from run_checks import fortran_command

ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path.home() / 'CLOUDY/c25.00/source/atmdat_adfa.cpp'
TARGETS = {1: 1, 2: 2, 6: 4, 7: 5, 8: 4, 10: 4, 12: 3,
           14: 4, 16: 6, 17: 6, 18: 6, 20: 5, 26: 6}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    path = parser.parse_args().work_dir.resolve()
    path.mkdir(parents=True, exist_ok=True)
    block = SOURCE.read_text().split('static const double DereRatio[30][30]', 1)[1].split('};', 1)[0]
    rows = [[float(v) for v in row.split(',')] for row in re.findall(r'\{([^{}]+)\}', block)]
    assert len(rows) == 30 and all(len(row) == 30 for row in rows)
    driver = '''program ionization_factors
  use define
  use recomb_mod, only : ci_dere_ratio
  implicit none
  integer :: z, stage
  par%ci_model = 'dere_hybrid'
  do z = 1, 26
     do stage = 1, 6
        write(*,'(2i4,es25.16)') z, stage, ci_dere_ratio(z,stage)
     end do
  end do
  par%ci_model = 'voronov'
  do z = 1, 26
     do stage = 1, 6
        if (ci_dere_ratio(z,stage) /= 1.0_wp) stop 1
     end do
  end do
end program ionization_factors
'''
    (path/'ionization_factors.f90').write_text(driver)
    objects = ['define', 'utility', 'recomb_mod', 'hydrogenic_recombination_mod']
    subprocess.run([*fortran_command(), '-I'+str(path), '-Isrc', '-module', str(path), '-o',
                    str(path/'check.x'), str(path/'ionization_factors.f90'),
                    *[str(ROOT/'src'/(name+'.o')) for name in objects]], cwd=ROOT, check=True)
    result = subprocess.run([str(path/'check.x')], check=True, text=True, capture_output=True)
    count = 0
    for line in result.stdout.splitlines():
        z, stage, value = line.split()
        z, stage, value = int(z), int(stage), float(value)
        if z in TARGETS and stage <= TARGETS[z]:
            assert abs(value-rows[z-1][stage-1]) < 1e-14, (z, stage, value)
            count += 1
    assert count == 56
    print(f'PASS: {count} H/He/metal HYBRID factors match Cloudy; raw Voronov factors remain one')


if __name__ == '__main__':
    main()
