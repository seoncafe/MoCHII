"""Exercise neutral metal initialization and the subsequent irradiated state."""
import argparse
import json
from pathlib import Path
from check_optional_metal_physics import fixture
from run_checks import link, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = link(work, 'tests/ion_stages/check_neutral_metal_initialization.f90')
    atoms = fixture(work/'atoms',7)
    run(['mpirun','-np','2',exe,atoms],log=work/'initialization.log',timeout=120)
    strict = link(work,'tests/ion_stages/check_adjacent_ion_fractions.f90')
    run([strict],log=work/'strict_equilibrium.log',timeout=120)
    result = {'passed': True, 'mpi_ranks': 2,
              'checks': ['neutral cache and direct populations', 'neutral cached electron balance',
                         'analytic illuminated populations', 'vacuum opacity convention',
                         'generic disconnected equilibrium still rejected']}
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
