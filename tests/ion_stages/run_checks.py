#!/usr/bin/env python3
"""Run the ion-stage checks against the built production objects.

Every Fortran check is linked to the objects the Makefile lists (all but
main.o) with the flags those objects were compiled with, read from
src/.build_configuration, plus -qmkl=sequential for the LAPACK solve of the
level populations (nlevel_mod::level_populations).  A check that needs a
private module entity compiles a copy of that module into its own work
directory and links the copy in place of the production object; the
production sources and objects are never changed.  Production atomic files
are never written: executables, fixtures and outputs live in --work-dir.

Steps:
  build     make MoCHII.x, then require every object to be newer than its
            source
  checks    the ion-stage solvers, the stage registry, rates, inputs and the
            thermal and charge-exchange consumers
  extended  every excitation atom of both models, and the atom generators
  full      checks, then extended
A failing check is recorded and the remaining checks still run; the step
writes <work-dir>/<step>_result.json and exits 1 if any check failed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT/'tests/ion_stages'
HDF5_PREFIX = Path(os.environ.get('HDF5_PREFIX', '/data/opt/hdf5_intel'))
SEDUST = ROOT/os.environ.get('SEDUST_LIBDIR', 'SEDust/sed/lib')
LIBRARIES = ['-lcfitsio', '-L/usr/local/lib', str(SEDUST/'libsedust.a'),
             str(HDF5_PREFIX/'lib/libhdf5_fortran.a'), str(HDF5_PREFIX/'lib/libhdf5.a'),
             '-lsz', '-ldl', '-lz', '-lm', '-qopenmp', '-qmkl=sequential']
SINGLE_THREAD = dict(os.environ, OMP_NUM_THREADS='1')


def command_text(command):
    return shlex.join([str(v) for v in command])


def run(command, cwd=ROOT, log=None, timeout=1800, expect_failure=None):
    """Run a command; with expect_failure, require a nonzero exit whose
    output contains that text."""
    print(command_text(command), flush=True)
    result = subprocess.run([str(v) for v in command], cwd=cwd, env=SINGLE_THREAD,
                            capture_output=True, text=True, timeout=timeout)
    output = result.stdout+result.stderr
    if log is not None:
        Path(log).write_text(command_text(command)+'\n'+output)
    if expect_failure is None:
        if result.returncode != 0:
            raise RuntimeError(f'exit {result.returncode}: {command_text(command)}\n{output[-4000:]}')
    elif result.returncode == 0 or expect_failure not in output:
        raise RuntimeError(f'expected a failure containing {expect_failure!r}: '
                           f'{command_text(command)}\n{output[-4000:]}')
    return output


def production_objects():
    """The Makefile's object list (without main.o), each newer than its source."""
    text = (ROOT/'Makefile').read_text()
    block = re.search(r'^OBJS\s*=\s*\\\n((?:.*\\\n)*.*)', text, re.M).group(1)
    objects = [ROOT/'src'/name for name in re.findall(r'(\S+\.o)', block) if name != 'main.o']
    stale = [o.name for o in objects
             if not o.exists() or o.stat().st_mtime < o.with_suffix('.f90').stat().st_mtime]
    if stale:
        raise RuntimeError('missing or stale production objects (rebuild with make): '+' '.join(stale))
    return objects


def fortran_command():
    """Compiler and flags the production objects were built with."""
    lines = (ROOT/'src/.build_configuration').read_text().splitlines()
    return [lines[0], *shlex.split(lines[1]), lines[3]]


def compile_private_module(source, work, replacements):
    """Compile a copy of a production module with textual replacements (each
    must occur exactly once) into work; return the object."""
    text = source.read_text()
    for old, new in replacements:
        if text.count(old) != 1:
            raise RuntimeError(f'{source.name}: expected one occurrence of {old!r}')
        text = text.replace(old, new)
    copy = work/source.name
    copy.write_text(text)
    obj = copy.with_suffix('.o')
    run([*fortran_command(), '-I'+str(work), '-Isrc', '-module', str(work),
         '-c', copy, '-o', obj], log=work/(source.stem+'_compile.log'))
    return obj


def link(work, source, private_objects=(), name=None):
    """Link a check driver to the production objects; a private object
    replaces the production object of the same name."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    replaced = {p.name for p in private_objects}
    objects = [o for o in production_objects() if o.name not in replaced]
    exe = work/((name or Path(source).stem)+'.x')
    run([*fortran_command(), '-I'+str(work), '-Isrc', '-module', str(work), '-o', exe,
         ROOT/source, *private_objects, *objects, *LIBRARIES],
        log=work/(exe.stem+'_link.log'))
    return exe


def fresh(work):
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    return work


def python_check(work, script, *arguments):
    output = run([sys.executable, HERE/script, '--work-dir', work, *arguments],
                 log=work/'run.log', timeout=7200)
    result = work/'result.json'
    if not result.exists():
        raise RuntimeError(script+' wrote no result.json')
    summary = json.loads(result.read_text())
    if not summary.get('passed'):
        raise RuntimeError(script+' did not pass: '+output[-2000:])
    return summary


# ---------------------------------------------------------------------------
# Checks.  Each takes its own fresh work directory and returns a summary.
CHECKS = []


def check(*steps):
    def register(function):
        CHECKS.append((steps, function))
        return function
    return register


@check('checks')
def adjacent_ion_fractions(work):
    exe = link(work, 'tests/ion_stages/check_adjacent_ion_fractions.f90')
    return run([exe], log=work/'run.log').strip()


@check('checks')
def nonadjacent_ion_equilibrium(work):
    # ion_stage_stationary_mod is not a production object yet (the general
    # stationary solve the multi-stage Auger matrix will need); compile it
    # into the work directory with the production flags.
    source = ROOT/'src/ion_stage_stationary_mod.f90'
    obj = work/'ion_stage_stationary_mod.o'
    run([*fortran_command(), '-I'+str(work), '-Isrc', '-module', str(work), '-c', source,
         '-o', obj], log=work/'ion_stage_stationary_mod_compile.log')
    exe = link(work, 'tests/ion_stages/check_nonadjacent_ion_equilibrium.f90', [obj])
    return run([exe], log=work/'run.log').strip()


SYNTHETIC_FE_TRANSITION = ['ETH {e}', 'PHOTO {e} 10 1 1 2 1 0 1', 'PHINNER 1000',
                           'PHSHELL 0 {e} 10 1 1 2 1', 'CI 1 0 0 1 0',
                           'RR 1e-11 0.5 1e4 1e4 0 1', 'DR 0', 'CXI 0', 'CXR 0',
                           'CXHEAT 0 0 0']


def synthetic_atomic(directory, element_text):
    """An atomic directory with one element file and an empty default
    excitation model (nlevel_mod::metal_excitation_data_file requires the
    model directory's manifest.json and reads nothing else from it)."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory/'element_fe.txt').write_text(element_text)
    model = directory/'excitation/chianti'
    model.mkdir(parents=True, exist_ok=True)
    (model/'manifest.json').write_text(json.dumps(
        {'model': 'chianti', 'files': [], 'scope': 'synthetic fixture without atoms'})+'\n')
    return directory


@check('checks')
def stage_registry(work):
    rows = ['# Synthetic adjacent-rate fixture; not Fe atomic data.', 'ELEMENT fe 26 7']
    for it in range(1, 7):
        rows.append(f'TRANSITION {it}')
        rows.extend(r.format(e=10*it) for r in SYNTHETIC_FE_TRANSITION)
    valid = '\n'.join(rows)+'\n'
    fixture = synthetic_atomic(work/'synthetic_atomic', valid)
    exe = link(work, 'tests/ion_stages/check_stage_registry.f90')
    summary = [run([exe, fixture, work/'seven_stage.fits'], log=work/'run.log').strip()]
    run([exe, fixture, work/'rejected_restart.fits', 'reject'], log=work/'rejected_restart.log',
        expect_failure='requires gas_niter >= 1')
    summary.append('incomplete fixed-state restart rejected')
    malformed = {
        'stage_bounds': (valid.replace('ELEMENT fe 26 7', 'ELEMENT fe 26 8'), 'stage count out of range'),
        'transition_bounds': (valid.replace('TRANSITION 1', 'TRANSITION 0', 1), 'TRANSITION index out of range'),
        'dr_bounds': (valid.replace('DR 0', 'DR 10', 1), 'DR term count out of range'),
        'negative_dr': (valid.replace('DR 0', 'DR -1', 1), 'DR term count out of range'),
        'before_transition': (valid.replace('TRANSITION 1\n', '', 1), 'without a valid TRANSITION'),
        'missing_ci': (valid.replace('CI 1 0 0 1 0\n', '', 1), 'incomplete atomic transition records'),
        'duplicate_dr': (valid.replace('DR 0', 'DR 0\nDR 0', 1), 'duplicate atomic rate record'),
        'nonfinite_rr': (valid.replace('RR 1e-11', 'RR NaN', 1), 'nonfinite atomic rate coefficients'),
        'negative_rr': (valid.replace('RR 1e-11', 'RR -1e-11', 1), 'nonpositive RR coefficient'),
        'identity': (valid.replace('ELEMENT fe 26 7', 'ELEMENT fe 25 7'), 'atomic number disagrees'),
        'unknown_record': (valid+'UNKNOWN 1\n', 'unknown atomic record'),
        'order': (valid.replace('TRANSITION 1', 'TRANSITION 2', 1), 'TRANSITION records must be consecutive'),
    }
    for case, (contents, diagnostic) in malformed.items():
        bad = synthetic_atomic(work/('malformed_'+case), contents)
        run([exe, bad, bad/'unused.fits'], log=bad/'rejection.log', expect_failure=diagnostic)
    summary.append(f'{len(malformed)} malformed element inputs rejected')
    return summary


@check('checks')
def physical_thresholds(work):
    return python_check(work, 'check_physical_thresholds.py')


@check('checks')
def metal_photoabsorption(work):
    return python_check(work, 'check_metal_photoabsorption.py')


@check('checks')
def charge_transfer_tables(work):
    return python_check(work, 'check_charge_transfer_tables.py')


@check('checks')
def langevin_charge_exchange(work):
    return python_check(work, 'check_langevin_charge_exchange.py')


@check('checks')
def sulfur_charge_exchange_model(work):
    return python_check(work, 'check_sulfur_charge_exchange_model.py')


@check('checks')
def fe_recombination(work):
    return python_check(work, 'check_fe_recombination.py')


@check('checks')
def optional_metal_physics(work):
    return python_check(work, 'check_optional_metal_physics.py')


@check('checks')
def optional_thermal_terms(work):
    return python_check(work, 'check_optional_thermal_terms.py')


@check('checks')
def neutral_hydrogen_excitation(work):
    return python_check(work, 'check_neutral_hydrogen_excitation.py')


@check('checks')
def hei_metastable_input(work):
    return python_check(work, 'check_hei_metastable_input.py')


@check('checks')
def excitation_model_selection(work):
    return python_check(work, 'check_excitation_model_selection.py')


@check('checks')
def cooling_temperature_interpolation(work):
    return python_check(work, 'check_cooling_temperature_interpolation.py')


@check('checks')
def gaunt_stages(work):
    exe = link(work, 'tests/ion_stages/check_gaunt_stages.f90')
    return run([exe, work/'gaunt_nodes.txt'], log=work/'run.log').strip().splitlines()[-1]


@check('checks')
def ci_stage_factors(work):
    return run([sys.executable, HERE/'check_ci_stage_factors.py', '--work-dir', work],
               log=work/'run.log').strip().splitlines()[-1]


# The charge-exchange gates of tests/charge_exchange, run with the same link.
for _name in ['check_charge_exchange', 'check_cell_convergence', 'check_detailed_balance',
              'check_helium_charge_exchange', 'check_charge_exchange_heating']:
    def _charge_exchange_gate(work, _name=_name):
        exe = link(work, 'tests/charge_exchange/'+_name+'.f90')
        return run([exe], log=work/'run.log').strip().splitlines()[-1]
    _charge_exchange_gate.__name__ = 'charge_exchange_'+_name.removeprefix('check_')
    check('checks')(_charge_exchange_gate)


@check('checks')
def selective_generation(work):
    return run([sys.executable, 'tests/fitting/test_selective_generation.py'],
               log=work/'run.log').strip().splitlines()[-1:]


@check('checks')
def python_rates(work):
    return run([sys.executable, 'tests/rates/check_python_rates.py'],
               log=work/'run.log').strip().splitlines()[-1:]


@check('extended')
def excitation_atoms(work):
    return python_check(work, 'check_excitation_atoms.py')


# The CHIANTI atom generator (tools/fitting/write_excitation_atom.py) on the
# ions whose radiative data it corrects; Mg III is compared byte for byte
# with the installed default-model atom.
@check('extended')
def mgiii_electric_dipole(work):
    return run([sys.executable, HERE/'check_mgiii_electric_dipole.py',
                ROOT/'data/atomic/excitation/chianti/nlevel_mg_3.txt'],
               log=work/'run.log').strip().splitlines()[-1]


@check('extended')
def nei_spontaneous_probabilities(work):
    return run([sys.executable, HERE/'check_nei_spontaneous_probabilities.py'],
               log=work/'run.log').strip().splitlines()[-1]


@check('extended')
def oxygen_electric_dipole(work):
    return run([sys.executable, HERE/'check_oxygen_electric_dipole.py'],
               log=work/'run.log').strip().splitlines()[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('step', choices=['build', 'checks', 'extended', 'full'])
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--only', nargs='*', help='run only the named checks')
    args = parser.parse_args()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    result = {'step': args.step, 'passed': False}
    if args.step == 'build':
        # A shared define change requires recompiling all consumers.
        run(['make', '-j2', 'MoCHII.x'], timeout=7200, log=work/'build.log')
        executable = ROOT/'MoCHII.x'
        if not executable.is_file() or executable.stat().st_size == 0:
            raise RuntimeError('Build produced no executable')
        production_objects()
        result['executable_sha256'] = hashlib.sha256(executable.read_bytes()).hexdigest()
        result['passed'] = True
    else:
        production_objects()
        steps = {'checks': ['checks'], 'extended': ['extended'],
                 'full': ['checks', 'extended']}[args.step]
        result['flags'] = fortran_command()
        result['checks'] = {}
        for step in steps:
            for check_steps, function in CHECKS:
                if step not in check_steps:
                    continue
                if args.only and function.__name__ not in args.only:
                    continue
                start = time.time()
                try:
                    summary = function(fresh(work/function.__name__))
                    entry = {'passed': True, 'summary': summary}
                    print('PASS', function.__name__, flush=True)
                except Exception as error:
                    entry = {'passed': False, 'error': str(error)[-6000:],
                             'traceback': traceback.format_exc()[-6000:]}
                    print('FAIL', function.__name__, str(error)[-2000:], flush=True)
                entry['seconds'] = round(time.time()-start, 1)
                result['checks'][function.__name__] = entry
        result['passed'] = all(e['passed'] for e in result['checks'].values())
        result['failed'] = [k for k, e in result['checks'].items() if not e['passed']]
    (work/(args.step+'_result.json')).write_text(json.dumps(result, indent=2, default=str)+'\n')
    print(('PASS' if result['passed'] else 'FAIL')+': '+args.step+
          ('' if result['passed'] else ' -- failed: '+' '.join(result['failed'])), flush=True)
    if not result['passed']:
        sys.exit(1)


if __name__ == '__main__':
    main()
