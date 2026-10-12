"""Measure the actual registered cooling table and validate its consumers.

The compiled driver registers the production table and limits exact-grid and
interior comparisons to the collision-table ceiling, at most 1,000,000 K.
Finite collision tables can shorten the temperature grid. Source-data accuracy
outside the independently accepted 1,000--100,000 K interval is not established
by a table interpolation check.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import shutil

import numpy as np
from atomic_output import atomic_write_text


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir', type=Path, required=True)
    p.add_argument('--element', required=True)
    p.add_argument('--stage', type=int, required=True)
    p.add_argument('--model', choices=['chianti', 'mocassin_v2'], default='chianti',
                   help='Excitation model directory under <source-dir>/excitation')
    p.add_argument('--driver', type=Path, required=True)
    p.add_argument('--work-dir', type=Path, required=True)
    p.add_argument('--timeout', type=int, default=3600)
    p.add_argument('--validate-existing', action='store_true',
                   help='Validate complete saved native outputs without rebuilding a table')
    a = p.parse_args()
    a.element = a.element.lower()
    directory = a.work_dir.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    output = directory/'registered.txt'
    if output.exists() and not a.validate_existing:
        raise FileExistsError('Preserve prior products; use a new validation directory')
    source = a.source_dir.resolve()
    filename = f'nlevel_{a.element}_{a.stage}.txt'
    selected = source/'excitation'/a.model/filename
    atom = selected if selected.exists() else source/filename
    atomic_input = source
    if not selected.exists() and not a.validate_existing:
        # A single candidate atom is placed in an explicit diagnostic model.
        # Preserve primary products; copying changes only their location.
        atomic_input = directory/'atomic_input'
        model = atomic_input/'excitation'/a.model
        model.mkdir(parents=True)
        for quantity in ('nlevel', 'hcoll', 'cooling'):
            original = source/f'{quantity}_{a.element}_{a.stage}.txt'
            if original.exists():
                shutil.copyfile(original, model/original.name)
        (model/'manifest.json').write_text(json.dumps(dict(
            model=a.model, scope='Registered single-atom diagnostic'))+'\n')
    command = ['/usr/bin/time', '-v', '-o', str(directory/'resources.txt'),
               str(a.driver.resolve()), str(atomic_input), a.model, a.element,
               str(a.stage), str(output)]
    if a.validate_existing:
        resources = (directory/'resources.txt').read_text()
        if '\tExit status: 0' not in resources:
            raise ValueError('Saved native calculation did not finish successfully')
    else:
        with (directory/'native.log').open('w') as stream:
            subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                           check=True, timeout=a.timeout)
    text = output.read_text()
    if '# PASS' not in text or '# failures 0' not in text:
        raise ValueError('Registered cooling checks did not pass')
    records = [line.split() for line in text.splitlines() if not line.startswith('#')]
    if len(records) != 40 or any(len(row) != 10 for row in records):
        raise ValueError('Incomplete registered profiles')
    # Fortran ES formats without an explicit exponent width omit E when the
    # exponent has three digits. Read those valid records without rerunning
    # a completed numerical calculation just to change its printed notation.
    values = np.array([[float(re.sub(r'(?<=[0-9.])([+-]\d{3})$', r'E\1', v))
                        for v in row[1:]] for row in records])
    if not np.all(np.isfinite(values)) or np.any(values[:, 2:7] < 0):
        raise ValueError('Invalid registered cooling output')
    summary = {}
    for line in text.splitlines():
        tokens = line.split()
        if len(tokens) == 3 and tokens[0] == '#' and tokens[1] in {
            'nlev', 'nrad', 'nups', 'registration_wall_seconds', 'source_ceiling_ln_K',
            'temperature_nodes', 'density_nodes', 'model_ceiling_ln_K',
            'worst_grid_relative_error', 'worst_middle_relative_error',
            'worst_apply_relative_error'}:
            summary[tokens[1]] = float(tokens[2])
    for line in (directory/'resources.txt').read_text().splitlines():
        if 'Maximum resident set size (kbytes)' in line:
            summary['native_peak_rss_KiB'] = int(line.rsplit(':', 1)[1])
    measured_shape = all(k in summary for k in ('temperature_nodes', 'density_nodes'))
    cells = int(summary['temperature_nodes']*summary['density_nodes']) if measured_shape else None
    result = dict(passed=True, element=a.element, stage=a.stage, samples=40,
                  table_cells=cells, table_cells_measured=measured_shape,
                  tested_temperature_minimum_K=float(values[:,0].min()),
                  tested_temperature_maximum_K=float(values[:,0].max()),
                  source_sha256=hashlib.sha256(atom.read_bytes()).hexdigest(),
                  executable_sha256=hashlib.sha256(a.driver.read_bytes()).hexdigest(),
                  scope='Cooling registration and interpolation; no spatial plasma model', **summary)
    atomic_write_text(directory/'result.json', json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
