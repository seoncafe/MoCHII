"""Check the excitation-model selection and the refusal of removed inputs.

par%metal_excitation_model selects data/atomic/excitation/<model>/ (default
'chianti'); nlevel_mod::metal_excitation_data_file names the files there.
Through the namelist reader and nlevel_load (check_excitation_model_
selection.f90): the default and 'mocassin_v2' load exactly the atoms their
manifest.json lists, every manifest file has its recorded sha256, an
unknown model ('extended', removed 2026-10-11) and a missing model directory
are refused, and the items removed or renamed on 2026-10-11 are refused with
the name of what replaced them.
"""
import argparse
import hashlib
import json
from pathlib import Path

from run_checks import ROOT, link, run

REMOVED = {
    'metal_inner_shell_absorption': '.false.',
    'additional_line_cooling': '.true.',
    'high_charge_exchange': '.true.',
    'fe_levels_full': '.true.',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    work = parser.parse_args().work_dir.resolve()
    exe = link(work, 'tests/ion_stages/check_excitation_model_selection.f90')
    summaries = []
    for option, setting in (('default', ''), ('mocassin_v2', " par%metal_excitation_model='mocassin_v2'\n")):
        inp = work/(option+'.in')
        inp.write_text('&parameters\n'+setting+'/\n')
        output = work/(option+'.txt')
        run([exe, inp, ROOT/'data/atomic', output], log=work/(option+'.log'))
        rows = output.read_text().splitlines()
        model = rows[0].strip()
        if model != ('chianti' if option == 'default' else option):
            raise RuntimeError(f'{option}: selected {model}')
        directory = ROOT/'data/atomic/excitation'/model
        manifest = json.loads((directory/'manifest.json').read_text())
        for record in manifest['files']:
            digest = hashlib.sha256((directory/record['file']).read_bytes()).hexdigest()
            if digest != record['sha256']:
                raise RuntimeError(f'{model}/{record["file"]}: sha256 differs from the manifest')
        loaded = {'_'.join(r.split()[:2]) for r in rows[1:]}
        listed = {r['file'].removeprefix('nlevel_').removesuffix('.txt')
                  for r in manifest['files'] if r['file'].startswith('nlevel_')}
        if loaded != listed:
            raise RuntimeError(f'{model}: loaded {sorted(loaded ^ listed)} differently from the manifest')
        if 'stages_with_atom' in manifest and set(manifest['stages_with_atom']) != loaded:
            raise RuntimeError(f'{model}: stages_with_atom differs from the files')
        summaries.append(dict(option=option, model=model, atoms=len(loaded),
                              levels=sum(int(r.split()[2]) for r in rows[1:]),
                              manifest_files=len(manifest['files'])))
    for option, setting, root, expected in (
            ('extended', " par%metal_excitation_model='extended'\n", ROOT/'data/atomic',
             'must be chianti or mocassin_v2'),
            ('missing', '', work/'missing', 'model directory is missing')):
        inp = work/(option+'.in')
        inp.write_text('&parameters\n'+setting+'/\n')
        run([exe, inp, root, work/(option+'.txt')], log=work/(option+'.log'), expect_failure=expected)
    for name, value in REMOVED.items():
        inp = work/(name+'.in')
        inp.write_text(f'&parameters\n par%{name}={value}\n/\n')
        run([exe, inp, ROOT/'data/atomic', work/(name+'.txt')], log=work/(name+'.log'),
            expect_failure=f'par%{name} (line 2 of')
    result = dict(passed=True, models=summaries, refused_models=2, refused_removed_inputs=list(REMOVED))
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
