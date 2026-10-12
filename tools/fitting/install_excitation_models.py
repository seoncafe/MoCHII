"""Install coherent approved metal excitation models without regenerating data.

Record of the October 11, 2026 installation of data/atomic/excitation/
mocassin_v2 and extended.  'extended' (the root CHIANTI subsets of HEAD plus
v2 copies) was removed the same day, superseded by the full-CHIANTI 'chianti'
model (tools/fitting/chianti_excitation_model_coverage.py), and the root
metal atoms it copied from were removed with it, so this script documents how
mocassin_v2 was assembled and cannot be rerun on the current tree."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'md/atomic_stage_expansion/validation_2026-10-11'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def numeric_records(text):
    return [s for s in text.splitlines() if s.strip() and not s.lstrip().startswith('#')]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    args=parser.parse_args()
    work=args.work_dir.resolve()
    target=ROOT/'data/atomic/excitation'
    if work.exists() or target.exists():
        raise FileExistsError('Preserve previous model installation')
    atom_manifest=ROOT/'md/atomic_stage_expansion/phase1_77_mocassin_v2/candidate_atoms/manifest.json'
    fit_manifest=BASE/'mocassin_cooling_native/manifest.json'
    native_fits=BASE/'mocassin_cooling_native_exact_endpoints/native/result.json'
    checks=json.loads(native_fits.read_text());assert checks['passed']
    fit_checks={r['ion']:r for r in checks['ions']}
    atoms=json.loads(atom_manifest.read_text())['atoms']
    fits={r['ion']:r for r in json.loads(fit_manifest.read_text())['fits']}
    for r in atoms:assert sha(Path(r['path']))==r['native_sha256']
    for tag,r in fits.items():
        assert sha(Path(r['path']))==r['sha256']==fit_checks[tag]['coefficients_sha256']
    fe=json.loads((BASE/'installed_feii_pair/result.json').read_text())
    assert fe['passed'] and fe['production_installed']
    for r in fe['records']:assert sha(Path(r['production_path']))==r['production_sha256']
    mgi=BASE/'mgi_constant_source/products/nlevel_mg_1.txt'
    mgi_review=json.loads((BASE/'mgi_constant_source/products/result.json').read_text())
    assert mgi_review['passed']
    work.mkdir(parents=True)
    staged=work/'excitation';staged.mkdir()
    extended=staged/'extended';extended.mkdir()
    reference=staged/'mocassin_v2';reference.mkdir()
    # Preserve the complete previously selected models and paired cooling.
    for pattern in ['nlevel_*.txt','cooling_*.txt','hcoll_*.txt']:
        for path in (ROOT/'data/atomic').glob(pattern):
            if path.name=='cooling_fit_parameters.txt' or path.name.startswith(('nlevel_h_','cooling_h_','hcoll_h_')) or '_full' in path.stem:continue
            shutil.copyfile(path,extended/path.name)
    hot={'c_1','s_3','ar_6'}
    for r in atoms:
        tag=f"{r['element']}_{r['stage']}"
        original=Path(r['path']).read_text();text=original
        if tag in hot:
            text=text.replace('# Cold continuation is an approximation; no continuation above the final node.',
                '# Approved constant final-node Upsilon through 100000 K; atomic accuracy remains unverified.')
            text=text.replace('DENSITY_RATIO 1','UPSILON_CONSTANT_UPPER_K 100000\nDENSITY_RATIO 1')
        (reference/f'nlevel_{tag}.txt').write_text(text)
        # Hot extension changes the numerical model. Previously qualified
        # cold-only fits must not be extrapolated as matching cooling.
        if tag in fits and tag not in hot and tag!='fe_2':
            shutil.copyfile(fits[tag]['path'],reference/f'cooling_{tag}.txt')
    for filename in ['nlevel_fe_2.txt','cooling_fe_2.txt']:
        shutil.copyfile(ROOT/'data/atomic'/filename,reference/filename)
    text=mgi.read_text().replace('hot extension pending user acceptance.','hot extension approved October 11, 2026.')
    (reference/'nlevel_mg_1.txt').write_text(text)
    # Mg I has no larger retained atom. Both selections must include its
    # approved source model rather than discard its excitation and radiation.
    (extended/'nlevel_mg_1.txt').write_text(text)
    # Optional hydrogen collisions use the same physical low-lying levels.
    for path in (ROOT/'data/atomic').glob('hcoll_*.txt'):
        if path.name.startswith('hcoll_h_'):continue
        tag=path.stem.removeprefix('hcoll_')
        old=ROOT/'data/atomic'/f'nlevel_{tag}.txt';new=reference/f'nlevel_{tag}.txt'
        if not new.exists():continue
        def levels(p):
            records=[s.split() for s in numeric_records(p.read_text())]
            i=next(i for i,r in enumerate(records) if r[0]=='NLEV')
            return [(float(r[1]),float(r[2])) for r in records[i+1:i+1+int(records[i][1])]]
        old_levels=levels(old);new_levels=levels(new)
        # These files address the ground fine-structure states only.
        addressed={int(r.split()[i]) for r in numeric_records(path.read_text()) if not r.startswith('NHCOL') for i in (0,1)}
        for index in sorted(addressed):
            a,b=old_levels[index-1],new_levels[index-1]
            assert a[0]==b[0] and abs(a[1]-b[1])<0.05,(tag,index,a,b)
        shutil.copyfile(path,reference/path.name)
    manifests={}
    for name, directory in [('mocassin_v2',reference),('extended',extended)]:
        records=[dict(file=p.name,sha256=sha(p)) for p in sorted(directory.glob('*.txt'))]
        manifest=dict(model=name,approved_date='2026-10-11',files=records,
                      atomic_accuracy_verified=False,complete_77_stage_goal=False,
                      missing_reference_atom=['ca_7'],
                      absolute_cooling_without_fit=sorted(p.stem.removeprefix('nlevel_') for p in directory.glob('nlevel_*.txt') if not (directory/p.name.replace('nlevel_','cooling_',1)).exists()))
        (directory/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        manifests[name]=manifest
    # Directory selection is a path change for existing models; primary
    # scientific products and completed diagnostics are preserved.
    shutil.copytree(staged,target)
    for name,manifest in manifests.items():
        for r in manifest['files']:assert sha(target/name/r['file'])==r['sha256']
    result=dict(passed=True,production_installed=True,models=manifests,
                regenerated_scientific_products=False,
                source_sha256={str(p):sha(p) for p in [Path(__file__),atom_manifest,fit_manifest,native_fits,mgi]})
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Installed coherent excitation directories; Ca VII remains physically unresolved.')


if __name__=='__main__':main()
