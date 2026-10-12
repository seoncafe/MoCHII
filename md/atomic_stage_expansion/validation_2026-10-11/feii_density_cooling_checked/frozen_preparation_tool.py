"""Prepare source-consistent cooling for the qualified Fe II radiation model.

Use the verified 142-level source atom with the 40 existing forbidden A values.
No production file or physical-default decision is changed. An exponential
fit is optional; rejection retains the absolute finite-level representation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import re
import shutil

import matplotlib.pyplot as plt
import numpy as np
from mocassin_excitation_profiles import source_atom
from mocassin_cooling_fits import exponential_fit, log_cooling
from mocassin_v2_atom_inventory import MOCASSIN

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'md/atomic_stage_expansion/validation_2026-10-11/feii_forbidden_effect/products'


def radiation_model():
    comparison=json.loads((SOURCE/'result.json').read_text())
    if not comparison['passed'] or comparison['added_radiative_pairs']!=40:
        raise ValueError('Incomplete forbidden-transition comparison')
    for name,digest in comparison['source_sha256'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Compared Fe II source changed')
    manifest=ROOT/'md/atomic_stage_expansion/phase1_77_mocassin_v2/candidate_atoms/manifest.json'
    record=next(r for r in json.loads(manifest.read_text())['atoms'] if r['element']=='fe' and r['stage']==2)
    external=MOCASSIN/record['relative_file']
    if hashlib.sha256(external.read_bytes()).hexdigest()!=record['source_sha256']:
        raise ValueError('Original Fe II collision/radiation source changed')
    levels,radiation,collisions=source_atom(record)
    installed=ROOT/'data/atomic/nlevel_fe_2.txt'
    rows=[s.split() for s in installed.read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
    start=next(i for i,r in enumerate(rows) if r[0]=='NRAD');n=int(rows[start][1])
    for row in rows[start+1:start+1+n]:
        l,u=int(row[0]),int(row[1]);a=float(row[2])
        if a<=0 or levels[u][1]<=levels[l][1]:raise ValueError('Invalid preserved forbidden transition')
        radiation.append((l,u,0.,0.,a))
    return levels,radiation,collisions


def prepare(work):
    if work.exists():raise FileExistsError('Preserve prior candidate products')
    work.mkdir(parents=True)
    levels,radiation,collisions=radiation_model()
    candidate=SOURCE/'nlevel_fe_2_with_forbidden.txt'
    try:
        amplitude,excitation,error,condition,T,power=exponential_fit(levels,radiation,collisions,100000.)
    except ValueError as reason:
        result=dict(passed=True,fit_accepted=False,representation='absolute_finite_level',reason=str(reason))
    else:
        fit=work/'cooling_fe_2.txt'
        lines=['# Qualified Fe II: 142 supplied source levels plus 40 existing forbidden A values.',
            '# This fit uses the same radiation and collision model as the candidate atom.',
            '# Model interval: 1000 <= T <= 100000 K; physical accuracy remains unverified.',
            '# Strict low-density statistical equilibrium, including radiatively closed states.',
            '# Signed terms represent total positive cooling, not individual emission rates.',
            '# Candidate only; no physical-default acceptance or production installation.',
            '# Source atom SHA-256: '+hashlib.sha256(candidate.read_bytes()).hexdigest(),str(len(amplitude))]
        lines += [f'{a:.17g} {e:.17g}' for a,e in zip(amplitude,excitation)]
        fit.write_text('\n'.join(lines)+'\n')
        computed=np.sum(amplitude[:,None]*np.exp(-excitation[:,None]/T),axis=0)/np.sqrt(T)
        np.savez_compressed(work/'low_density_profiles.npz',temperature_K=T,source=power,fitted=computed)
        fig,axes=plt.subplots(2,1,figsize=(7,7),sharex=True,constrained_layout=True)
        axes[0].loglog(T,power,label='Matching source model');axes[0].loglog(T,computed,':',label='Signed exponential fit')
        axes[0].legend();axes[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
        axes[1].semilogx(T,computed/power-1);axes[1].set_xlabel(r'$T$ [K]');axes[1].set_ylabel('Relative representation difference')
        fig.suptitle('Fe II: 142 levels and 40 preserved forbidden transitions')
        fig.savefig(work/'low_density_profiles.pdf');plt.close(fig)
        result=dict(passed=True,fit_accepted=True,representation='matching_low_density_fit',terms=len(amplitude),
            maximum_relative_source_error=error,maximum_cancellation=condition,
            fit_path=str(fit),fit_sha256=hashlib.sha256(fit.read_bytes()).hexdigest())
    result.update(production_modified=False,physical_default_accepted=False,atomic_accuracy_verified=False,
        source_atom_path=str(candidate),source_atom_sha256=hashlib.sha256(candidate.read_bytes()).hexdigest(),
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),SOURCE/'result.json',candidate,ROOT/'data/atomic/nlevel_fe_2.txt']})
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Prepared Fe II matching cooling representation:',result['representation'],flush=True)


def native(work,prepared):
    if work.exists():raise FileExistsError('Preserve previous native products')
    work.mkdir(parents=True)
    metadata=json.loads((prepared/'result.json').read_text())
    if not metadata['passed']:raise ValueError('Candidate preparation incomplete')
    if not metadata['fit_accepted']:
        result=dict(passed=True,native_fit_checked=False,fit_accepted=False,
            reason='No accepted exponential fit; absolute atom remains the calculated representation',production_modified=False)
    else:
        radiation_model()
        fit=Path(metadata['fit_path'])
        if hashlib.sha256(fit.read_bytes()).hexdigest()!=metadata['fit_sha256']:raise ValueError('Prepared fit changed')
        # Compile the actual linked reader, preserving all production source.
        sys.path.insert(0,str(ROOT/'tests/ion_stages'))
        from run_checks import production_objects
        objects=[p for p in production_objects() if p.name!='main.o'];hdf=Path('/data/opt/hdf5_intel/lib');exe=work/'cooling_fit.x'
        driver=ROOT/'tests/ion_stages/check_cooling_fit.f90'
        subprocess.run(['mpiifort','-cpp','-DMPI','-DHDF5','-O0','-no-ipo','-check','all,noarg_temp_created',
            '-fpe0','-traceback','-qopenmp','-Isrc','-module',str(work),'-o',str(exe),str(driver),*map(str,objects),
            '-lcfitsio',str(ROOT/'SEDust/sed/lib/libsedust.a'),str(hdf/'libhdf5_fortran.a'),str(hdf/'libhdf5.a'),
            '-lsz','-ldl','-lz','-lm','-qopenmp'],cwd=ROOT,check=True,timeout=120)
        subprocess.run(['python3','tools/fitting/validate_native_cooling_coefficients.py','--coefficients',str(fit),
            '--driver',str(exe),'--work-dir',str(work/'reader')],cwd=ROOT,check=True,timeout=180)
        rows=np.loadtxt(work/'reader/native_cooling.txt');log_source=log_cooling(*radiation_model(),rows[:,0])
        error=float(max(abs(np.expm1(np.log(rows[:,1])-log_source))))
        if error>=.005:raise ValueError('Actual native fit differs from matching Fe II source')
        result=dict(passed=True,native_fit_checked=True,fit_accepted=True,samples=801,
            maximum_relative_source_error=error,production_modified=False,
            source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                           [Path(__file__),prepared/'result.json',fit,driver,exe,work/'reader/result.json']})
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Completed Fe II native representation check; fit checked:',result['native_fit_checked'],flush=True)


def feii_density_cooling(work,prepared):
    if work.exists():raise FileExistsError('Preserve previous finite-density calculation')
    work.mkdir(parents=True)
    metadata=json.loads((prepared/'result.json').read_text())
    if not metadata['passed'] or not metadata['fit_accepted']:
        raise ValueError('Finite-density fit check requires an accepted matching fit')
    radiation_model()
    atom=Path(metadata['source_atom_path']);fit=Path(metadata['fit_path'])
    if hashlib.sha256(atom.read_bytes()).hexdigest()!=metadata['source_atom_sha256']:
        raise ValueError('Prepared radiation model changed')
    if hashlib.sha256(fit.read_bytes()).hexdigest()!=metadata['fit_sha256']:
        raise ValueError('Prepared cooling fit changed')
    atomic=work/'atomic';atomic.mkdir();shutil.copyfile(atom,atomic/'nlevel_fe_2.txt')
    atom_source=ROOT/'src/nlevel_mod.f90';module_source=ROOT/'src/nlevel_cooling_mod.f90'
    private_atom=work/'nlevel_mod.f90';private_atom.write_bytes(atom_source.read_bytes())
    text=module_source.read_text();marker='    limited_grid = reg(ie)%source_ceiling_ln_K < log(10.0_wp)*gT(NT)'
    if text.count(marker)!=1:raise ValueError('Cooling interval declaration changed')
    text=text.replace(marker,'    ! Diagnostic only: build the unchanged grid prefix through the fitted 100000 K bound.\n    reg(ie)%source_ceiling_ln_K = min(reg(ie)%source_ceiling_ln_K,log(100000.0_wp))\n'+marker)
    public='  public :: nlevel_cooling_lam_se, nlevel_cooling_report\n'
    text=text.replace(public,public+'  public :: reg\n')
    private_module=work/'nlevel_cooling_mod.f90';private_module.write_text(text)
    flags=['mpiifort','-cpp','-DMPI','-DHDF5','-O2','-no-ipo','-fp-model','precise',
           '-fpe0','-traceback','-qopenmp',
           '-I'+str(work),'-Isrc','-module',str(work)]
    for source in (private_atom,private_module):
        subprocess.run([*flags,'-c',str(source),'-o',str(source.with_suffix('.o'))],cwd=ROOT,check=True,timeout=120)
    original_driver=ROOT/'tests/ion_stages/check_registered_excitation_atom.f90';text=original_driver.read_text()
    text=text.replace('nlevel_cooling_apply, nlevel_cooling_report','nlevel_cooling_apply, nlevel_cooling_report, reg')
    text=text.replace('  implicit none','  use cooling_mod, only: cooling_fit_type, cooling_load, cooling_eval\n  implicit none\n  type(cooling_fit_type) :: matching_fit\n  character(len=1024) :: fit_file',1)
    if text.count('command_argument_count() /= 4') != 1:
        raise ValueError('Diagnostic argument declaration changed')
    text=text.replace('command_argument_count() /= 4','command_argument_count() /= 5')
    text=text.replace('ATOMIC_DIR ELEMENT STAGE OUTPUT','ATOMIC_DIR ELEMENT STAGE OUTPUT MATCHING_FIT')
    text=text.replace('  call get_command_argument(4, output)','  call get_command_argument(4, output)\n  call get_command_argument(5, fit_file)\n  call cooling_load(trim(fit_file),matching_fit)')
    start=text.index('  source_ceiling_ln_K = log(1.0e6_wp)');end=text.index('  open(newunit=unit',start)
    text=text[:start]+'  source_ceiling_ln_K = log(1.0e5_wp)\n'+text[end:]
    text=text.replace('grid_t(5) = [100.0_wp, 1000.0_wp, &\n       1.0e4_wp, 1.0e5_wp, 1.0e6_wp]',
        'grid_t(5) = [1000.0_wp, 10.0_wp**3.5_wp, &\n       1.0e4_wp, 10.0_wp**4.5_wp, 1.0e5_wp]')
    text=text.replace('middle_log_t(5) = [2.005_wp, 3.005_wp, &\n       4.005_wp, 5.005_wp, 5.995_wp]',
        'middle_log_t(5) = [3.005_wp, 3.505_wp, &\n       4.005_wp, 4.505_wp, 4.995_wp]')
    text=text.replace('reference = nlevel_cooling_at(atom, t, grid_ne(1))','reference = cooling_eval(matching_fit,t)')
    text=text.replace('pass_index == 1 .and. apply_error > 1.0e-10_wp','pass_index == 1 .and. apply_error > 5.0e-3_wp')
    text=text.replace("'# source_ceiling_ln_K '","'# model_ceiling_ln_K '")
    marker="  write(unit, '(a,i0)') '# nlev ', atom%nlev"
    text=text.replace(marker,marker+"\n  write(unit, '(a,i0)') '# temperature_nodes ', size(reg(idx)%logT)\n  write(unit, '(a,i0)') '# density_nodes ', size(reg(idx)%lambda,2)")
    driver=work/'feii_density_cooling.f90';driver.write_text(text);exe=work/'feii_density_cooling.x'
    sys.path.insert(0,str(ROOT/'tests/ion_stages'));from run_checks import production_objects
    objects=[p for p in production_objects() if p.name not in ('main.o','nlevel_mod.o','nlevel_cooling_mod.o')]
    hdf=Path('/data/opt/hdf5_intel/lib')
    subprocess.run([*flags,'-o',str(exe),str(driver),str(private_atom.with_suffix('.o')),
        str(private_module.with_suffix('.o')),*map(str,objects),'-lcfitsio',
        str(ROOT/'SEDust/sed/lib/libsedust.a'),str(hdf/'libhdf5_fortran.a'),str(hdf/'libhdf5.a'),
        '-lsz','-ldl','-lz','-lm','-qopenmp'],cwd=ROOT,check=True,timeout=120)
    output=work/'registered.txt'
    with (work/'native.log').open('w') as stream:
        subprocess.run([str(exe),str(atomic),'fe','2',str(output),str(fit)],cwd=ROOT,
                       stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=1200)
    text=output.read_text()
    if '# PASS' not in text or '# failures 0' not in text:raise ValueError('Finite-density checks incomplete')
    rows=np.array([[float(re.sub(r'(?<=[0-9.])([+-]\d{3})$',r'E\1',v)) for v in line.split()[1:]]
                   for line in text.splitlines() if not line.startswith('#')])
    if rows.shape!=(40,9) or not np.isfinite(rows).all():raise ValueError('Incomplete native temperature/density profiles')
    error=float(max(abs(rows[:,4]/rows[:,2]-1)))
    if error>.005:raise ValueError('Matching thermal/direct line energy discrepancy exceeds criterion')
    summary={s.split()[1]:float(s.split()[2]) for s in text.splitlines() if s.startswith('# ') and len(s.split())==3
             and s.split()[1] in ('temperature_nodes','density_nodes','worst_middle_relative_error')}
    if summary['temperature_nodes']!=301 or summary['density_nodes']!=161:raise ValueError('Unmeasured or unexpected table shape')
    result=dict(passed=True,samples=40,production_modified=False,physical_default_accepted=False,
        temperature_K=[1000,100000],table_cells=301*161,maximum_relative_thermal_direct_energy_difference=error,
        maximum_relative_table_interpolation_difference=summary['worst_middle_relative_error'],
        compilation='Isolated O2 precise diagnostics; original source atom equations, grid prefix limited to fit interval',
        scope='Actual matching fit read and nlevel_cooling_apply versus direct source atom cooling; no species thermal balance, line writer, or spatial model claim',
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__),atom,fit,atom_source,module_source,private_atom,private_module,driver,original_driver,exe,output]})
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: 40 Fe II matching finite-density cooling checks',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('step',choices=('prepare','native','density'))
    p.add_argument('--work-dir',type=Path,required=True);p.add_argument('--prepared',type=Path)
    a=p.parse_args()
    if a.step=='prepare':prepare(a.work_dir.resolve())
    elif a.prepared and a.step=='native':native(a.work_dir.resolve(),a.prepared.resolve())
    elif a.prepared:feii_density_cooling(a.work_dir.resolve(),a.prepared.resolve())
    else:raise ValueError('Native check needs candidate preparation')


if __name__=='__main__':main()
