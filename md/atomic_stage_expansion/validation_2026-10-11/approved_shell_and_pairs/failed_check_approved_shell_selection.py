"""Check the approved shell input, installed atom pairs, and actual packet consumers."""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
import numpy as np
from run_checks import production_objects
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/fitting'))
from check_outer_shell_photoionization import rates,cross_section,ELEMENTS
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work-dir',type=Path,required=True);a=p.parse_args()
w=a.work_dir.resolve()
if w.exists():raise FileExistsError('Preserve prior approved-selection check')
w.mkdir(parents=True)
flags=['mpiifort','-cpp','-DMPI','-DHDF5','-O0','-no-ipo','-check','all,noarg_temp_created','-fpe0','-traceback','-qopenmp','-I'+str(w),'-Isrc','-module',str(w)]
objects=[x for x in production_objects() if x.name!='main.o'];hdf=Path('/data/opt/hdf5_intel/lib')
libs=['-lcfitsio',str(ROOT/'SEDust/sed/lib/libsedust.a'),str(hdf/'libhdf5_fortran.a'),str(hdf/'libhdf5.a'),'-lsz','-ldl','-lz','-lm','-qopenmp']
source=ROOT/'src/species_mod.f90';private=w/'species_mod.f90'
private.write_text(source.read_text().replace('  private\n','  private\n  public :: outer_shell_photoionization, flat_ie, flat_it, stage_frac, metal_rate_path\n',1))
subprocess.run([*flags,'-c',str(private),'-o',str(w/'species_mod.o')],cwd=ROOT,check=True,timeout=120)
driver=ROOT/'tests/ion_stages/check_outer_shell_photoionization.f90';exe=w/'shell_selection.x'
subprocess.run([*flags,'-o',str(exe),str(driver),str(w/'species_mod.o'),*[str(o) for o in objects if o.name!='species_mod.o'],*libs],cwd=ROOT,check=True,timeout=120)
# Fixtures preserve rate files exactly and bypass unrelated line cooling.
checks=[];guard_checks=0
for roster,location,nt in [('production',ROOT/'data/atomic',32),('staged77',ROOT/'md/atomic_stage_expansion/phase1_77_mocassin_v2/staged_rates',66)]:
 atomic=w/roster;atomic.mkdir()
 for path in location.glob('element_*.txt'):
  shutil.copyfile(path,atomic/path.name)
  count=int(next(s.split()[3] for s in path.read_text().splitlines() if s.startswith('ELEMENT ')))
  el=path.stem.removeprefix('element_')
  for i in range(1,count+1):(atomic/f'cooling_{el}_{i}.txt').write_text('# Zero cooling in cross-section diagnostic only.\n1\n0 0\n')
 shutil.copyfile(ROOT/'data/atomic/nlevel_s_2.txt',atomic/'nlevel_s_2.txt')
 reference=None
 for option,settings,mode,n in [('default','','',18),('off',' par%metal_inner_shell_absorption=.false.\n','',18),('on_below',' par%metal_inner_shell_absorption=.true.\n','below',3)]:
  inp=w/(roster+'_'+option+'.in');inp.write_text('&parameters\n'+settings+'/\n');output=w/(roster+'_'+option+'.txt')
  run=subprocess.run([str(exe),str(atomic),str(output),str(inp),mode],cwd=ROOT,capture_output=True,text=True,timeout=300,env=dict(os.environ,OMP_NUM_THREADS='1'))
  (w/(roster+'_'+option+'.log')).write_text(run.stdout+run.stderr)
  if run.returncode:raise RuntimeError(roster+' '+option+' failed')
  rows=np.loadtxt(output);assert rows.shape==(nt*n,7) and np.isfinite(rows).all()
  if option=='default':reference=rows
  if option=='off':assert np.array_equal(reference,rows)
  data={i:rates(atomic/f'element_{el}.txt') for i,el in enumerate(ELEMENTS,1)}
  worst=0.
  for ie,it,E,sig,heat,selected,selected_heat in rows:
   expected=np.array(cross_section(data[int(ie)][int(it)],E));got=np.array([sig,heat]);positive=expected>0
   assert np.all(got[~positive]==0)
   if positive.any():worst=max(worst,float(max(abs(got[positive]/expected[positive]-1))))
   assert sig==selected and heat==selected_heat
  assert worst<1e-12
  checks.append(dict(roster=roster,option=option,samples=len(rows),maximum_formula_relative_error=worst))
 inp=w/(roster+'_on_above.in');inp.write_text('&parameters\n par%metal_inner_shell_absorption=.true.\n/\n')
 run=subprocess.run([str(exe),str(atomic),str(w/(roster+'_on_above.txt')),str(inp)],cwd=ROOT,capture_output=True,text=True,timeout=300,env=dict(os.environ,OMP_NUM_THREADS='1'))
 (w/(roster+'_on_above.log')).write_text(run.stdout+run.stderr)
 assert run.returncode!=0 and 'unqualified vacancy charge and electron-energy model' in run.stdout+run.stderr
 guard_checks+=1
# Actual installed file readers, without rebuilding unrelated density tables.
installation=ROOT/'md/atomic_stage_expansion/validation_2026-10-11/installed_shared_atom_pairs/result.json'
installed=json.loads(installation.read_text());assert installed['passed'] and installed['production_installed']
pairdriver=ROOT/'tests/ion_stages/check_approved_atomic_data.f90';pairexe=w/'installed_pairs.x'
subprocess.run([*flags,'-o',str(pairexe),str(pairdriver),*map(str,objects),*libs],cwd=ROOT,check=True,timeout=120)
pairs=[]
for row in installed['atoms']:
 atom=Path(row['production_path']);fit=Path(row['selected_cooling_path'])
 assert hashlib.sha256(atom.read_bytes()).hexdigest()==row['native_sha256']
 assert hashlib.sha256(fit.read_bytes()).hexdigest()==row['cooling_sha256']
 output=w/(row['element']+'_'+str(row['stage'])+'_installed.txt')
 subprocess.run([str(pairexe),str(ROOT/'data/atomic'),row['element'],str(row['stage']),str(output)],cwd=ROOT,check=True,timeout=120)
 values=np.loadtxt(output);assert values.shape==(3,3) and np.isfinite(values).all() and np.all(values[:,1:]>0)
 err=float(max(abs(values[:,1]/values[:,2]-1)));assert err<.005
 pairs.append(dict(element=row['element'],stage=row['stage'],nlev=row['nlev'],samples=3,maximum_low_density_fit_emission_difference=err))
result=dict(passed=True,shell_samples=sum(c['samples'] for c in checks),unsupported_vacancy_rejections=guard_checks,installed_atom_pairs=6,installed_reader_samples=18,cases=checks,pairs=pairs,
 scope='Actual namelist, direct cross section, packet builder, cached opacity, scored ionization and heat; installed atom/fit readers. No spatial model or vacancy cascade.',
 source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),ROOT/'src/define.f90',source,driver,pairdriver,installation,exe,pairexe,private]})
(w/'result.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS: approved shell selection and installed paired-data consumers')
