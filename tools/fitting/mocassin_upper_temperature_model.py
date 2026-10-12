"""Diagnose qualified constant upper collision strengths in isolated modules.

Historical record (2026-10-11 diagnosis of the upper collision-table
ceiling). It cannot run on the current tree: the nlevel_mod expression it
replaces ('temperature exceeds collision table') no longer exists, because
holding Upsilon at the final table node is now the production rule
(nlevel_mod::nlevel_upsilon).

Original nodes and their log-PCHIP interpolation are retained. Only Upsilon
is held at the final node; Maxwellian rates use actual T. The 100000 K model
ceiling is a proposed numerical scope, not a verified atomic accuracy bound.
No production module, data format, input parameter, or atom is changed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.interpolate import PchipInterpolator
import chianti_cooling as C
from mocassin_excitation_profiles import source_atom
from validate_source_excitation_atom import compile_driver, grid_file, DENSITY

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'md/atomic_stage_expansion/phase1_77_mocassin_v2/candidate_atoms'
TARGETS = {('c',1), ('s',3), ('ar',6)}


def electron_excitation_balance_constant_endpoints(levels,radiation,collisions,T):
    ids=sorted(levels);position={k:i for i,k in enumerate(ids)};n=len(ids)
    coll=np.zeros((len(T),n,n));rad=np.zeros((n,n));power=np.zeros(n)
    for l,u,w,g,a in radiation:
        i,j=position[l],position[u];gap=(levels[u][1]-levels[l][1])*C.HC_ERG_CM
        if gap<=0:raise ValueError('Nonpositive photon energy')
        rad[i,j]+=a;rad[j,j]-=a;power[j]+=a*gap
    for row in collisions:
        l,u=sorted((row['ll'],row['ul']),key=lambda k:(levels[k][1],k))
        i,j=position[l],position[u];gap=(levels[u][1]-levels[l][1])*C.HC_ERG_CM
        strength=C.upsilon(row,np.minimum(T,row['temperature_grid'][-1]))
        common=C.COLL_PREF*strength/np.sqrt(T)
        up=common/levels[l][0]*np.exp(-gap/(C.K_B_ERG*T));down=common/levels[u][0]
        coll[:,j,i]+=up;coll[:,i,i]-=up;coll[:,i,j]+=down;coll[:,j,j]-=down
    result={}
    for ne in DENSITY:
        matrix=ne*coll+rad
        scale=np.max(abs(matrix),axis=2)
        matrix/=np.where(scale>0,scale,1)[:,:,None]
        matrix[:,0,:]=1
        rhs=np.zeros((len(T),n));rhs[:,0]=1
        pop=np.linalg.solve(matrix,rhs[...,None])[...,0]
        if not np.isfinite(pop).all() or pop.min()<-64*np.finfo(float).eps:
            raise ValueError('Nonphysical matrix populations')
        pop=np.maximum(pop,0);pop/=pop.sum(axis=1)[:,None]
        result[ne]=(pop@power/ne,pop)
    return result


def compile_module(source,directory):
    obj=directory/(source.stem+'.o')
    subprocess.run(['mpiifort','-O0','-no-ipo','-check','all,noarg_temp_created','-fpe0',
        '-traceback','-I'+str(directory),'-Isrc','-module',str(directory),'-c',str(source),
        '-o',str(obj)],cwd=ROOT,check=True,timeout=120)
    return obj


def profiles(work):
    if work.exists():raise FileExistsError('Preserve previous model calculations')
    work.mkdir(parents=True)
    original=ROOT/'src/nlevel_mod.f90';text=original.read_text()
    previous="""       st=log(T)
       if (st > atom%node_x(atom%offset(k)+atom%nc(k)-1)) &
          error stop 'temperature exceeds collision table'
"""
    replacement="""       ! Diagnostic only: retain original PCHIP within source nodes and
       ! hold only Upsilon at the final node over T <= 100000 K.
       ! The hot model is unaccepted and its atomic accuracy is unverified.
       if (log(T) > log(100000.0_wp)) error stop 'temperature exceeds proposed endpoint model'
       st=min(log(T),atom%node_x(atom%offset(k)+atom%nc(k)-1))
"""
    if text.count(previous)!=1:raise ValueError('Temperature-table expression changed')
    private=work/'nlevel_mod.f90';private.write_text(text.replace(previous,replacement))
    obj=compile_module(private,work);exe=work/'endpoint_atom.x'
    driver=ROOT/'tests/ion_stages/check_source_excitation_atom.f90'
    subprocess.run(['mpiifort','-O0','-no-ipo','-check','all,noarg_temp_created','-fpe0',
        '-traceback','-I'+str(work),'-Isrc','-module',str(work),'-o',str(exe),str(driver),
        str(ROOT/'src/define.o'),str(obj)],cwd=ROOT,check=True,timeout=120)
    baseline_dir=work/'original';baseline_dir.mkdir();baseline=compile_driver(baseline_dir)
    manifest=BASE/'manifest.json';records=[r for r in json.loads(manifest.read_text())['atoms']
                                         if (r['element'],r['stage']) in TARGETS]
    if len(records)!=3:raise ValueError('Incomplete three-ion source scope')
    sources=[Path(__file__),original,private,obj,exe,driver,baseline,manifest,ROOT/'src/nlevel_mod.o',ROOT/'src/define.o']
    results=[]
    with PdfPages(work/'endpoint_profiles.pdf') as pdf:
        for record in records:
            name=record['element']+'_'+str(record['stage']);directory=work/name;directory.mkdir()
            atom=Path(record['path'])
            if hashlib.sha256(atom.read_bytes()).hexdigest()!=record['native_sha256']:
                raise ValueError('Source candidate changed')
            levels,radiation,collisions=source_atom(record);limit=min(row['temperature_grid'][-1] for row in collisions)
            cold=np.geomspace(1000.,limit,101);hot=np.geomspace(limit,100000.,101)
            T=np.r_[cold,hot[1:]]
            for kind,t,executable in [('cold_original',cold,baseline),('cold_endpoint',cold,exe),('combined',T,exe)]:
                grid=directory/(kind+'.in');output=directory/(kind+'.txt');grid_file(grid,2,t)
                subprocess.run([str(executable),str(atom),str(grid),str(output)],check=True,timeout=120)
            old=np.loadtxt(directory/'cold_original.txt');unchanged=np.loadtxt(directory/'cold_endpoint.txt')
            if not np.array_equal(old,unchanged):raise ValueError('Upper continuation changed in-range populations or cooling')
            data=np.loadtxt(directory/'combined.txt')
            if data.shape!=(804,3+record['nlev']) or not np.isfinite(data).all():raise ValueError('Incomplete model profiles')
            independent=electron_excitation_balance_constant_endpoints(levels,radiation,collisions,T)
            error=pop_error=0.;appended_node_error=0.
            for row in collisions:
                grid=row['temperature_grid'];strength=row['upsilon_grid']
                t=np.geomspace(grid[-2],grid[-1],201)
                original_strength=C.upsilon(row,t)
                appended=np.exp(PchipInterpolator(np.log(np.r_[grid,100000.]),
                    np.log(np.r_[strength,strength[-1]]))(np.log(t)))
                appended_node_error=max(appended_node_error,float(max(abs(appended/original_strength-1))))
            fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
            saved={'temperature_K':T,'native_profiles':data}
            for j,ne in enumerate(DENSITY):
                power,pop=independent[ne];measured=data[j::4]
                error=max(error,float(max(abs(measured[:,2]/power-1))))
                pop_error=max(pop_error,float(max(abs(measured[:,3:]-pop).ravel())))
                saved[f'independent_power_{j}']=power;saved[f'independent_populations_{j}']=pop
                label=r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$'
                line=axes.flat[0].loglog(T,measured[:,2],label=label)[0]
                axes.flat[0].loglog(T,power,':',color=line.get_color(),lw=.7)
                for k in range(record['nlev']):
                    axes.flat[j+1].loglog(T,np.maximum(measured[:,3+k],1e-30),lw=.7)
                    axes.flat[j+1].loglog(T,np.maximum(pop[:,k],1e-30),'k:',lw=.3)
                axes.flat[j+1].set_title(label);axes.flat[j+1].set_ylabel('Level populations')
                axes.flat[5].semilogx(T,measured[:,2]/power-1)
            if error>1e-8 or pop_error>1e-8:raise ValueError('Model differs from independent matrix')
            axes.flat[0].legend(fontsize=7);axes.flat[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
            axes.flat[5].set_ylabel('Relative matrix cooling difference')
            for ax in axes.flat:ax.set_xlabel(r'$T$ [K]');ax.axvline(limit,color='.5',ls='--',lw=.6)
            fig.suptitle(name+': proposed constant upper strength; native solid / matrix dotted')
            pdf.savefig(fig);plt.close(fig);np.savez_compressed(directory/'profiles.npz',**saved)
            grid=directory/'above_model.in';grid_file(grid,1,np.array([100000.1]))
            rejected=subprocess.run([str(exe),str(atom),str(grid),str(directory/'above_model.txt')],
                                    capture_output=True,text=True,timeout=30)
            (directory/'above_model.log').write_text(rejected.stdout+rejected.stderr)
            if not rejected.returncode or 'temperature exceeds proposed endpoint model' not in rejected.stdout+rejected.stderr:
                raise ValueError('Proposed model ceiling not enforced')
            results.append(dict(element=record['element'],stage=record['stage'],nlev=record['nlev'],
                source_maximum_K=float(limit),model_maximum_K=100000.,cold_samples=404,combined_samples=804,
                in_range_profiles_identical=True,model_ceiling_rejected=True,
                maximum_relative_cooling_difference=error,maximum_absolute_population_difference=pop_error,
                maximum_in_range_strength_change_if_synthetic_upper_node_appended=appended_node_error))
            sources.extend((atom,directory/'cold_original.txt',directory/'cold_endpoint.txt',directory/'combined.txt',directory/'profiles.npz'))
    result=dict(passed=True,production_modified=False,hot_model_accepted=False,atomic_accuracy_verified=False,
        models=results,model_count=3,combined_samples=2412,cold_identity_samples=1212,figure_pages=3,
        scope='Private evaluation and independent matrices with constant Upsilon endpoints and actual Maxwellian T; no production input or data change',
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: three upper-endpoint models; 2412 combined samples; 1212 original-range profiles unchanged.')


def cooling(work,prepared):
    if work.exists():raise FileExistsError('Preserve previous cooling checks')
    work.mkdir(parents=True)
    r=json.loads((prepared/'result.json').read_text())
    if not r['passed']:raise ValueError('Direct model validation incomplete')
    for path,digest in r['source_sha256'].items():
        source=Path(path)
        if source.resolve()==Path(__file__).resolve():
            source=prepared/'frozen_analysis_tool.py'
        if hashlib.sha256(source.read_bytes()).hexdigest()!=digest:
            raise ValueError('Completed model source or profile changed')
    # Keep completed profiles. The new diagnostic differs only in the upper
    # guard's floating-point coordinate: log T matches registration's guard.
    model_text=(prepared/'nlevel_mod.f90').read_text().replace(
        'if (T > 100000.0_wp)', 'if (log(T) > log(100000.0_wp))')
    model_module=work/'nlevel_mod.f90';model_module.write_text(model_text)
    model_obj=compile_module(model_module,work)
    module=ROOT/'src/nlevel_cooling_mod.f90';text=module.read_text()
    marker='    limited_grid = reg(ie)%source_ceiling_ln_K < log(10.0_wp)*gT(NT)'
    if text.count(marker)!=1:raise ValueError('Cooling source-ceiling expression changed')
    text=text.replace(marker,'    ! Diagnostic model only: unaccepted constant-Upsilon upper endpoints through 100000 K.\n    reg(ie)%source_ceiling_ln_K = log(100000.0_wp)\n'+marker)
    public='  public :: nlevel_cooling_lam_se, nlevel_cooling_report\n'
    text=text.replace(public,public+'  public :: reg\n')
    private=work/'nlevel_cooling_mod.f90';private.write_text(text)
    # Import the identical atom type from the previously compiled private module.
    flags=['mpiifort','-O0','-no-ipo','-check','all,noarg_temp_created','-fpe0','-traceback',
           '-I'+str(work),'-I'+str(prepared),'-Isrc','-module',str(work)]
    obj=work/'nlevel_cooling_mod.o'
    subprocess.run([*flags,'-c',str(private),'-o',str(obj)],cwd=ROOT,check=True,timeout=120)
    original_driver=ROOT/'tests/ion_stages/check_registered_excitation_atom.f90';driver_text=original_driver.read_text()
    start=driver_text.index('  source_ceiling_ln_K = log(1.0e6_wp)')
    end=driver_text.index('  if (source_ceiling_ln_K < log(1.0e6_wp)) then',start)
    driver_text=driver_text[:start]+'  source_ceiling_ln_K = log(1.0e5_wp)\n'+driver_text[end:]
    driver_text=driver_text.replace('nlevel_cooling_apply, nlevel_cooling_report','nlevel_cooling_apply, nlevel_cooling_report, reg')
    driver_text=driver_text.replace("'# source_ceiling_ln_K '","'# model_ceiling_ln_K '")
    marker="  write(unit, '(a,i0)') '# nlev ', atom%nlev"
    driver_text=driver_text.replace(marker,marker+"\n  write(unit, '(a,i0)') '# temperature_nodes ', size(reg(idx)%logT)\n  write(unit, '(a,i0)') '# density_nodes ', size(reg(idx)%lambda,2)")
    driver=work/'registered_model.f90';driver.write_text(driver_text);exe=work/'registered_model.x'
    subprocess.run([*flags,'-I'+str(work),'-o',str(exe),str(driver),str(ROOT/'src/define.o'),
        str(model_obj),str(obj)],cwd=ROOT,check=True,timeout=120)
    records=[]
    for model in r['models']:
        name=model['element']+'_'+str(model['stage']);directory=work/name
        subprocess.run(['python3','tools/fitting/validate_registered_excitation_atom.py',
            '--source-dir',str(BASE),'--element',model['element'],'--stage',str(model['stage']),
            '--driver',str(exe),'--work-dir',str(directory),'--timeout','120'],cwd=ROOT,check=True,timeout=180)
        checked=json.loads((directory/'result.json').read_text())
        if not checked['passed'] or checked['table_cells']!=301*161:raise ValueError('Private table shape or checks incomplete')
        records.append(checked)
    sources=[Path(__file__),module,private,model_module,model_obj,driver,original_driver,exe,
             prepared/'result.json',prepared/'frozen_analysis_tool.py']
    result=dict(passed=True,model_count=3,samples=120,production_modified=False,hot_model_accepted=False,
        models=records,scope='Private 100000 K endpoint model registration and actual cooling consumers; no production model adoption',
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: 120 private endpoint-model cooling registration samples.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step',choices=('profiles','cooling'))
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--prepared',type=Path)
    a=parser.parse_args()
    if a.step=='profiles':profiles(a.work_dir.resolve())
    elif a.prepared:cooling(a.work_dir.resolve(),a.prepared.resolve())
    else:raise ValueError('Cooling checks require completed private profiles')


if __name__=='__main__':main()
