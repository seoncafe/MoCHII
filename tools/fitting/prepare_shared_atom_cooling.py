"""Prepare matched source atoms and cooling fits without production changes.

Reuse the measured registry density ratios and independently validated native
fit files. Derived thermal profiles are predictions, not a new invoked
production calculation or a physical accuracy certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

ROOT=Path(__file__).resolve().parents[2]


def coefficient(path,temperature):
    rows=[s.split() for s in path.read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
    values=np.array(rows[1:],float)
    if values.shape!=(int(rows[0][0]),2):raise ValueError('Incomplete cooling fit')
    return np.sum(values[:,0,None]*np.exp(-values[:,1,None]/temperature),axis=0)/np.sqrt(temperature)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    work=parser.parse_args().work_dir.resolve()
    if work.exists():raise FileExistsError('Preserve previous review products')
    work.mkdir(parents=True)
    base=ROOT/'md/atomic_stage_expansion/validation_2026-10-11'
    preparation=base/'shared_atom_expansion_candidate/installation/result.json'
    atoms=json.loads(preparation.read_text())['atoms']
    fit_manifest=base/'mocassin_cooling_native/manifest.json'
    fit_native=base/'mocassin_cooling_native_exact_endpoints/native/result.json'
    checks=json.loads(fit_native.read_text())
    if not checks['passed']:raise ValueError('Native fit reader checks incomplete')
    verified={r['ion']:r for r in checks['ions']}
    fits={r['ion']:r for r in json.loads(fit_manifest.read_text())['fits']}
    native=base/'cooling_temperature_pchip/native'
    if not json.loads((native/'result.json').read_text())['passed']:
        raise ValueError('No measured registry consumer results')
    previous=work/'previous_cooling';previous.mkdir()
    selected=work/'matched_data';selected.mkdir()
    sources=[Path(__file__),preparation,fit_manifest,fit_native,native/'result.json']
    records=[]
    with PdfPages(work/'matched_cooling_profiles.pdf') as pdf:
        for atom in atoms:
            name=atom['element']+'_'+str(atom['stage'])
            fit=fits[name];candidate=Path(fit['path'])
            digest=hashlib.sha256(candidate.read_bytes()).hexdigest()
            if digest!=fit['sha256'] or digest!=verified[name]['coefficients_sha256']:
                raise ValueError('Candidate fit changed after native checks')
            old=ROOT/'data/atomic'/candidate.name
            shutil.copyfile(old,previous/old.name)
            shutil.copyfile(candidate,selected/candidate.name)
            atom_path=Path(atom['selected_path'])
            shutil.copyfile(atom_path,selected/atom_path.name)
            data_path=native/(name+'_local_ne.txt')
            data=np.loadtxt(data_path)
            if data.shape!=(132,6):raise ValueError('Incomplete measured consumer profile')
            t=data[::4,0]
            if t.min()<1000 or t.max()>fit['temperature_maximum_K']:
                raise ValueError('Prediction exceeds fitted interval')
            old_fit=coefficient(old,t);new_fit=coefficient(candidate,t)
            if np.any(old_fit<=0) or np.any(new_fit<=0):raise ValueError('Nonpositive fitted cooling')
            predictions=data.copy();mismatch=0.;old_mismatch=0.
            fig,axes=plt.subplots(2,1,figsize=(7,7),sharex=True,constrained_layout=True)
            for j,ne in enumerate((1e-6,100.,1e6,1e10)):
                predicted=data[j::4,4]*(new_fit/old_fit)
                direct=data[j::4,5]
                error=predicted/direct-1
                mismatch=max(mismatch,float(max(abs(error))))
                old_mismatch=max(old_mismatch,float(max(abs(data[j::4,4]/direct-1))))
                predictions[j::4,4]=predicted
                predictions[j::4,2]=data[j::4,2]*(new_fit/old_fit)
                label=r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$'
                line=axes[0].loglog(t,predicted,label=label)[0]
                axes[0].loglog(t,direct,':',color=line.get_color(),lw=.8)
                axes[1].semilogx(t,error,color=line.get_color())
            # Match each source atom's thermal and radiative model to within
            # numerical representation errors; do not hide a source mismatch.
            if mismatch>.005:raise ValueError('Matched cooling prediction exceeds 0.5 percent')
            axes[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]');axes[0].legend(fontsize=8)
            axes[0].set_title(name.replace('_',' stage ')+': predicted cooling / direct emitted energy')
            axes[1].set_xlabel(r'$T$ [K]');axes[1].set_ylabel('Relative energy difference')
            pdf.savefig(fig);plt.close(fig)
            profile=work/(name+'_derived_profiles.txt');np.savetxt(profile,predictions,fmt='%.17e')
            records.append(dict(element=atom['element'],stage=atom['stage'],
                nlev=atom['nlev'],source_atom_sha256=atom['native_sha256'],
                selected_atom_path=str(selected/atom_path.name),selected_cooling_path=str(selected/candidate.name),
                cooling_sha256=digest,previous_cooling_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
                cooling_maximum_K=fit['temperature_maximum_K'],
                maximum_mixed_source_relative_energy_difference=old_mismatch,
                maximum_matched_source_relative_energy_difference=mismatch))
            sources.extend((old,candidate,atom_path,data_path,profile))
    result=dict(passed=True,production_modified=False,profile_values_derived=True,
        atoms=records,selected_atoms=6,selected_cooling_fits=6,predicted_samples=792,
        scope='Measured registry ratios times separately validated matching fits; no new invoked combined production calculation',
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Prepared six matched atom/fit pairs and derived full cooling profiles; production unchanged.')


if __name__=='__main__':main()
