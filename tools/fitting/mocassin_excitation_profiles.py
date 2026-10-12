"""Save complete independent excitation profiles and vector comparison figures."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from mocassin_v2_atom_inventory import read_atom, MOCASSIN
from validate_excitation_density_profiles import full_excitation_balance


def source_atom(record):
    source=read_atom(MOCASSIN/record['relative_file'],include_tables=True)
    levels={r['index']:(r['statistical_weight'],r['energy_cm_inverse']) for r in source['levels']}
    radiation=[(r['lower'],r['upper'],0.,0.,r['A_s_inverse'])
               for r in source['spontaneous_probabilities'] if r['A_s_inverse']>0]
    collisions=[dict(ll=r['lower'],ul=r['upper'],temperature_grid=np.array(source['temperatures_K']),
                     upsilon_grid=np.array(r['upsilon']))
                for r in source['collision_tables'] if any(r['upsilon'])]
    return levels,radiation,collisions


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate-dir',required=True,type=Path)
    p.add_argument('--native-dir',required=True,type=Path)
    p.add_argument('--output-dir',required=True,type=Path)
    a=p.parse_args();out=a.output_dir.resolve()
    if out.exists():raise FileExistsError('Preserve earlier figures')
    out.mkdir(parents=True)
    manifest=json.loads((a.candidate_dir/'manifest.json').read_text())
    densities=[1e-6,100.,1e6,1e10];pdf=out/'excitation_profiles.pdf';products=[]
    with PdfPages(pdf) as pages:
        for index,record in enumerate(manifest['atoms'],1):
            name=f"{record['element']}_{record['stage']}"
            data=np.loadtxt(a.native_dir/name/'profile.txt');T=data[::4,0]
            levels,radiation,collisions=source_atom(record)
            reference=full_excitation_balance(levels,radiation,collisions,T,densities)
            fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
            stored={'temperature_K':T,'density_cm_inverse_3':np.array(densities)}
            for j,ne in enumerate(densities):
                power,pop=reference[ne];measured=data[j::4]
                stored[f'source_cooling_{j}']=power;stored[f'source_populations_{j}']=pop
                stored[f'native_cooling_{j}']=measured[:,2];stored[f'native_populations_{j}']=measured[:,3:]
                label=r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$'
                axes.flat[0].loglog(T,measured[:,2],label=label)
                axes.flat[0].loglog(T,power,'k:',linewidth=.6)
                ax=axes.flat[j+1]
                # The display threshold hides only zero/tiny population tails;
                # all original values remain in the saved NPZ and native text.
                for k in range(record['nlev']):
                    ax.loglog(T,np.maximum(measured[:,3+k],1e-30),linewidth=.65)
                    ax.loglog(T,np.maximum(pop[:,k],1e-30),'k:',linewidth=.3)
                ax.set_title(label);ax.set_ylim(1e-30,1.1)
                ax.set_ylabel('All level populations')
                axes.flat[5].semilogx(T,(measured[:,2]/power-1),label=label)
            axes.flat[0].set_ylabel(r'$Lambda$ [erg cm$^3$ s$^{-1}$]');axes.flat[0].legend(fontsize=7)
            axes.flat[5].set_ylabel('Relative cooling difference')
            for ax in axes.flat:ax.set_xlabel(r'$T$ [K]')
            fig.suptitle(name+f": {record['nlev']} levels; native solid / independent dotted")
            pages.savefig(fig);plt.close(fig)
            path=out/(name+'.npz');np.savez_compressed(path,**stored)
            products.append(path)
            print(f'Profile page {index}/54: {name}',flush=True)
    products.append(pdf)
    result={'passed':True,'pages':54,'visual_review_completed':False,
            'population_display_floor':1e-30,'products_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in products},
            'scope':'All source/native level populations at four densities and full sampled temperature curves; figure generation does not mean visual review'}
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
