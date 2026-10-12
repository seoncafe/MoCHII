"""Plot preserved isolated thermal/emission checks without rerunning native code."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from mocassin_excitation_profiles import source_atom
from validate_excitation_density_profiles import full_excitation_balance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--installation', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    result = json.loads((work/'result.json').read_text())
    if not result['passed']:
        raise ValueError('No completed validation')
    output = work/'thermal_emission_profiles.pdf'
    if output.exists():
        raise FileExistsError('Preserve reviewed figures')
    records = json.loads(args.installation.read_text())['atoms']
    sources = [args.installation,work/'result.json',Path(__file__)]
    densities = [1e-6,100.,1e6,1e10]
    with PdfPages(output) as pdf:
        for record in records:
            name=record['element']+'_'+str(record['stage'])
            local_path=work/(name+'_local_ne.txt')
            low_path=work/(name+'_low_density.txt')
            local=np.loadtxt(local_path);low=np.loadtxt(low_path)
            sources.extend((local_path,low_path))
            t=local[::4,0]
            reference=full_excitation_balance(*source_atom(record),t,densities)
            fit=work/(name+'_atomic')/('cooling_'+name+'.txt')
            rows=[s.split() for s in fit.read_text().splitlines()
                  if s.strip() and not s.lstrip().startswith('#')]
            coefficients=np.array(rows[1:],float)
            fitted=np.sum(coefficients[:,0,None]*np.exp(-coefficients[:,1,None]/t),axis=0)/np.sqrt(t)
            fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
            for j,ne in enumerate(densities):
                color=f'C{j}';label=r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$'
                power=reference[ne][0]
                expected=fitted*power/reference[1e-6][0]
                axes[0,0].loglog(t,local[j::4,4],color=color,label=label)
                axes[0,0].loglog(t,expected,':',color=color,lw=.8)
                axes[1,0].semilogx(t,local[j::4,4]/expected-1,color=color)
                axes[0,1].loglog(t,local[j::4,5],color=color)
                axes[0,1].loglog(t,power,':',color=color,lw=.8)
                axes[1,1].loglog(t,local[j::4,3],color=color)
                axes[1,1].loglog(t,low[j::4,3],':',color=color,lw=.8)
            axes[0,0].loglog(t,low[::4,4],'k--',label='Low-density fit')
            axes[0,0].legend(fontsize=7)
            axes[0,0].set_ylabel(r'Thermal coefficient [erg cm$^3$ s$^{-1}$]')
            axes[0,1].set_ylabel(r'Emitted coefficient [erg cm$^3$ s$^{-1}$]')
            axes[1,0].set_ylabel('Relative thermal interpolation difference')
            axes[1,1].set_ylabel('Selected ion fraction')
            for ax in axes.flat:ax.set_xlabel(r'$T$ [K]')
            fig.suptitle(name.replace('_',' stage ')+'; native solid / independent dotted')
            pdf.savefig(fig)
            plt.close(fig)
    (work/'figure_record.json').write_text(json.dumps(dict(
        pages=6,native_profiles_regenerated=False,
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        figure_sha256=hashlib.sha256(output.read_bytes()).hexdigest()),indent=2)+'\n')
    print('Created six-page figure from preserved native profiles.')


if __name__ == '__main__':
    main()
