"""Full temperature and density profiles against supplied excitation atoms.

Independent full level matrix: original spline collision strengths, physical
energy gaps, spontaneous emission and Maxwellian detailed balance. No
photoexcitation, radiation trapping or ionization from excited levels.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import chianti_cooling as C
from validate_icft_atoms import atom
from atomic_output import atomic_write_text


def full_excitation_balance(levels,radiation,collisions,temperature,density):
    ids=sorted(levels);pos={i:j for j,i in enumerate(ids)};n=len(ids)
    rc=np.zeros((len(temperature),n,n));rr=np.zeros((n,n));power=np.zeros(n)
    for l,u,w,g,a in radiation:
        if a<=0:continue
        i,j=pos[l],pos[u];gap=(levels[u][1]-levels[l][1])*C.HC_ERG_CM
        if gap<=0:raise ValueError('Nonphysical emission energy')
        rr[i,j]+=a;rr[j,j]-=a;power[j]+=a*gap
    for tr in collisions:
        l,u=sorted((tr['ll'],tr['ul']),key=lambda k:(levels[k][1],k));i,j=pos[l],pos[u]
        gap=C.transition_energy_erg(levels,tr)
        common=C.COLL_PREF*C.upsilon(tr,temperature)/np.sqrt(temperature)
        upward=common/levels[l][0]*np.exp(-gap/(C.K_B_ERG*temperature))
        downward=common/levels[u][0]
        rc[:,j,i]+=upward;rc[:,i,i]-=upward
        rc[:,i,j]+=downward;rc[:,j,j]-=downward
    result={}
    for ne in density:
        matrix=ne*rc+rr
        scale=np.max(np.abs(matrix),axis=2)
        matrix/=np.where(scale>0,scale,1)[:,:,None]
        matrix[:,0,:]=1;rhs=np.zeros((len(temperature),n));rhs[:,0]=1
        pop=np.linalg.solve(matrix,rhs[...,None])[...,0]
        if not np.all(np.isfinite(pop)) or np.any(pop < -64*np.finfo(float).eps):
            raise ValueError('Invalid full excitation populations')
        pop=np.maximum(pop,0);pop/=pop.sum(axis=1)[:,None]
        result[ne]=(pop@power/ne,pop)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work-dir',required=True,type=Path)
    args=p.parse_args();work=args.work_dir.resolve();result={}
    import matplotlib.pyplot as plt
    for name,el,st in [('OV','o',5),('ClVI','cl',6)]:
        rows=np.loadtxt(work/'products'/name/'atom_native.txt');temperature=rows[::3,0]
        density=[1e-6,100.,1e6];lev,rad,coll=atom(el,st)
        source=full_excitation_balance(lev,rad,coll,temperature,density)
        fig,axes=plt.subplots(2,1,figsize=(6,6),sharex=True,constrained_layout=True)
        differences=[];saved={'temperature_K':temperature}
        for j,ne in enumerate(density):
            full,pop=source[ne];measured=rows[j::3,3];difference=measured/full-1
            if not np.all(np.isfinite(full)) or np.any(full<=0):raise ValueError('Invalid full source cooling')
            differences.append(float(np.max(np.abs(difference))))
            label=r'$n_e = '+f'{ne:g}'+r'$ cm$^{-3}$'
            line=axes[0].loglog(temperature,full,label=label)[0]
            axes[0].loglog(temperature,measured,'--',color=line.get_color())
            axes[1].semilogx(temperature,100*difference,color=line.get_color())
            saved[f'cooling_full_density_{j}']=full;saved[f'populations_full_density_{j}']=pop
        if max(differences)>.005:raise ValueError(f'{name}: full temperature/density convergence {differences}')
        axes[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]');axes[0].set_title(name+' (solid: full source; dashed: 29 states)')
        axes[0].legend(fontsize=8);axes[1].set_xlabel(r'$T$ [K]');axes[1].set_ylabel(r'Difference [\%]')
        fig.savefig(work/'products'/name/(name+'_density_profiles.pdf'));plt.close(fig)
        np.savez_compressed(work/'products'/name/'full_density_profiles.npz',**saved)
        result[name]=dict(full_levels=len(lev),selected_levels=29,temperatures=801,densities_cm3=density,maximum_relative_differences=differences,passed=True)
    atomic_write_text(work/'density_profiles_result.json',json.dumps(dict(passed=True,ions=result),indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
