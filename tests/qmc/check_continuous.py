#!/usr/bin/env python3
"""Validate the fixed-state ensemble and measure error, variance, and cost."""
import json
from pathlib import Path
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from continuous_reference import cell_average

plt.rcParams.update({'text.usetex':True,'font.family':'serif'})
OUT=Path(sys.argv[1]).resolve()


def read(path):
    with h5py.File(path) as f:
        return {k:f[k]['data'][:] for k in ('LeafXYZ','LeafSize','Gamma_HI','T_e','n_e','x_HI','x_HeI')}


def main():
    first=read(OUT/'mc_p17_11_rates.h5')
    xyz=first['LeafXYZ']; size=first['LeafSize']; r=np.linalg.norm(xyz,axis=0)
    select=(r>.15)&(r<.85)
    args=(xyz[:,select],size[select])
    reference=cell_average(*args)
    refined_energy=cell_average(*args,energy_order=192)
    refined_space=cell_average(*args,spatial_order=8,radial_nodes=16385)
    errors={'energy_max_relative':float(np.max(abs(refined_energy/reference-1))),
            'space_and_interpolation_max_relative':float(np.max(abs(refined_space/reference-1)))}
    assert max(errors.values()) < 2.e-6, errors
    reference=refined_space
    edges=np.linspace(.15,.85,29); center=(edges[:-1]+edges[1:])/2
    bin_index=np.digitize(r[select],edges)-1
    def profile(a):
        return np.array([np.mean(a[...,bin_index==j],axis=-1) for j in range(len(center))]).T
    reference_profile=profile(reference)
    report={'ok':False,'reference_checks':errors,'selected_cells':int(select.sum()),
        'conditions':{'grid':'32^3 octree, cell size 0.0625 pc; sphere selected by cell centers',
        'nH_cm3':1.,'sphere_radius_pc':1.,'T_K':1.e4,'Tstar_K':1.e5,
        'luminosity_erg_s':1.e38,'energy_eV':[13.598,100.],
        'xHI':.1,'xHeI':1.,'xHeII':0.,'yHe':.1,'velocity':'static, no velocity field',
        'gravity':'not modeled','source':'isotropic central point, Planck energy spectrum',
        'boundary':'open cube [-1,1]^3 pc','dust':False,'metals':False,'diffuse':False,
        'secondary_ionization':False,'gas_iterations':0,'mpi_ranks':8,
        'threads_each_rank':1,'diagnostic_energy_bins':16,'comparison_shell_pc':[.15,.85]},
        'statistics_note':'Four independent seeds at each count. Mean spatial bias is not a cell-wise bias bound.',
        'groups':{}}
    profiles={}
    for power in range(17,21):
        groups={}
        for mode,seeds in [('mc',(11,22,33,44)),('qmc',(101,202,303,404))]:
            values=[]; times=[]; process_times=[]; hashes=[]
            for seed in seeds:
                tag=f'{mode}_p{power}_{seed}'
                meta=json.loads((OUT/f'{tag}.json').read_text());assert meta['ok']
                current=read(OUT/f'{tag}_rates.h5')
                for field in ('LeafXYZ','LeafSize','T_e','n_e','x_HI','x_HeI'):
                    assert np.array_equal(current[field],first[field]), (tag,field)
                values.append(current['Gamma_HI'][select]);times.append(meta['elapsed_seconds'])
                process_times.append(meta['program_seconds']);hashes.append(meta['binary_sha256'])
            assert len(set(hashes))==1
            values=np.array(values); relative=values/reference-1
            mean_errors=relative.mean(axis=1)
            mse=float(np.mean(relative**2));variance=float(np.mean(np.var(relative,axis=0,ddof=1)))
            radial=profile(values)/reference_profile-1
            groups[mode]={'N':2**power,'seeds':list(seeds),'rms_relative':mse**.5,
                'mse_relative':mse,'sampling_variance_relative':variance,
                'ensemble_mean_error_rms':float(np.sqrt(np.mean(relative.mean(axis=0)**2))),
                'spatial_mean_bias':float(mean_errors.mean()),
                'spatial_mean_bias_standard_error':float(mean_errors.std(ddof=1)/2),
                'replicate_spatial_mean_bias':mean_errors.tolist(),
                'radial_rms_relative':float(np.sqrt(np.mean(radial**2))),
                'elapsed_seconds_mean':float(np.mean(times)),
                'elapsed_seconds':times,'program_seconds_mean':float(np.mean(process_times))}
            profiles[f'{mode}_{power}']=profile(values)
        groups['mse_ratio']=groups['mc']['mse_relative']/groups['qmc']['mse_relative']
        groups['variance_ratio']=groups['mc']['sampling_variance_relative']/groups['qmc']['sampling_variance_relative']
        groups['equal_time_mse_efficiency']=groups['mse_ratio']*groups['mc']['elapsed_seconds_mean']/groups['qmc']['elapsed_seconds_mean']
        report['groups'][str(power)]=groups
    for mode in ('mc','qmc'):
        report[f'{mode}_rms_log_slope']=float(np.polyfit(np.log([2**p for p in range(17,21)]),
            np.log([report['groups'][str(p)][mode]['rms_relative'] for p in range(17,21)]),1)[0])
    np.savez(OUT/'profiles.npz',radius_pc=center,reference=reference_profile,**profiles)
    fig,axes=plt.subplots(2,3,figsize=(12,7.5))
    for mode,color in [('mc','C0'),('qmc','C3')]:
        p=profiles[f'{mode}_20'];mean=p.mean(axis=0);scatter=p.std(axis=0,ddof=1)
        axes[0,0].plot(center,mean,label=mode.upper(),color=color)
        axes[0,0].fill_between(center,mean-scatter,mean+scatter,alpha=.2,color=color)
        axes[0,1].plot(center,100*(mean/reference_profile-1),color=color,label=mode.upper())
        axes[0,1].fill_between(center,100*((mean-scatter)/reference_profile-1),
            100*((mean+scatter)/reference_profile-1),color=color,alpha=.2)
    axes[0,0].plot(center,reference_profile,'k--',label='Continuous reference')
    axes[0,0].set_yscale('log');axes[0,0].set_ylabel(r'$\Gamma_{\rm HI}\ [{\rm s}^{-1}]$')
    axes[0,0].legend();axes[0,1].set_ylabel('Relative error [percent]')
    axes[0,1].axhline(0,color='k',lw=.5)
    full_edges=np.linspace(0,1.75,57); full_r=(full_edges[:-1]+full_edges[1:])/2
    bi=np.digitize(r,full_edges)-1
    for ax,key,label,scale in [(axes[0,2],'T_e',r'$T_e\ [{\rm K}]$',1.),
                               (axes[1,0],'n_e',r'$n_{\rm H}\ [{\rm cm}^{-3}]$',1/.9),
                               (axes[1,1],'x_HI',r'$x_{\rm HI}$',1.),
                               (axes[1,2],'x_HeI',r'$x_{\rm HeI}$',1.)]:
        avg=np.array([np.mean(first[key][bi==i])*scale if np.any(bi==i) else np.nan for i in range(len(full_r))])
        ax.plot(full_r,avg,'k-');ax.set_ylabel(label)
    for ax in axes.ravel():ax.set_xlabel(r'$r\ [{\rm pc}]$');ax.grid(alpha=.25)
    fig.suptitle(r'Fixed state, $N=2^{20}$; shaded bands: scatter of four realizations')
    fig.tight_layout();fig.savefig(OUT/'profiles.pdf');plt.close(fig)
    report['ok']=True
    (OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    for p,g in report['groups'].items():
        print(p, 'MC/RQMC RMS [%]:',100*g['mc']['rms_relative'],100*g['qmc']['rms_relative'],
              'MSE ratio:',g['mse_ratio'],'variance ratio:',g['variance_ratio'],flush=True)


if __name__=='__main__':main()
