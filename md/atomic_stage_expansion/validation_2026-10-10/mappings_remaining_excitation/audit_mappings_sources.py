"""Inspect MAPPINGS inputs and isolate its Cl I statistical-weight defect.

External files are read unchanged. The comparison holds the source's unknown
constant effective strength fixed; it measures only the weight substitution,
not the accuracy of Cl I electron excitation or a complete plasma model.
"""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
MAPPINGS = Path.home() / 'RT_Codes/MAPPINGS/mappings_V-521'
LINES = MAPPINGS / 'lab/data/lines'
STOUT = Path.home() / 'CLOUDY/c25.00/data/stout/cl/cl_1/cl_1.nrg'
TARGET = {(17,1),(17,5),(20,1),(20,3),(18,1),(26,1),(12,1),(14,1)}


def records(path):
    return [line.strip() for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('%')]


def main():
    output = HERE / 'result.json'
    if output.exists():
        raise FileExistsError('Preserve existing results')
    multi = iter(records(LINES/'MULTIDAT.txt'))
    title, count = next(multi), int(next(multi))
    ions = []
    for expected in range(1, count+1):
        index = int(next(multi))
        assert index == expected
        z, stage = map(int, next(multi).split())
        nlevels = int(next(multi))
        levels = [next(multi) for _ in range(nlevels)]
        ntrans = int(next(multi))
        types = []
        for _ in range(ntrans):
            row = next(multi).split()
            kind, nodes = int(row[2]), int(row[3])
            types.append(kind)
            rows = {0:0,1:1,2:1,3:1,6:1,13:1,4:2,7:2,5:3,8:3}[kind]
            for _ in range(rows):
                assert len(next(multi).split()) == nodes
        ions.append({'Z':z,'stage':stage,'nlevels':nlevels,'ntrans':ntrans,
                     'types':sorted(set(types))})
    assert not list(multi)
    inter = records(LINES/'INTERDAT2.txt')
    inter_count = int(inter[1])
    assert len(inter[2:]) == inter_count
    selected = [line for line in inter[2:] if tuple(map(int,line.split()[:2])) in TARGET]
    cl = [line.split() for line in selected if line.startswith('17 01')]
    cl = [r for r in cl if float(r[2]) < 1e-12]
    assert len(cl) == 1
    energy,g1,g2,A,strength = map(float,cl[0][2:])
    observed = {}
    for line in STOUT.read_text().splitlines():
        fields = line.split()
        if len(fields)>=3 and fields[0] in ('1','2'):
            observed[int(fields[0])] = (float(fields[1]),float(fields[2]))
    assert observed == {1:(0.0,4.0),2:(882.351,2.0)}
    assert (g1,g2) == (2.0,4.0)
    T = np.geomspace(1000,100000,501)
    density = np.asarray([1e-6,100.,1e6,1e10])
    def cooling(weights):
        q12 = 8.629133680865427904126e-6/np.sqrt(T)*strength*np.exp(-energy/(1.380649e-16*T))/weights[0]
        q21 = 8.629133680865427904126e-6/np.sqrt(T)*strength/weights[1]
        population = density[:,None]*q12/(A+density[:,None]*(q12+q21))
        return population*A*energy
    original,physical = cooling((g1,g2)),cooling((4.,2.))
    ratio = original/physical
    assert np.all(np.isfinite(ratio)) and np.all((ratio>=2-1e-12)&(ratio<=4+1e-12))
    np.savez(HERE/'cl_i_weight_profiles.npz',temperature_K=T,electron_density_cm3=density,
             source_cooling_erg_s_ion=original,physical_weight_cooling_erg_s_ion=physical,
             source_to_physical_weight_ratio=ratio)
    fig,axes = plt.subplots(2,2,figsize=(9,6),sharex=True)
    for axis,n,curve in zip(axes.flat,density,ratio):
        axis.semilogx(T,curve)
        axis.set_title(r'$n_e = '+f'{n:.0e}'+r'\ {\rm cm}^{-3}$')
        axis.set_xlabel('Temperature (K)')
        axis.set_ylabel('Source / physical weights')
        axis.set_ylim(1.9,4.1)
    fig.tight_layout();fig.savefig(HERE/'cl_i_weight_profiles.pdf');plt.close(fig)
    paths = [Path(__file__),STOUT,LINES/'MULTIDAT.txt',LINES/'INTERDAT2.txt',
             LINES/'RESONDAT2.txt',MAPPINGS/'src/mapinit.f',MAPPINGS/'src/inter.f',
             MAPPINGS/'src/cool.f',MAPPINGS/'src/reson.f',MAPPINGS/'src/const.inc']
    result = {'multilevel_title':title,'multilevel_ions':count,
              'target_ions_in_multilevel':[r for r in ions if (r['Z'],r['stage']) in TARGET],
              'intercombination_rows':inter_count,'target_rows':selected,
              'Cl_I_source_weights':[g1,g2],'Cl_I_physical_weights':[4,2],
              'Cl_I_constant_strength':strength,'Cl_I_A_s_inverse':A,
              'Cl_I_energy_erg':energy,'ratio_min':float(ratio.min()),'ratio_max':float(ratio.max()),
              'assumptions':['Maxwellian electrons; kappa correction disabled.',
                             'Identical source constant strength and energy in both calculations.',
                             'No acceptance of source strength accuracy or temperature coverage.',
                             'Isolated two-level model; no photons or other collisions.'],
              'external_sources_changed':False,'production_changed':False,
              'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('multilevel_ions','target_ions_in_multilevel','ratio_min','ratio_max')}))


if __name__ == '__main__':
    main()
