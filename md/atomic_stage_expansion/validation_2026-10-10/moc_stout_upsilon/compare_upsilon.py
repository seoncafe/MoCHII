"""Compare invoked MOCASSIN v2 Upsilon interpolation with selected Stout data."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator


HERE = Path(__file__).resolve().parent
BUILD = Path('/tmp/mochii_moc_stout_upsilon')
MOC_ROOT = Path('/home/kiseon/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2')
STOUT_ROOT = Path('/home/kiseon/CLOUDY/c25.00/data/stout')
IONS = {'Ar II': ('arii', 'ar/ar_2/ar_2'), 'Ca IV': ('caiv', 'ca/ca_4/ca_4')}


def source_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_moc(path):
    rows = path.read_text().splitlines()
    n_comments = int(rows[0])
    n_levels, n_temperature, _ = map(int, rows[1+n_comments].split())
    assert (n_levels, n_temperature) == (2, 11)
    start = 1+n_comments+1+n_levels
    temperature = np.array([float(row) for row in rows[start:start+11]])
    assert rows[start+11].strip() == '0'
    strength = np.array([float(row.split()[-1]) for row in rows[start+12:start+23]])
    assert len(strength) == 11 and rows[start+23].split()[-1] == '0'
    a_fields = rows[start+24].split()
    assert a_fields[:2] == ['1', '2']
    a_value = float(a_fields[2])
    level_fields = [rows[start+25+i].split() for i in range(2)]
    weights = [float(fields[1]) for fields in level_fields]
    energies = [float(fields[2]) for fields in level_fields]
    return temperature, strength, weights, energies, a_value


def parse_stout(stem):
    coll = Path(f'{stem}.coll')
    nrg = Path(f'{stem}.nrg')
    tp = Path(f'{stem}.tp')
    rows = [row.split() for row in coll.read_text().splitlines()]
    temperature = np.array([float(x) for x in next(row for row in rows if row[0] == 'TEMP')[1:]])
    strength = np.array([float(x) for x in next(row for row in rows if row[:4] == ['CS','ELECTRON','1','2'])[4:]])
    levels = [row.split() for row in nrg.read_text().splitlines()[1:3]]
    energies = [float(row[1]) for row in levels]
    weights = [float(row[2]) for row in levels]
    a_rows = [row.split() for row in tp.read_text().splitlines() if row.split()[:3] == ['A','1','2']]
    assert len(a_rows) == 1
    a_value = float(a_rows[0][3])
    return temperature, strength, weights, energies, a_value, (nrg,tp,coll)


def log_pchip(temperature, strength, query):
    curve = PchipInterpolator(np.log(temperature), np.log(strength), extrapolate=False)
    value = np.exp(curve(np.log(np.maximum(query, temperature[0]))))
    return value


def run_moc(executable, temperature, strength, query, name):
    node_file = BUILD/f'{name}_moc_nodes.txt'
    query_file = BUILD/f'{name}_queries.txt'
    result_file = BUILD/f'{name}_moc_upsilon.txt'
    node_file.write_text(''.join(f'{t:.9g} {y:.9g}\n' for t,y in zip(temperature,strength)))
    query_file.write_text(''.join(f'{t:.17g}\n' for t in query))
    subprocess.run([str(executable),str(node_file),str(query_file),str(result_file)],check=True)
    return np.loadtxt(result_file)[:,1]


def maximum_difference(temperature, reference, comparison, mask):
    relative = np.abs(comparison[mask]/reference[mask]-1)
    return {'percent':float(100*np.max(relative)), 'temperature_K':float(temperature[mask][np.argmax(relative)])}


if __name__ == '__main__':
    os.environ.setdefault('MPLCONFIGDIR','/tmp/mochii_mpl_moc_stout')
    BUILD.mkdir(parents=True,exist_ok=True)
    module_file = MOC_ROOT/'source/interpolation_mod.f90'
    command = ['mpiifort','-O0','-no-ipo','-module',str(BUILD),str(module_file),
               str(HERE/'fortran_spline_values.f90'),'-o',str(BUILD/'moc_spline.x')]
    compile_result = subprocess.run(command,text=True,capture_output=True,check=True)
    (HERE/'compile.log').write_text(compile_result.stdout+compile_result.stderr)
    np.set_printoptions(precision=8)
    figure, axes = plt.subplots(2,2,figsize=(11,8.5),sharex='col')
    results = {'method': {'moc':'Original MOCASSIN v2 interpolation_mod.f90 natural cubic on log10(T), real32, constant outside nodes',
                          'stout':'MoCHII log-PCHIP on log(T) and log(Upsilon), constant below first node'},
               'compiler_command':command, 'source_sha256':{str(module_file):source_hash(module_file)}}
    for column,(ion,(moc_name,stout_name)) in enumerate(IONS.items()):
        moc_path = MOC_ROOT/'data'/f'{moc_name}.dat'
        stout_stem = STOUT_ROOT/stout_name
        moc_t,moc_y,moc_g,moc_e,moc_a = parse_moc(moc_path)
        stout_t,stout_y,stout_g,stout_e,stout_a,stout_files = parse_stout(stout_stem)
        assert np.array_equal(moc_y,stout_y)
        assert moc_g == stout_g == [4.,2.]
        assert moc_a == stout_a
        all_t = np.r_[np.geomspace(1e3,1e5,801),moc_t,stout_t,
                      moc_t[0]*(1+np.array([-1e-5,0,1e-5])),
                      stout_t[0]*(1+np.array([-1e-5,0,1e-5]))]
        temperature = np.unique(all_t[(all_t>=1e3)&(all_t<=1e5)])
        moc_actual = run_moc(BUILD/'moc_spline.x',moc_t,moc_y,temperature,moc_name)
        stout_actual = log_pchip(stout_t,stout_y,temperature)
        same_coordinate_pchip = log_pchip(moc_t,moc_y,temperature)
        actual_error = 100*(stout_actual/moc_actual-1)
        interpolation_error = 100*(same_coordinate_pchip/moc_actual-1)
        coordinate_error = 100*(stout_actual/same_coordinate_pchip-1)
        sample = np.column_stack((temperature,moc_actual,stout_actual,same_coordinate_pchip,
                                  actual_error,interpolation_error,coordinate_error))
        np.savetxt(HERE/f'{moc_name}_full_profile.csv',sample,delimiter=',',
                   header='T_K,MOC_v2_Upsilon,Stout_logPCHIP_Upsilon,logPCHIP_on_MOC_nodes,actual_percent,interpolation_percent,node_coordinate_percent',comments='')
        ax=axes[0,column]
        ax.plot(temperature,moc_actual,label='MOCASSIN v2 natural cubic',lw=1.4)
        ax.plot(temperature,stout_actual,label='Stout log-PCHIP',lw=1.1)
        plotted_nodes=moc_t<=1e5
        ax.scatter(moc_t[plotted_nodes],moc_y[plotted_nodes],s=12,color='black',label='MOCASSIN nodes')
        ax.axvline(moc_t[0],ls=':',color='gray')
        ax.axvline(stout_t[0],ls='--',color='gray')
        ax.set_xscale('log');ax.set_xlim(1e3,1e5)
        ax.set_ylabel(r'Effective collision strength $\Upsilon$')
        ax.set_title(ion);ax.legend(fontsize=7)
        ax=axes[1,column]
        ax.plot(temperature,actual_error,label='Actual source formats',lw=1.4)
        ax.plot(temperature,interpolation_error,label='Same node temperatures',lw=1.0)
        ax.plot(temperature,coordinate_error,label='Node coordinate difference',lw=1.0)
        ax.axhline(0,color='black',lw=.5)
        ax.axvline(moc_t[0],ls=':',color='gray')
        ax.axvline(stout_t[0],ls='--',color='gray')
        ax.set_xscale('log');ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('Relative difference (percent)');ax.legend(fontsize=7)
        in_range=(temperature>=moc_t[0])&(temperature<=min(1e5,moc_t[-1]))
        cold=(temperature>=1e3)&(temperature<moc_t[0])
        moc_at_nodes=run_moc(BUILD/'moc_spline.x',moc_t,moc_y,moc_t,moc_name+'_nodes')
        stout_at_nodes=log_pchip(stout_t,stout_y,stout_t)
        common_moc_nodes=moc_t<=stout_t[-1]
        stout_on_moc_nodes=log_pchip(stout_t,stout_y,moc_t[common_moc_nodes])
        results[ion]={'MOC':{'temperature_K':moc_t.tolist(),'strength':moc_y.tolist(),'weights':moc_g,
                            'energies_cm_1':moc_e,'A_s_1':moc_a},
                      'Stout':{'temperature_K':stout_t.tolist(),'strength':stout_y.tolist(),'weights':stout_g,
                               'energies_cm_1':stout_e,'A_s_1':stout_a},
                      'samples':len(temperature),'sample_range_K':[float(temperature[0]),float(temperature[-1])],
                      'max_actual_percent':maximum_difference(temperature,moc_actual,stout_actual,np.ones(len(temperature),dtype=bool)),
                      'max_in_MOC_nodes_percent':maximum_difference(temperature,moc_actual,stout_actual,in_range),
                      'max_interpolation_only_percent':maximum_difference(temperature,moc_actual,same_coordinate_pchip,in_range),
                      'max_node_coordinate_percent':maximum_difference(temperature,same_coordinate_pchip,stout_actual,np.ones(len(temperature),dtype=bool)),
                      'max_cold_percent':maximum_difference(temperature,moc_actual,stout_actual,cold),
                      'node_checks':{'MOC_max_absolute':float(np.max(np.abs(moc_at_nodes-moc_y))),
                                     'Stout_max_absolute':float(np.max(np.abs(stout_at_nodes-stout_y))),
                                     'max_cross_source_percent_at_common_MOC_nodes':maximum_difference(
                                         moc_t[common_moc_nodes],moc_at_nodes[common_moc_nodes],
                                         stout_on_moc_nodes,np.ones(np.count_nonzero(common_moc_nodes),dtype=bool))},
                      'source_sha256':{str(p):source_hash(p) for p in (moc_path,*stout_files)}}
    figure.tight_layout()
    figure.savefig(HERE/'upsilon_full_profiles.pdf')
    plt.close(figure)
    (HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n')
