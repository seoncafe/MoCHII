"""Audit the exact terminal reverse-rate records used by the staged diagnostic."""

import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT/'tools/fitting'))
import make_element_data as sources
import validate_terminal_ionization as diagnostic

HERE = Path(__file__).resolve().parent
MASTER = {'RR': Path('/tmp/badnell_rr_cfout_K.gz'),
          'DR': Path('/tmp/badnell_dr_cfout_K.gz')}
PAPERS = {
    1: ('H-like', 'Badnell 2006, A&A 447, 389', '10.1051/0004-6361:20053269'),
    3: ('Li-like', 'Colgan et al. 2004, A&A 417, 1183', '10.1051/0004-6361:20034197'),
    5: ('B-like', 'Altun et al. 2004, A&A 420, 775', '10.1051/0004-6361:20040964'),
    8: ('O-like', 'Zatsarinny et al. 2003, A&A 412, 587', '10.1051/0004-6361:20031462'),
    9: ('F-like', 'Zatsarinny et al. 2006, A&A 447, 379', '10.1051/0004-6361:20053737'),
    10: ('Ne-like', 'Zatsarinny et al. 2004, A&A 426, 699', '10.1051/0004-6361:20040463'),
    11: ('Na-like', 'Altun et al. 2006, A&A 447, 1165', '10.1051/0004-6361:20053278'),
    14: ('Si-like', 'Kaur et al. 2018, A&A 610, A41', '10.1051/0004-6361/201731243'),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_comparison(path):
    """Map (Z,N,M) to the authors' original/fit table and DIFF column."""
    text = gzip.decompress(path.read_bytes()).decode('ascii')
    blocks = {}
    key = None
    in_table = False
    for line in text.splitlines():
        match = re.match(r'^Z=\s*(\d+) N=\s*(\d+)(?: M=\s*(\d+))?', line)
        if match:
            z,n,m=match.groups()
            key=(int(z),int(n),int(m) if m is not None else 1)
            blocks[key] = []
            in_table = False
            continue
        if 'TEMP(K)' in line and 'DIFF' in line:
            in_table = True
            continue
        if in_table:
            fields = line.split()
            if len(fields) == 4:
                try:
                    values = tuple(float(re.sub(r'(?<=\d)([+-]\d{2,3})$', r'E\1', value))
                                   for value in fields)
                except ValueError:
                    continue
                blocks[key].append(values)
    return blocks


def dielectronic_component(row, temperature):
    if 'DR' in row:
        fields=row['DR']
        count=int(fields[0])
        c=np.asarray(fields[1:1+count])
        energy=np.asarray(fields[1+count:1+2*count])
        return temperature**(-1.5)*sum(c*np.exp(-energy/temperature))
    a,b,t0,t1=row['DR2']
    return a*temperature**(-1.5)*np.exp(-t0/temperature)*(1+b*np.exp(-t1/temperature))


def main():
    rr_table, dr_table, ci_table = (sources.parse_badnell_rr(),
                                  sources.parse_badnell_dr(),sources.parse_cf())
    originals = {kind: source_comparison(path) for kind,path in MASTER.items()}
    points = np.array([1000.,10000.,100000.])
    temperature=np.geomspace(1e3,1e5,801)
    profiles={'temperature_K':temperature}
    figure,axes=plt.subplots(4,3,figsize=(12,13),sharex=True)
    rows = []
    for axis,(element,(z,terminal)) in zip(axes.flat,diagnostic.ROSTER.items()):
        row,rr_origin,dr_origin = diagnostic.source_terminal(
            element,z,terminal,rr_table,dr_table,ci_table)
        key=(z,z-terminal,1)
        rr_key='RR' if 'RR' in row else 'RR2'
        rr = np.array([diagnostic.recombination({rr_key:row[rr_key]},t) for t in points])
        total=np.array([diagnostic.recombination(row,t) for t in points])
        dr=np.array([dielectronic_component(row,t) for t in points])
        assert np.allclose(rr+dr,total,rtol=1e-12,atol=0)
        rr_profile=np.array([diagnostic.recombination({rr_key:row[rr_key]},t) for t in temperature])
        dr_profile=np.array([dielectronic_component(row,t) for t in temperature])
        total_profile=np.array([diagnostic.recombination(row,t) for t in temperature])
        assert np.allclose(rr_profile+dr_profile,total_profile,rtol=1e-12,atol=0)
        profiles[element+'_RR_cm3_s']=rr_profile
        profiles[element+'_DR_cm3_s']=dr_profile
        profiles[element+'_total_cm3_s']=total_profile
        profiles[element+'_DR_fraction']=dr_profile/total_profile
        profiles[element+'_log_terminal_CI_over_reverse']=diagnostic.ionization_log(row,temperature)-np.log(total_profile)
        axis.loglog(temperature,rr_profile,label='RR')
        axis.loglog(temperature,np.where(dr_profile>0,dr_profile,np.nan),label='DR')
        axis.loglog(temperature,total_profile,label='RR + DR',lw=1)
        axis.set_ylim(0.3*min(np.min(rr_profile),np.min(total_profile)),
                      3.0*max(np.max(rr_profile),np.max(total_profile)))
        axis.set_title(element.upper())
        axis.set_xlabel('Temperature (K)')
        axis.set_ylabel(r'Rate coefficient (cm$^3$ s$^{-1}$)')
        record={'element':element,'terminal_stage':terminal,'recombining_stage':terminal+1,
                'source_key_Z_N_M':list(key),'rr_source':rr_origin,'dr_source':dr_origin,
                'CI_source_key_Z_N':[z,z-terminal+1],
                'CI_source':'Voronov cfit.f via parse_cf',
                'CI_coefficients':row['CI'],
                'RR_coefficients':row[rr_key],
                'DR_coefficients':row['DR'] if 'DR' in row else row['DR2'],
                'temperature_K':points.tolist(),'RR_cm3_s':rr.tolist(),'DR_cm3_s':dr.tolist(),
                'total_cm3_s':total.tolist(),'DR_fraction':(dr/total).tolist(),
                'RR_table_record':key[:2] in rr_table,'DR_table_record':key[:2] in dr_table}
        if element!='fe':
            record['DR_sequence_paper']=PAPERS[key[1]]
            record['RR_paper']=('Badnell 2006, ApJS 167, 334', '10.1086/508465') if key[1]<=11 else PAPERS[key[1]][1:]
            charge=z-key[1]
            record['source_native_tabulation_range_K']=[10*charge**2,1e7*charge**2]
            for kind in ('RR','DR'):
                entries=np.array(originals[kind][key])
                valid=entries[(entries[:,0]>=1000)&(entries[:,0]<=100000)]
                record[kind+'_master']={'number_of_nodes':int(len(entries)),
                    'first_temperature_K':float(entries[0,0]),
                    'last_temperature_K':float(entries[-1,0]),
                    'max_tabulated_fit_residual_fraction_full':float(np.max(entries[:,3])),
                    'max_tabulated_fit_residual_fraction_1000_100000K':float(np.max(valid[:,3]))}
            rr_entries=np.array(originals['RR'][key])
            dr_entries=np.array(originals['DR'][key])
            common=np.intersect1d(rr_entries[:,0],dr_entries[:,0])
            sampled=common[(common>=1000)&(common<=100000)]
            rr_common=rr_entries[np.isin(rr_entries[:,0],sampled)]
            dr_common=dr_entries[np.isin(dr_entries[:,0],sampled)]
            assert np.array_equal(rr_common[:,0],dr_common[:,0])
            total_original=rr_common[:,1]+dr_common[:,1]
            total_fit=rr_common[:,2]+dr_common[:,2]
            combined=np.abs(total_fit/total_original-1)
            record['source_printed_total_fit_comparison']={
                'max_absolute_relative_difference':float(np.max(combined)),
                'temperature_K':float(sampled[np.argmax(combined)]),
                'number_of_common_nodes_1000_100000K':int(len(sampled)),
                'qualification':'Computed from rounded printed original/fit RR and DR columns; not a theoretical error estimate'}
        rows.append(record)
    result={'case':'terminal outgoing electron-impact ionization and reverse RR+DR source audit',
            'rate_evaluation':'diagnostic source_terminal and recombination functions; no production rate activation',
            'master_table_interpretation':'DIFF is the source-native fit residual against original tabulated RR or DR values, not a theoretical uncertainty',
            'source_sha256':{str(p):digest(p) for p in [Path(sources.BADNELL_RR),Path(sources.BADNELL_DR),
                                Path(diagnostic.__file__),Path(sources.__file__),*MASTER.values()]},
            'channels':rows}
    for axis in axes.flat[len(diagnostic.ROSTER):]:
        axis.axis('off')
    handles,labels=axes.flat[0].get_legend_handles_labels()
    axes.flat[-1].legend(handles,labels,loc='center')
    figure.tight_layout()
    figure.savefig(HERE/'reverse_rate_profiles.pdf')
    plt.close(figure)
    np.savez_compressed(HERE/'full_profiles.npz',**profiles)
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    for x in rows:
        print(x['element'],x['source_key_Z_N_M'],x['rr_source'],x['dr_source'],
              'RR100k',x['RR_cm3_s'][-1],'DR100k',x['DR_cm3_s'][-1])


if __name__=='__main__':
    main()
