"""Compare the actual cooling interpolator with SciPy PCHIP and zero limits."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.interpolate import PchipInterpolator
from run_checks import compile_private_module, link

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    work=parser.parse_args().work_dir.resolve()
    original=ROOT/'src/nlevel_cooling_mod.f90'
    marker='  public :: nlevel_cooling_lam_se, nlevel_cooling_report\n'
    obj=compile_private_module(original,work,[(marker,marker+'  public :: reg,col_cooling\n')])
    exe=link(work,'tests/ion_stages/check_cooling_temperature_interpolation.f90',[obj])
    temperatures=np.array([100.,1000.,3000.,10000.,32000.])
    cases=[('boltzmann',temperatures,np.exp(-10000/temperatures)),
           ('constant',temperatures,np.full(5,1e-30)),
           ('thermal_peak',temperatures,np.array([1e-30,1e-22,3e-21,2e-21,1e-22])),
           ('cold_positive',temperatures,np.array([1e-250,1e-220,1e-200,1e-190,1e-188])),
           ('two_temperatures',temperatures[:2],np.array([1e-25,1e-21])),
           ('zero_endpoints',temperatures,np.array([0.,1e-25,1e-23,1e-24,0.])),
           ('isolated_positive_pair',temperatures,np.array([0.,0.,1e-23,1e-24,0.])),
           ('zero_cooling',temperatures,np.zeros(5))]
    lines=[str(len(cases))];expected=[];kinds=[];bounds=[]
    for index,(name,t,power) in enumerate(cases,1):
        query=np.unique(np.concatenate((t,np.geomspace(t[0],t[-1],97))))
        query=np.clip(query,t[0],t[-1])
        lines.append(f'{len(t)} {len(query)}')
        lines.extend(f'{a:.17g} {b:.17g}' for a,b in zip(t,power))
        lines.extend(f'{a:.17g}' for a in query)
        for q in query:
            i=min(np.searchsorted(t,q,side='right')-1,len(t)-2)
            i=max(i,0);positive=(power[i]>0 and power[i+1]>0)
            if positive:
                left=i;right=i+1
                while left>0 and power[left-1]>0:left-=1
                while right<len(t)-1 and power[right+1]>0:right+=1
                value=np.exp(PchipInterpolator(-1/t[left:right+1],np.log(power[left:right+1]))(-1/q))
            else:
                weight=(np.log10(q)-np.log10(t[i]))/(np.log10(t[i+1])-np.log10(t[i]))
                value=(1-weight)*power[i]+weight*power[i+1]
            expected.append(value);kinds.append(name);bounds.append((min(power[i:i+2]),max(power[i:i+2])))
    inp=work/'profiles.in';inp.write_text('\n'.join(lines)+'\n');output=work/'profiles.txt'
    run=subprocess.run([str(exe),str(inp),str(output)],cwd=ROOT,capture_output=True,text=True,timeout=60)
    (work/'native.log').write_text(run.stdout+run.stderr)
    if run.returncode:raise RuntimeError('Native interpolation failed')
    measured=np.loadtxt(output)[:,2];expected=np.array(expected)
    error=abs(measured-expected)/np.maximum(expected,1e-300)
    if len(measured)!=len(expected) or max(error)>1e-10:
        raise ValueError('Native PCHIP or zero limit differs from independent reference')
    for value,(lo,hi) in zip(measured,bounds):
        if value<lo*(1-1e-10) or value>hi*(1+1e-10):raise ValueError('Endpoint bound exceeded')
    result=dict(passed=True,samples=len(expected),profiles=len(cases),
        maximum_relative_reference_difference=float(max(error)))
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS:',len(expected),'actual cooling interpolation samples')


if __name__=='__main__':main()
