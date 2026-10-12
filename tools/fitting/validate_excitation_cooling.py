"""Check selected excitation cooling coefficients through native evaluation.

Below the normal floating-point range, source accuracy is checked separately
in logarithmic space. Numerical underflow is permitted without a rate floor.
"""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
from atomic_output import atomic_write_text
from validate_added_charge_transfer import compile_physics_check
from validate_stage_cooling import coefficients
from fit_cooling import evaluate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--products',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--ion',action='append',metavar='NAME:ELEMENT:STAGE',
                        help='Validate specified staged ions; repeat as needed')
    args=parser.parse_args()
    work=args.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
    exe=compile_physics_check('tests/ion_stages/check_cooling_fit.f90',work)
    result={'passed':True,'build_mode':'debug','ions':{}}
    selected=[('CV','c',5),('NVI','n',6),('SVII','s',7),('NeIV','ne',4)]
    if args.ion:
        selected=[]
        for item in args.ion:
            name,el,stage=item.split(':')
            selected.append((name,el,int(stage)))
    for name,el,stage in selected:
        path=args.products.resolve()/name/f'cooling_{el}_{stage}.txt'
        output=work/(name+'_native.txt')
        subprocess.run([str(exe),str(path),str(output)],check=True)
        sample=np.loadtxt(output)
        if sample.shape!=(801,2):raise ValueError('Incomplete native profile')
        a,t=coefficients(path);expected=evaluate(a,t,sample[:,0])
        normal=expected>1e12*np.finfo(float).tiny
        error=float(np.max(np.abs(sample[normal,1]/expected[normal]-1)))
        # This bound validates floating-point behavior; it is not a physical
        # cooling cutoff and never changes a computed coefficient.
        small=~normal
        if (error>1e-10 or np.any(sample[:,1]<0) or
            not np.all(np.isfinite(sample[:,1])) or
            np.any(np.abs(sample[:,1]-expected) >
                   1e-10*np.abs(expected)+len(a)*np.finfo(float).tiny)):
            raise ValueError(f'{name}: native cooling error {error}')
        result['ions'][name]=dict(points=801,normal_relative_check_points=int(normal.sum()),
                                  maximum_relative_difference=error,
                                  small_values_check='Absolute roundoff bound n_terms * minimum normal value, plus relative tolerance',
                                  native_zero_points=int(np.sum(sample[:,1]==0)),
                                  source_accuracy='See dense-grid logarithmic source audit in ion JSON')
    atomic_write_text(work/'native_result.json',json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
