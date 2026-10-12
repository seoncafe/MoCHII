"""Validate final-journal Fe fits in the linked equilibrium calculations.

CI and charge transfer are zero only in the isolated coefficient fixture.
The seven-stage fixture does not authorize a physical seven-stage Fe roster.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import make_element_data as G
from validate_added_charge_transfer import compile_physics_check
from zhang_fe_recombination import DR_FITS, RR_FITS, dr_fit, rr_fit, compare, ARCHIVE

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step',choices=['prepare','fortran'])
    parser.add_argument('--work-dir',required=True,type=Path)
    args=parser.parse_args()
    work=args.work_dir.resolve(); work.mkdir(parents=True,exist_ok=True)
    if args.step=='prepare':
        atomic=work/'atomic'; atomic.mkdir()
        rows=['# Isolated RR+DR coefficient fixture; CI and CX disabled.',
              'ELEMENT fe 26 7']
        for q in range(1,7):
            rows += [f'TRANSITION {q}',f'ETH {10*q}',
                     f'PHOTO {10*q} 10 1000 1 1 0 0 1',
                     'PHINNER 1000',f'PHSHELL 0 {10*q} 10 1 1 2 1',
                     'CI 1 0 0 1 0']
            rr=RR_FITS.get(q,(1e-11,.5,1e4,1e4,0,1))
            rows.append('RR '+' '.join(f'{x:.6e}' for x in rr))
            c,e=DR_FITS.get(q,((),()))
            rows.append('DR '+str(len(c))+' '+
                        ' '.join(f'{x:.6e}' for x in (*c,*e)))
            rows += ['CXI 0','CXR 0','CXHEAT 0 0 0']
        (atomic/'element_fe.txt').write_text('\n'.join(rows)+'\n')
        # Exercise the actual generator, including q=3 correction and q=4 addition.
        G.generate([('fe',26,5)],work/'generated')
        audit=compare(ARCHIVE)
        (work/'adf09_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
        result={'passed':True,'fixture_stages':7,'generated_stages':5,
                'fixture_sha256':hashlib.sha256((atomic/'element_fe.txt').read_bytes()).hexdigest()}
    else:
        exe=compile_physics_check('tests/ion_stages/check_fe_recombination.f90',work)
        subprocess.run([str(exe),str(work/'atomic'),str(work/'native_rates.txt')],
                       cwd=ROOT,check=True)
        errors=[]; counts={q:0 for q in DR_FITS}
        for line in (work/'native_rates.txt').read_text().splitlines():
            q,t,a=line.split(); q=int(q); t=float(t); a=float(a)
            expected=rr_fit(q,t)+dr_fit(q,t)
            errors.append(abs(a/expected-1)); counts[q]+=1
        assert all(n==301 for n in counts.values()), counts
        assert max(errors)<1e-12, max(errors)
        result={'passed':True,'samples':len(errors),'maximum_relative_error':max(errors),
                'source':'Published ground-level S1 DR plus RR, not suppressed at finite density',
                'domain_K':'20*q^2 through 1e7*q^2 for each q=3,4,5,6',
                'limitation':'Compact DR accuracy relative to ADF09 exceeds 2% at some nodes'}
    (work/(args.step+'_result.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':
    main()
