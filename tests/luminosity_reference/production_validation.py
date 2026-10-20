"""Check namelist validation, saved normalization, and default gas transport."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import h5py
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--work-dir', type=Path, required=True)
    ap.add_argument('--baseline', type=Path)
    args = ap.parse_args()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    density = work/'density.fits'
    fits.writeto(density, np.full((3,3,3), 100., dtype=np.float32))
    common = f"""&parameters
 par%no_photons=256
 par%iseed=1234
 par%grid_type='car'
 par%density_file='{density}'
 par%xmax=1
 par%ymax=1
 par%zmax=1
 par%distance_unit='pc'
 par%rmin=0
 par%gas_niter=1
 par%solve_te=.false.
 par%te_fixed=10000
 par%xHI_init=0.01
 par%abund_He=0.1
 par%nnu_ion=8
 par%tstar=40000
 par%luminosity=1.e36
 par%add_fuv=.true.
 par%nnu_fuv=4
 par%data_dir='{ROOT/'data'}'
 par%atomic_dir='{ROOT/'data/atomic'}'
 par%out_file='result.h5'
 par%file_format='hdf5'
 par%emis_output=.true.
"""
    env = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               OPENBLAS_NUM_THREADS='1', I_MPI_FABRICS='shm')
    records = {}
    def run(name, extra='', expected_error=None, binary=ROOT/'MoCHII.x'):
        directory = work/name
        directory.mkdir()
        inp = directory/'input.in'
        inp.write_text(common+extra+'\n/\n')
        command = ['mpirun','-np','2',str(binary),str(inp)]
        with (directory/'run.log').open('w') as stream:
            completed = subprocess.run(command,cwd=directory,env=env,stdout=stream,
                                       stderr=subprocess.STDOUT,timeout=180)
        text = (directory/'run.log').read_text()
        records[name] = {'return_code': completed.returncode, 'command': command,
                         'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest()}
        if expected_error:
            assert completed.returncode != 0 and expected_error in text, text[-3000:]
            return
        assert completed.returncode == 0 and '>>> STOP' in text, text[-3000:]
        with h5py.File(directory/'result_rates.h5') as f:
            assert all(np.isfinite(f[k]['data'][:]).all() for k in f if 'data' in f[k])
        return directory/'result_rates.h5'
    run('invalid'," par%luminosity_reference='invalid'",'par%luminosity_reference must be')
    run('slab'," par%luminosity_reference='full_spectrum'\n par%source_geometry='slab'\n par%xy_periodic=.true.\n par%slab_top_source=1",'requires internal point sources')
    run('external'," par%luminosity_reference='full_spectrum'\n par%nsource=0\n par%ext_intensity=1",'requires internal point sources')
    full = run('full'," par%luminosity_reference='full_spectrum'")
    spectrum = work/'absolute_spectrum.dat'
    np.savetxt(spectrum, [[1.,2.e34],[200.,2.e34]])
    absolute = run('absolute_unset',f" par%luminosity_reference='full_spectrum'\n"
                   f" par%spectrum_type='per_ev'\n par%ion_spectrum='{spectrum}'\n"
                   " par%luminosity=-999")
    with h5py.File(absolute) as f:
        np.testing.assert_allclose(f['SourceLuminosity/data'][0,2],199*2.e34,rtol=2e-14)
        np.testing.assert_allclose(f['Gamma_HI'].attrs['TOT_LUM'],199*2.e34,rtol=2e-14)
    run('mixed_unset',f" par%luminosity_reference='full_spectrum'\n"
        f" par%spectrum_type='per_ev'\n par%ion_spectrum='{spectrum}'\n"
        " par%nsource=2\n par%src_tstar(1)=50000", 'mixed absolute-file and Planck')
    columns = work/'absolute_columns.dat'
    np.savetxt(columns, [[1.,2.e34,4.e34],[200.,2.e34,4.e34]])
    multi = run('multi_absolute',f" par%luminosity_reference='full_spectrum'\n"
                f" par%spectrum_type='per_ev'\n par%src_spectrum_file='{columns}'\n"
                " par%nsource=2")
    with h5py.File(multi) as f:
        np.testing.assert_allclose(f['SourceLuminosity/data'][:,2],[199*2.e34,199*4.e34],rtol=2e-14)
        np.testing.assert_allclose(f['Gamma_HI'].attrs['TOT_LUM'],199*6.e34,rtol=2e-14)
    default = run('default')
    explicit = run('explicit'," par%luminosity_reference='ionizing_band'")
    with h5py.File(full) as f:
        values = f['SourceLuminosity/data'][:]
        assert values.shape == (1,7), values.shape
        v=values[0]
        assert v[0]==0 and v[1]==-1
        assert v[2]==1e36
        np.testing.assert_allclose(v[3]+v[4],v[6],rtol=2e-14)
        np.testing.assert_allclose(v[5]+v[6],v[2],rtol=2e-14)
        ref=np.asarray(f['Gamma_HI'].attrs['LUMREF']).ravel()[0]
        assert (ref.decode() if isinstance(ref,bytes) else ref)=='full_spectrum'
        records['full']['source_luminosity']=v.tolist()
    comparisons=[('default_vs_explicit',default,explicit)]
    if args.baseline:
        baseline=run('baseline',binary=args.baseline.resolve())
        comparisons.append(('baseline_vs_default',baseline,default))
    for label,a,b in comparisons:
        with h5py.File(a) as fa,h5py.File(b) as fb:
            keys=[k for k in fa if 'data' in fa[k]]
            for k in keys:
                np.testing.assert_array_equal(fa[k]['data'][:],fb[k]['data'][:],err_msg=k)
        records[label]={'all_numeric_datasets_identical':True,'datasets':keys}
    (work/'result.json').write_text(json.dumps({'passed':True,'records':records},indent=2)+'\n')
    print(json.dumps({'passed':True,'result':str(work/'result.json')}))

if __name__=='__main__':
    main()
