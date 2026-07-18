#!/usr/bin/env bash
# V3 variance study: MC (vary iseed) vs RQMC (vary qmc_seed) on the fixed-state
# G0 attenuation gate.  ion_align_edges off so the analytic log-grid reference
# (check_v3.py, reused from the g0 checker) is exact.  Outputs to tests/qmc/v3/.
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"
OUT=tests/qmc/v3
mkdir -p "$OUT"
FITS="$ROOT/tests/g0_gamma/amr_sphere_L5.fits"
NP=8

gen_input () {   # $1=tag $2=nphot $3=mode(random|sobol) $4=seedval
  cat > "$OUT/$1.in" <<EOF
&parameters
 par%no_photons    = $2
 par%no_print      = $2
 par%iseed         = $( [ "$3" = random ] && echo "$4" || echo 1234 )
 par%grid_type     = 'amr'
 par%amr_file      = '$FITS'
 par%dust_model    = 'none'
 par%distance_unit = 'pc'
 par%use_ion_band  = .true.
 par%nnu_ion       = 16
 par%ion_align_edges = .false.
 par%eion_min      = 13.598
 par%eion_max      = 100.0
 par%tstar         = 1.0e5
 par%luminosity    = 1.0e38
 par%xHI_init      = 0.1
 par%xHeI_init     = 1.0
 par%xHeII_init    = 0.0
 par%He_abund      = 0.1
 par%ion_add_dust  = .false.
 par%source_geometry = 'point'
 par%launch_sequence = '$3'
 par%qmc_seed        = $( [ "$3" = sobol ] && echo "$4" || echo 12345 )
 par%file_format     = 'hdf5'
 par%out_file        = '$1.h5'
/
EOF
}

MC_SEEDS="11 22 33 44"
QMC_SEEDS="101 202 303 404"
for NPHOT_TAG in "131072 p17" "1048576 p20"; do
  set -- $NPHOT_TAG; NPHOT=$1; PT=$2
  for s in $MC_SEEDS; do
    tag="mc_${PT}_${s}"; gen_input "$tag" "$NPHOT" random "$s"
    mpirun -np $NP ./MoCHII.x "$OUT/$tag.in" > "$OUT/$tag.log" 2>&1
    mv "${tag}_rates.h5" "$OUT/${tag}_rates.h5"
  done
  for s in $QMC_SEEDS; do
    tag="qmc_${PT}_${s}"; gen_input "$tag" "$NPHOT" sobol "$s"
    mpirun -np $NP ./MoCHII.x "$OUT/$tag.in" > "$OUT/$tag.log" 2>&1
    mv "${tag}_rates.h5" "$OUT/${tag}_rates.h5"
  done
done
echo "V3 runs complete -> $OUT"
