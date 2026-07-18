#!/usr/bin/env bash
# V3(e): diffuse-field variance study.  Same diffuse Stromgren (case A +
# diffuse, aligned grid = default), noise-limited photon budget: 3 ordinary-MC
# replicates (vary iseed, launch='random') vs 3 RQMC replicates (vary qmc_seed,
# launch='sobol' -> stellar AND diffuse on the two Sobol streams).  The checker
# reports the replicate scatter (std) of R_eff and V_ion for each set.
set -e
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
cd "$ROOT"
OUT=tests/qmc/v3sobol
BIN=./MoCHII.x
NP=8

gen () {  # $1=tag $2=mode(random|sobol) $3=seed
  cat > "$OUT/$1.in" <<EOF
&parameters
 par%no_photons    = 5.0e5
 par%no_print      = 5.0e5
 par%iseed         = $( [ "$2" = random ] && echo "$3" || echo 1234 )
 par%grid_type     = 'car'
 par%car_walk      = 'dda'
 par%nx            = 32
 par%ny            = 32
 par%nz            = 32
 par%xmax          = 2.0
 par%ymax          = 2.0
 par%zmax          = 2.0
 par%nH_const      = 200.0
 par%rmax          = 1.9
 par%dust_model    = 'none'
 par%distance_unit = 'pc'
 par%use_ion_band  = .true.
 par%nnu_ion       = 16
 par%eion_min      = 13.598
 par%eion_max      = 100.0
 par%tstar         = 4.0e4
 par%luminosity    = 1.5e37
 par%xHI_init      = 1.0
 par%xHeI_init     = 1.0
 par%xHeII_init    = 0.0
 par%He_abund      = 0.1
 par%ion_add_dust  = .false.
 par%gas_niter     = 40
 par%gas_tol       = 5.0e-3
 par%te_fixed      = 1.0e4
 par%ion_relax     = 1.0
 par%case_ab       = 'A'
 par%diffuse_field = .true.
 par%conv_crit     = 'vol'
 par%source_geometry = 'point'
 par%launch_sequence = '$2'
 par%qmc_seed        = $( [ "$2" = sobol ] && echo "$3" || echo 12345 )
 par%file_format     = 'hdf5'
 par%out_file        = '$1.h5'
/
EOF
}

for s in 11 22 33; do
  t="dv_mc_$s"; gen "$t" random "$s"
  mpirun -np $NP "$BIN" "$OUT/$t.in" > "$OUT/$t.log" 2>&1
  mv "${t}_rates.h5" "$OUT/${t}_rates.h5"; echo "  ran $t"
done
for s in 101 202 303; do
  t="dv_qmc_$s"; gen "$t" sobol "$s"
  mpirun -np $NP "$BIN" "$OUT/$t.in" > "$OUT/$t.log" 2>&1
  mv "${t}_rates.h5" "$OUT/${t}_rates.h5"; echo "  ran $t"
done
echo "V3(e) diffuse-variance runs complete."
