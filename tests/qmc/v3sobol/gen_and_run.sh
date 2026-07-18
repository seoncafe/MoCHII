#!/usr/bin/env bash
# V3: physics gates under launch_sequence='sobol'.  Builds sobol variants of the
# recorded gate inputs and runs them with the production binary.  Rates files
# land in tests/qmc/v3sobol/<tag>_rates.h5 for check_v3sobol.py.
set -e
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
cd "$ROOT"
OUT=tests/qmc/v3sobol
BIN=./MoCHII.x
NP=${NP:-8}

# make a sobol variant of an input: strip the final '/', append launch lines,
# repoint out_file to <tag>, close the namelist.
mk () {  # $1=src_input  $2=tag  $3=qmc_seed
  local src="$1" tag="$2" seed="$3"
  sed -e "s#par%out_file.*#par%out_file      = '${tag}.h5'#" \
      -e '/^\s*\/\s*$/d' "$src" > "$OUT/${tag}.in"
  cat >> "$OUT/${tag}.in" <<EOF
 par%launch_sequence = 'sobol'
 par%qmc_seed        = ${seed}
/
EOF
}

run () {  # $1=tag  $2=np
  mpirun -np "$2" "$BIN" "$OUT/$1.in" > "$OUT/$1.log" 2>&1
  mv "$1_rates.h5" "$OUT/$1_rates.h5"
  echo "  ran $1 (np=$2)"
}

# (a) two-point superposition ; (b) ext_rec + ext_sph ; (c) mixed
mk tests/multi_src/two_point.in  v3s_twop    101
mk tests/ext_field/ext_rec.in    v3s_extrec  101
mk tests/ext_field/ext_sph.in    v3s_extsph  101
mk tests/multi_src/mixed.in      v3s_mixed   101
run v3s_twop   "$NP"
run v3s_extrec "$NP"
run v3s_extsph "$NP"
run v3s_mixed  "$NP"

# (d) MPI independence: multi-source AND external sobol at np=4 vs np=8
mk tests/multi_src/two_point.in  v3s_twop_np4   101 ; run v3s_twop_np4   4
mk tests/multi_src/two_point.in  v3s_twop_np8   101 ; run v3s_twop_np8   8
mk tests/ext_field/ext_rec.in    v3s_extrec_np4 101 ; run v3s_extrec_np4 4
mk tests/ext_field/ext_rec.in    v3s_extrec_np8 101 ; run v3s_extrec_np8 8
echo "V3 sobol runs complete."
