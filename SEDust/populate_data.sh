#!/bin/bash
# Populate the SEDust optics data this tree does not keep in git.
#
# Large generated optics files are excluded by .gitignore.  MoCHII's default
# photoionization path needs the EUV companion Q table; the base Q table is
# copied as well so the embedded SEDust tree remains complete:
#
#   tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400_euv.dat   (33 MB)
#     The random-orientation T-matrix Q table of the b/a = 1.4 astrodust
#     spheroid, 169 radii x 1762 wavelengths spanning 1.0e-4 to 3.981e4 um.
#     It is par%sed_qtable's default for MoCHII, and both the 'astrodust' and
#     'dl07' models read it.  ('zubko' is on its own DustEM tables under
#     data/zubko/, which are committed, and needs nothing from here.)
#
#   data/dielectric/qlib_gra_D16MGemt_1.400             (17 MB)
#     The D16 graphite spheroid Q-library, read only by q_graphite_d16_mod,
#     which nothing selects by default.
#
# Everything else the dust path reads is committed.  This script re-copies it
# as well, so it also serves as a full refresh from a canonical SEDust tree.
#
# >>> EDIT the path below to your own SEDust location, or run with
# >>>   SEDUST_SRC=/your/path/to/SEDust ./populate_data.sh
set -e
SED=${SEDUST_SRC:-/home/kiseon/MoCafe/Grain/SEDust_v1.00}
if [ ! -d "$SED" ]; then
  echo "ERROR: SEDust source tree not found at: $SED" >&2
  echo "       Set it to your own location, e.g.:" >&2
  echo "         SEDUST_SRC=/path/to/SEDust $0" >&2
  exit 1
fi
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HERE/data/dielectric" "$HERE/data/release" "$HERE/data/zubko" "$HERE/tmatrix/output"
# Dielectric functions and PAH cross sections.  index_DH21Ad_P0.20_0.00_1.400
# is the astrodust refractive index the Q table above was computed from; SEDust
# reads it when a run asks for wavelengths below the table's own short end.
cp "$SED"/data/dielectric/{DH21_aeff,index_CpaD03,index_CpeD03,index_DH21Ad_P0.20_0.00_1.400,index_silD03,PAHion.31,PAHneu.31,q_D16graphite.dat,qlib_gra_D16MGemt_1.400} "$HERE/data/dielectric/"
# Size-integrated extinction curves dust_extinction serves to the transport,
# one for each par%dust_model (calc_kext.x writes them in the SEDust tree).
# The three below are the tables SEDust falls to when par%sed_kext is blank.
cp "$SED"/data/{kext_astrodust_MW_euv.dat,kext_dl07_MW_euv.dat,kext_zubko_BARE_GR_S.dat} "$HERE/data/"
# The plain astrodust curve, written by `calc_kext.x astrodust`: the same 1762
# wavelengths as the _euv one, taking the whole range straight off the Q table
# instead of recomputing below 0.0912 um from the dielectric function.  The two
# agree to 5e-8 on C_ext, C_abs and C_sca; nothing selects this one on its own,
# and it is here for a run that names it through par%sed_kext.
cp "$SED"/data/kext_astrodust_MW.dat "$HERE/data/"
cp "$SED"/data/kext_dl07_MW.dat "$HERE/data/"
cp "$SED"/data/release/{extinction.dat,kext_albedo_WD_MW_3.1_60_D03.all_2003,scattering.dat,size_distribution.dat} "$HERE/data/release/"
cp "$SED"/tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat "$HERE/tmatrix/output/"
cp "$SED"/tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400_euv.dat "$HERE/tmatrix/output/"
# Zubko (ZDA 2004 BARE-GR-S) optics/calorimetry/config for par%dust_model='zubko'
cp "$SED"/data/zubko/* "$HERE/data/zubko/"
echo "Populated SEDust data under $HERE (from $SED)"
