#!/bin/bash
# Populate SEDust optics data that is NOT committed to git.
#
# The minimal data set the default 'astrodust'+PAH dust-emission path needs is
# already in the repository (see .gitignore), so a fresh checkout runs
# dust emission out of the box -- you normally do NOT need this script.
#
# The only file kept out of git is the 17 MB D16 spheroid Q-library
# (qlib_gra_D16MGemt_1.400), read only by q_graphite_d16_mod, which nothing uses
# by default.  Run this script to fetch it if you enable that non-default
# graphite path.  It also re-copies the rest, so it doubles as a full refresh
# from a canonical SEDust tree.
#
# >>> EDIT the path below to your own SEDust location, or run with
# >>>   SEDUST_SRC=/your/path/to/SEDust ./populate_data.sh
set -e
SED=${SEDUST_SRC:-/home/kiseon/MoCafe/Grain/SEDust}
if [ ! -d "$SED" ]; then
  echo "ERROR: SEDust source tree not found at: $SED" >&2
  echo "       Set it to your own location, e.g.:" >&2
  echo "         SEDUST_SRC=/path/to/SEDust $0" >&2
  exit 1
fi
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HERE/data/dielectric" "$HERE/data/release" "$HERE/data/zubko" "$HERE/tmatrix/output"
cp "$SED"/data/dielectric/{DH21_aeff,index_CpaD03,index_CpeD03,index_silD03,PAHion.31,PAHneu.31,q_D16graphite.dat,qlib_gra_D16MGemt_1.400} "$HERE/data/dielectric/"
cp "$SED"/data/kext_astrodust_MW.dat "$HERE/data/"
cp "$SED"/data/release/{extinction.dat,kext_albedo_WD_MW_3.1_60_D03.all_2003,scattering.dat,size_distribution.dat} "$HERE/data/release/"
cp "$SED"/tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat "$HERE/tmatrix/output/"
# Zubko (ZDA 2004 BARE-GR-S) optics/calorimetry/config for dust_model_sed='zubko'
cp "$SED"/data/zubko/* "$HERE/data/zubko/"
echo "Populated SEDust data under $HERE (from $SED)"
