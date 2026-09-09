#!/bin/bash
# Refresh this tree's SEDust data from a canonical SEDust tree.
#
# You normally do NOT need this script: everything MoCHII opens at run time is
# kept in git, including each model's optics product (.gitignore un-ignores
# SEDust/data/*/*.h5), so a fresh checkout runs the dust path as it stands.
# Run it to bring the DATA of this copy up to date with a newer SEDust, or to
# fetch the one file deliberately kept out of git -- the 17 MB D16 spheroid
# Q-library data/dielectric/qlib_gra_D16MGemt_1.400, read only by
# q_graphite_d16_mod, which no par%dust_model selects.  It copies data and
# nothing else: the library sources under sed/src/ are a separate step, done by
# hand against the same tree, and an optics product that came out of a newer
# SEDust than the sources here is exactly the pairing to avoid.  See the
# paragraph on sed/src/ below for what that step covers.
#
# WHAT THE LAYOUT IS.  SEDust keeps one directory per dust model, holding what
# that model owns -- its optics product sedust_<model>.h5 and, where the model
# IS a set of files, the files that define it -- and keeps what the models SHARE
# beside those directories:
#
#   data/astrodust/  data/dl07/  data/mrn/  data/zubko/   one per par%dust_model
#   data/dielectric/                             optical constants, PAH cross sections
#   data/release/                                published tables + the size distribution
#
# WHAT IS COPIED, AND WHAT IS NOT.  sedust_<model>.h5 is the form the optics
# take here: it carries the model's wavelength axis, its cross-section tables
# and its size-integrated extinction curve (/kext) together, and MoCHII reads
# the optics from it alone -- par%sed_data_dir names the directory, build_dust
# resolves the file.  SEDust also writes the same optics as text (the Q tables
# q_*.dat and the extinction curves kext_*.dat), and those are NOT copied: they
# are ~250 MB of a second copy of what the .h5 already holds, and nothing in a
# MoCHII run opens one.  Every other file in a model directory is an INPUT
# rather than a product -- Zubko's ZDA config, size distributions, DustEM optics
# and calorimetry; astrodust's DH21 size/wavelength tabulations -- so the rule
# below is "the directory, minus SEDust's own text products".  MRN's directory
# comes out of that rule with the product alone, which is right: that model is
# a power law and two D03 dielectric functions, both of them stated in code and
# in dielectric/, so it owns no input file of its own.
#
# The sources under sed/src/ are byte-identical to their SEDust v1.00
# counterparts, so `cmp` against that tree is the whole check that this copy is
# current.  sed/build_lib.sh is the one exception: it differs from upstream's by
# a single entry, dustem_io, which SEDust v1.00 added as a library module (see
# `use dustem_io` in sed_astrodust.f90) without adding it to its own build_lib.sh
# source list.  The one file of the library that reaches into the
# T-matrix, src/euv_astrodust_tmatrix.f90, is deliberately not copied: it serves
# the astrodust EUV band, which only a lam_min shorter than the model's own axis
# asks for, and MoCHII passes no lam_min.  Without it the archive links with no
# T-matrix at all, and euv_tmatrix = .true. is refused (build_dust status 8;
# build_lib.sh's own message quotes build_astrodust's 6) instead of being
# answered with a different grain.
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

# Shared material data.  A dielectric function belongs to no one model -- DL07
# and Zubko read the same D03 astrosilicate -- so it sits beside the model
# directories rather than inside one of them.  build_dust resolves them from the
# data directory it is given, the same one the model directories come from, so
# they have to be here and nowhere else.
mkdir -p "$HERE/data/dielectric" "$HERE/data/release"
cp "$SED"/data/dielectric/{index_CpaD03,index_CpeD03,index_silD03,index_DH21Ad_P0.20_0.00_1.400,PAHion.31,PAHneu.31,q_D16graphite.dat,qlib_gra_D16MGemt_1.400} "$HERE/data/dielectric/"
cp "$SED"/data/release/{extinction.dat,kext_albedo_WD_MW_3.1_60_D03.all_2003,scattering.dat,size_distribution.dat} "$HERE/data/release/"

# One directory per model, minus SEDust's own text optics products.  The two
# name patterns are what SEDust writes and what the .h5 already contains:
# q_*.dat from make_qtable.x, kext_*.dat from calc_kext.x.  Neither pattern
# matches an input -- Zubko's DustEM tables are Gra_/suvSil_/PAH_*.dat and
# astrodust's DH21 refractive-index tabulation is q_DH21Ad_*.dat.gz.
#
# SEDust v1.00 also builds two further models, themis and g18d, and this loop
# deliberately leaves them out.  Their wavelength axis runs 0.04 to 1e5 um (800
# points), so the shortest wavelength they carry is 31.0 eV, while MoCHII
# transports 13.598-100 eV.  A model that cannot cover the transported band
# would be refused by the band-coverage check in gas_opacity_mod anyway, so
# offering it here would only move the failure later.
for m in astrodust dl07 mrn zubko; do
  mkdir -p "$HERE/data/$m"
  find "$SED/data/$m" -maxdepth 1 -type f \
       ! -name 'q_*.dat' ! -name 'kext_*.dat' \
       -exec cp {} "$HERE/data/$m/" \;
done

echo "Populated SEDust data under $HERE (from $SED)"
