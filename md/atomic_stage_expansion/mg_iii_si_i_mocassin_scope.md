# Mg III and Si I line-cooling scope in the local MOCASSIN tree

Date: October 10, 2026. The MOCASSIN audit is read-only. The original
Mg III subset assessment below is historical; MoCHII subsequently installed
the corrected full 167-state atom and its cooling fit. External MOCASSIN
source and data were not changed.

The local MOCASSIN `data/fileNames.dat` lists `data/mgiii.dat` and
`data/si1.dat` in their element-stage positions. Those two numerical files
are absent from the checked source data directory and installed
`/home/kiseon/share/mocassin/data/` directory, although the two `fileNames.dat`
copies have identical SHA-256
`eea996e8546d3400d475be062eeacbb9e76f3ac7ef7131951438d599ad0df80b`.
Neither source release archive contains these two files. The `/opt` and
`/usr/local` installed data routes checked also lack them; the actual
compiled executable's `PREFIX` could not be established from those paths
alone. Therefore the local data do not establish any historical Mg III or
Si I numerical atom to reproduce.

The source invocation is explicit. `hydro_mod::makeElements` reads each
entry in `fileNames.dat`, allocates and reads its atomic arrays only if
`PREFIX/share/mocassin/<entry>` exists, then sets `lgDataAvailable(elem,ion)`
by opening the same file. `emission_mod::forLines` calls `equilibrium` only
when that flag is true. `update_mod` sums these `forbiddenLines` into the
heavy-element collisional cooling. A listed filename alone therefore does
not cause either ion to contribute line cooling. When a file exists,
`equilibrium` interpolates its tabulated effective collision strengths,
forms excitation/de-excitation rates with actual `Te`, includes radiative
decay, solves level populations, and sums `A_ul (E_u-E_l) n_u`. Its equations
also include a recombination source proportional to the next ion's density;
the standalone electron-only MoCHII cooling atom does not model that source.
The local `nForLevels=17` limit in `constants_mod` means a 27-level numerical
file would require a separate MOCASSIN limit change. The current installed
binary was not executed as part of this audit.

The first 27 bound Mg III CHIANTI levels form a useful **candidate source
subset**, independent of the absent MOCASSIN file: 2p^6 and the complete
2p^5 3s, 3p, 3d configurations. All 27 have observed `.elvlc` energies;
their 117 positive-A radiative rows have strictly positive observed gaps and
all 351 electron pairs have collision strengths. The `.wgfa` header assigns
the first-27 radiative rates to Jönsson et al. (2014), *Atomic Data and
Nuclear Data Tables* **100**, 1, DOI
[10.1016/j.adt.2013.06.001](https://doi.org/10.1016/j.adt.2013.06.001),
while `.scups` identifies the Liang and Badnell (2010) Ne-like R-matrix
calculation. This subset avoids the 16 high-level A orientations inverted
by the current mixed observed/theoretical energy vector. Its one
radiatively closed state is the ground. Its radiative rows and photon
energies are internally oriented, but the published accuracy of a
27-level truncation has not been established.

The 27-level subset omits 140 bound levels, 3,247 radiative rows, and
13,510 electron pairs from the current 167-level bound source. It also
omits direct ground excitation into those 140 higher levels. On 401
temperatures from the published Mg III source minimum 1,800 K to 100,000 K,
the direct ground-excitation cooling contribution retained by the first 27
falls from effectively 100% at the cold end to 82.80% at 100,000 K. The
omitted contribution first exceeds 0.3% near 31,824 K. These direct
ground-excitation sums do not validate full finite-density populations:
high-level cascades and collisional redistribution had not yet been compared at the time of that initial assessment.
The 17-level prefix permitted by the MOCASSIN source retains only 52
positive-A rows and 136 electron pairs, cutting through the 2p^5 3p
configuration rather than retaining the complete n=3 complex. Its direct
ground-excitation contribution is 73.28% of the full source at 100,000 K.
Neither smaller atom meets the 0.3% full-source cooling target over the
whole 1,800–100,000 K interval merely because it avoids high-level
orientation defects.

For Si I, the same absent-file and conditional-call result holds. The
[final Dufresne, Del Zanna, and Storey (2021) methods](https://academic.oup.com/mnras/article/505/3/3968/6283727)
state that Si I radiative data came from MCHF and ATSP2K calculations, while
its electron collision strengths came from Autostructure distorted waves
used to estimate metastable populations for ionization. The authors
explicitly caution against treating those neutral distorted-wave values
as reliable line-intensity predictions. This is independent of the anomalous
Si I 15→14 wavelength/gf/A row described in the cooling-source audit.

## Subsequent full Mg III implementation

The completed MoCHII model uses all 167 bound states, 3,364 spontaneous
rows, 13,861 collision pairs, and 207,915 original source nodes. The first
27-state radiation values are retained. Higher positive-gf electric-dipole
rates use physical level orientation and preserve line strength through
A proportional to nu^3/g_upper. Higher gf=0 rates remain explicitly
identified effective-rate approximations rather than assigned multipoles.
An observed energy incorrectly assigned to the J=0 source state 167 was
rejected after comparison with the official NIST level table; that state
uses its own theoretical energy of 634430 cm^-1. Other selected energies
are unchanged.

On 801 temperatures from 1,000 to 100,000 K at four electron densities,
the corrected full source and actual compiled level solver agree in cooling
to a maximum relative difference of 2.58e-14. Truncation to the first 27
states instead differs from the full atom by up to 17.19%, failing the
0.3% criterion. The installed cooling fit has maximum measured error
0.2107046%; the full registered table's sampled application error is
0.0005348%. These are numerical representation checks, not estimates of
atomic accuracy. Original source collision strengths are held constant
below their 1,800 K floor with actual-temperature Maxwellian factors.

Installation and source limitations are documented in
[the corrected-level validation](validation_2026-10-10/mgiii_corrected_level_identity/main_review.json)
and [the installation manifest](validation_2026-10-10/mgiii_corrected_level_identity/installation_manifest.json).
This Mg III calculation extends beyond the locally invoked MOCASSIN v2
line-cooling data, for which the Mg III numerical file is absent. It does
not establish whole-code parity or reproduce MOCASSIN's next-ion
recombination feeding.
