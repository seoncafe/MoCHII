# Cooling-file gaps in the registered metal stages

Date: October 10, 2026. This audit follows the current `species_setup` and
`metal_cooling` code path, not a prior stage list. For every registered stage,
`species_setup` constructs `data/atomic/cooling_<element>_<stage>.txt` and sets
`has_cool` only if that file exists. `metal_cooling` sums `cooling_eval` only
for stages with `has_cool`. A missing file therefore contributes no line
cooling through that routine; an existing excitation file alone is not read
there. These are source and numerical-fit findings, not permission to install
candidate coefficients.

| Stage without a production cooling file | CHIANTI v11 `.elvlc/.wgfa/.scups` | Bound, collision-connected levels / radiation rows / electron rows | Radiatively closed states | Source assessment |
| --- | --- | ---: | --- | --- |
| C IV | All three exist | 24 / 74 / 66 | 1 | The `.scups` header gives a 3,200–32,000,000 K collision-strength range; a 1,000 K fit would extrapolate below it. |
| Ne I | All three exist | 31 / 142 / 465 | 1, 15 | Source-connected low-density atom; two closed states require stationary collisional balance. Numerical fit staged separately. |
| Mg I | All three exist | 15 / 45 / 96 | 1, 2 | The source states that LS R-matrix strengths were split into fine structure with distorted-wave ratios and cautions against spectroscopic diagnostics. Two type-2 ordinate arrays have a negative *infinite-temperature endpoint*; direct Upsilon on the finite 1,000–100,000 K grid remained positive. The present source-node atom writer rejects those negative nodes. The numerical cooling fit does not certify source accuracy. |
| Mg III | All three exist | 167 / 3,364 / 13,861 | 1 | Current selected energies reverse 16 positive-A source orientations and 122 electron source index orders. The Python low-density calculation already reorients symmetric electron pairs; one closed state makes total cooling independent of higher-state radiative branching. Source A assignments remain unresolved for a native excitation atom. The source article's abstract gives a 1,800 K minimum. |
| Si I | All three exist | 35 / 229 / 595 | 1 | The `.wgfa` row `15 14` has a positive A although observed E15 exceeds E14 by 67.493 cm^-1; its 1,481,635.13 Å wavelength matches the *magnitude* of this gap. Wavelength alone does not establish decay direction or reconcile its listed A and gf. Source provenance needs review before fitting. |
| Si V | All three exist | 169 / 3,778 / 14,196 | 1 | Current selected energies reverse 15 electron source index orders and make 11 pairs equal in energy. The source-node writer now retains these symmetric electron pairs without invented splitting, and the private full-atom native comparison passes. The source article's abstract gives a 5,000 K minimum. |
| Cl I | None of the three | Not applicable | Not applicable | Cloudy c25.00 Stout `cl_1.coll` has one `RATE ELECTRON 1 2` row, nine identical 4.7e-8 cm^3/s entries at 30–3,000 K, attributed by the file to Hollenbach and McKee (1989). Its `.nrg` has a 0/882.351 cm^-1 ground fine-structure pair with g=4/2, and `.tp` lists A(2→1)=0.0123 s^-1. The rate coefficient's direction and temperature convention need verification before deriving cooling; the source stops at 3,000 K. |
| Cl V | None of the three | Not applicable | Not applicable | Cloudy Stout `cl_5.coll` has no electron row. |
| Ca I | None of the three | Not applicable | Not applicable | Cloudy Stout `ca_1.coll` says only `Use gbar`; no tabulated electron strength or rate was found there. |

The counts and reversal tests above come from direct calls to
`chianti_cooling.bound_excitation_atom` and from its selected level energies.
The six CHIANTI atoms have no nonpositive selected radiative or electron
source-index gaps except the identified Mg III, Si I, and Si V cases. The
Python cooling reader already orders each symmetric electron pair by selected
energy and uses the absolute physical gap. An exactly degenerate pair has
zero excitation energy but can still exchange populations; its Burgess-Tully
scaling energy remains a separate positive source quantity. The native atom
writer retains such a pair with its zero physical gap and separate positive
Burgess-Tully scaling energy. A radiative A value
depends on its true upper state and cannot be reassigned by index without
source review. For Mg III, replacing all
selected energies with its theoretical `.elvlc` vector removes every negative
gap but leaves two *zero* collision gaps, 100–101 and 112–113, where the
theoretical energies are printed equal. A vector with observed Mg III levels
1–27 and theoretical energies above 27 removes all radiative inversions,
but it disagrees by more than 1% with 580 of the 1,081 positive source
wavelengths; the current mixed vector matches every positive wavelength
within 2.45e-6 relative. Thus absence of inversions alone is insufficient
to adopt that hybrid photon-energy vector. The 16 inverted rows instead have
negative theoretical wavelengths consistent with theoretical gaps and their
listed A and gf values. No direct use of the `.scups` scaling energy
as a physical splitting is justified. For Si V, a private hybrid with
observed levels 1–23 and theoretical levels above 23 has no reversed
radiative or electron gap. However, observed level 23 is above observed
level 24 whereas their theoretical order is opposite. Selecting theoretical
energy for level 24 solely to remove a reversal would not justify the
physical energy assignment. No row was reordered in the production source
files.

The 11 exact Si V zero-gap collision pairs occur in five groups:
39–40, 93–94, 96–97, 98–99, 102–103, and the six pairs among 112–115.
Their theoretical `.elvlc` energies have small positive splittings, but
these are not substitutes for the measured equal energies without an
independent physical identity and accuracy assessment. Retaining the zero
physical gaps is valid for collision balance; no photon energy is assigned
to a zero-gap exchange.

The Si I anomalous row connects `3s2 3p 3d 1D2` (level 15, g=5) and
`3s2 3p 4p 1P1` (level 14, g=3). It lists A=132.8 s^-1 and gf=0.03558.
The standard E1 gf-to-A conversion with the listed wavelength and g=5 gives
about 21.6 s^-1; other strong rows usually agree with their conversion to
within a few percent. Simply swapping the indices would therefore not
resolve this row's physical inconsistency. The energy, wavelength, gf, and
A provenance must be reconciled together. The `.wgfa` comments cite Tachiev
and Froese Fischer (2002), but that article concerns nitrogen-like and
oxygen-like sequences and does not establish Si I radiative rates. The
[final Dufresne, Del Zanna, and Storey (2021) methods](https://academic.oup.com/mnras/article/505/3/3968/6283727)
attribute the Tachiev addition to N I; for Si I they describe MCHF rates
from Fischer (2005), additional ATSP2K rates, and Autostructure distorted-wave
electron strengths. They caution that neutral distorted-wave strengths are
not usually reliable for line intensities and used these rates to estimate
metastable populations for ionization. The `.wgfa` file also states that
`wgfa_tidy` processed it in 2024, while `.elvlc` identifies NIST ASD
observed energies. Neither source identifies the calculation or
normalization of the anomalous row specifically.

The C IV source header attributes its data to Liang and Badnell (2011),
*A&A* **528**, A69. That paper's
[ADS abstract](https://ui.adsabs.harvard.edu/abs/2011A%26A...528A..69L/abstract)
prints the Li-like interval as `200(z+1)^2` to `2e6(z+1)^2` K, with residual
charge `z`; C IV has `z=3`, hence 3,200–32,000,000 K. The Mg III and Si V source headers attribute their data
to Liang and Badnell (2010), *A&A* **518**, A64. The latter paper's
[ADS abstract](https://ui.adsabs.harvard.edu/abs/2010A%26A...518A..64L/abstract)
prints the effective-strength interval as
`200(q+1)^2` to `2e6(q+1)^2` K, with residual charge `q`; thus the two
source intervals are 1,800–18,000,000 K for Mg III and
5,000–50,000,000 K for Si V. The publisher's final PDF was not
accessible during this audit, so detailed methods and accuracy claims remain
unverified. The Ne I source cites Zatsarinny and Bartschat (2012), *Physical
Review A* **86**, 022717; its
[ADS abstract](https://ui.adsabs.harvard.edu/abs/2012PhRvA..86b2717Z/abstract)
describes electron energies from excitation threshold to 300 eV, not a
Maxwellian-temperature accuracy interval. Mg I and Si I source comments do
not establish a 1,000–100,000 K physical accuracy interval either.

The private [existing_remaining_cooling](validation_2026-10-10/existing_remaining_cooling)
controller completed Ne I and Mg I with the current Python source reader and
an immutable compiled `cooling_load`/`cooling_eval` driver. On 801 temperatures
from 1,000 to 100,000 K, the largest source-fit differences were 0.03177%
and 0.02037%, respectively; the native reader agreed with serialized
coefficients within 8.16e-12 and 6.22e-12 at resolved points. Their full
cooling and error PDFs were inspected. The checked PDF, serialized
coefficients, source hashes, and native-reader output are numerical
diagnostics, not evidence of published physical accuracy over that entire
interval. Si I remains excluded pending the identified radiation check.
Mg III and Si V source-valid restricted-temperature low-density diagnostics
also completed on 801 temperatures: Mg III on 1,800–100,000 K has 0.01136%
maximum fit difference and 6.42e-12 maximum resolved native-reader
difference; Si V on 5,000–100,000 K has 0.37174% maximum fit difference
and 4.45e-16 maximum resolved native-reader difference. The Si V reader
comparison resolved 695 of its 801 points; the rest lie below its declared
subnormal comparison threshold and are not a physical rate floor. Both full
source-versus-fit and error PDFs were inspected. The private full Si V
source-node atom contains 169 levels, 3,778 positive-A rows, 14,196
electron pairs, and 212,940 source nodes. On 401 temperatures from 5,000
to 100,000 K and four electron densities (1e-6, 100, 1e6, and 1e10 cm^-3),
the compiled native evaluator differs from direct source collision strengths
by at most 3.34e-15 relative, from independent cooling profiles by at most
2.84e-14 relative, and from independent level populations by at most
1.23e-15 absolute. Its full cooling and difference PDFs in
`validation_2026-10-10/existing_remaining_cooling/SiV_source_node/products/`
were visually inspected. All positive-A source rows have a strictly
positive selected physical gap. These numerical checks do not establish
published accuracy below 5,000 K or authorize production use at 1,000 K.
Mg III's reversed A assignments still require physical source review.
A separate C IV source-valid 3,200–100,000 K
fit completed on 801 points, with 0.16338% maximum difference and
4.45e-16 maximum resolved native-reader difference; its full profile PDF
was inspected. This does not cover 1,000–3,200 K. A below-source
temperature treatment would need a separate physical decision and
implementation; the current Burgess-Tully reader has no constant-Upsilon
clamp for C IV.

S V is a separate missing production file with an existing staged fit at
`validation_2026-10-09/stationary_cooling_exponentials/products/SV/`;
this audit did not regenerate that product.
