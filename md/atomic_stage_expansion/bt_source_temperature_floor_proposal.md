# Physical temperature limits for Burgess–Tully source atoms

Date: October 10, 2026. This is an implementation proposal, not a change to
the source evaluator or production data. The adopted low-temperature rule is
to hold the **effective collision strength** constant below the lowest
published source temperature. The Maxwellian rate coefficient, Boltzmann
factor, detailed balance, and energy loss must still use the actual electron
temperature and selected physical level energies.

| Ion | Collision calculation | Published source interval | Proposed lower bound |
| --- | --- | ---: | ---: |
| C IV | Liang and Badnell (2011), Li-like R-matrix with Auger and radiation damping, *A&A* **528**, A69, DOI [10.1051/0004-6361/201016417](https://doi.org/10.1051/0004-6361/201016417) | `200(z+1)^2` to `2e6(z+1)^2` K, residual charge `z=3` | 3,200 K |
| Mg III | Liang and Badnell (2010), Ne-like intermediate-coupling frame-transformation R-matrix, *A&A* **518**, A64, DOI [10.1051/0004-6361/201014170](https://doi.org/10.1051/0004-6361/201014170) | `200(q+1)^2` to `2e6(q+1)^2` K, residual charge `q=2` | 1,800 K |
| Si V | Same Ne-like paper, `q=4` | Same formula | 5,000 K |

These intervals are stated in the published-paper abstracts and, for C IV,
the local CHIANTI `.scups` header explicitly prints 3,200–32,000,000 K.
The local Mg III and Si V `.scups` headers identify the Ne-like paper but do
not print temperature limits. An abscissa of zero among the scaled
Burgess–Tully nodes is a mathematical interpolation endpoint, not evidence
that the scattering calculation establishes Upsilon at zero kelvin. The
chosen constant continuation below the paper's minimum is an approximation,
not an accuracy claim from either paper. The use of each source above its
published maximum would remain a separate issue.

## Minimal evaluator and file representation

Keep every original `.scups` `xs`, `ys`, Burgess–Tully type, scaling constant,
and positive scaling energy unchanged. Associate the physical lower bound
with each identified ion in the Python reader, e.g.
`BT_SOURCE_TMIN_K = {('c', 4): 3200., ('mg', 3): 1800.,
('si', 5): 5000.}`. Put the selected value on each returned collision record
as `bt_tmin_K`, with absent metadata meaning the current evaluator behavior.
Within `chianti_cooling.upsilon`, for Burgess–Tully records only, form
`T_source = maximum(T_actual, bt_tmin_K)` and use `T_source` in **both** the
scaled-coordinate transformation and type-dependent descaling. In particular,
type 5 divides by the `et` computed from `T_source`; clamping only its scaled
coordinate would not produce constant Upsilon. The temperature-table source
path retains its separate first-node behavior. `full_excitation_balance` and
other rate builders continue to pass `T_actual` to the Maxwellian factors.

For version-2 native files, add one optional atom-wide data record
`BT_TMIN_K value` after all `NUPS` rows, before `DENSITY_RATIO`. A version-2
file without that record retains its existing meaning; legacy files never
accept it. The writer emits this record only for the three identified atoms
and states the paper interval and adopted cold continuation in adjacent
comments. `nlevel_load` stores one finite, nonnegative atom value and rejects
duplicates or a `BT_TMIN_K` record in a source without BT rows. In
`nlevel_upsilon`, use `max(T_actual, atom%bt_tmin_K)` only when
`source_model(k)==BT_SOURCE`, including both the scaled-coordinate and
descaling calculations. `nlevel_cooling_prepT` retains `T_actual` for
`sqrt(T)`, the Boltzmann factor, and detailed balance. No public plasma
parameter changes. An explicit file record is preferable to inferring an
ion from a filename or silently changing the meaning of all version-2 atoms.

At the source minimum, the old and new Upsilon curves must coincide exactly.
At lower temperatures, compare Python and native Upsilon at multiple points
for each Burgess–Tully type present, especially type 5 and type 6. Check
Maxwellian upward/downward rates independently at actual temperatures and
inspect full cooling and population profiles over 1,000–100,000 K at the
four established electron densities. Existing version-2 and legacy fixtures
without the new record should retain their current rates. These numerical
checks establish implementation agreement, not the cold continuation's
physical accuracy.

## Source access and radiative-data findings

An ADS API search verified the bibcodes, titles, DOIs, and explicit
open-access publisher PDF links for both Liang–Badnell articles. The sandbox
could not resolve `api.adsabs.harvard.edu`, so the API request was repeated
successfully with authorized network access. Both A&A final publisher PDF
links returned HTTP 403 even with that access, as did the OUP final PDF for
Dufresne et al.; an EDP Sciences host did not resolve. No final journal PDF
for either Liang–Badnell article was in the local references directory.
The University of Strathclyde repository records for the
[Ne-like article](https://strathprints.strath.ac.uk/28926/) and
[Li-like article](https://strathprints.strath.ac.uk/31247/) explicitly say
that their full text is unavailable there. The author's
[APAP publication list](https://amdpp.phys.strath.ac.uk/UK_RmaX/APAP_pub.html)
links these entries to the DOI/publisher rather than to a deposited final
PDF. A similarly named Strathprints PDF under record 19077 is the different
2009 **Na-like** paper and was not substituted.
Thus the detailed published methods and any
accuracy claims beyond the abstracts remain **unverified**; the final
articles are Liang and Badnell (2010), *A&A* **518**, A64,
DOI 10.1051/0004-6361/201014170, and Liang and Badnell (2011), *A&A*
**528**, A69, DOI 10.1051/0004-6361/201016417. Their methods sections and
any source-data exceptions must be checked in the final publisher articles
before claims about precision are made.

The Mg III `.wgfa` header attributes radiative rows within levels 1–27 to
Jönsson et al. (2014), *Atomic Data and Nuclear Data Tables* **100**, 1,
DOI [10.1016/j.adt.2013.06.001](https://doi.org/10.1016/j.adt.2013.06.001),
and the remaining radiative data to Liang and Badnell (2010). All 16
positive-A rows reversed by the current mixed selected energies have indices
above 27 and negative printed wavelengths, which the CHIANTI format marks
as theoretical. For example, row 128→131 has theoretical energies
620,751→623,762 cm^-1 and wavelength −33,211.558 Å, while the selected
observed energy of 131 is 618,601.312 cm^-1 and the selected lower energy is
theoretical 620,751 cm^-1. Across all 16 inverted rows, the printed negative
wavelengths reproduce theoretical energy gaps to within 3.2e-8 relative,
and the standard E1 gf-to-A relation using the listed upper statistical
weight agrees within 0.36%. The A, gf, and wavelength tuple is internally
theoretical, even though mixed selected energies invert those rows.

Using observed energies for levels 1–27 and theoretical energies above 27
would remove every reversed radiative and collision row while retaining two
exactly degenerate collision pairs. It is **not** a coherent physical photon
energy assignment for the full atom: the present selected energies reproduce
all 1,081 positive-wavelength radiative rows within 2.45e-6 relative, whereas
the proposed hybrid would disagree by more than 1% for 580 of those rows
(median 1.10%; maximum 96.9%). Consequently, source A orientation and the
choice of physical photon energies need row-level assessment before a native
Mg III atom is adopted. Neither the A values nor source level identities
were changed.

The Si I `.wgfa` header cites Tachiev and Froese Fischer (2002), *A&A*
**385**, 716, DOI
[10.1051/0004-6361:20011816](https://doi.org/10.1051/0004-6361:20011816),
but that article explicitly covers **nitrogen-like and oxygen-like**
sequences; it cannot identify a neutral Si (Si-like) radiative calculation.
The [final published methods](https://academic.oup.com/mnras/article/505/3/3968/6283727)
of Dufresne, Del Zanna, and Storey (2021), *MNRAS* **505**, 3968,
DOI [10.1093/mnras/stab1498](https://doi.org/10.1093/mnras/stab1498),
attribute the Tachiev radiative addition to **N I** and describe the **Si I**
decay model as MCHF rates from Froese Fischer (2005), with additional ATSP2K
calculations. The Si I electron collision model uses Autostructure distorted
waves; the authors say that neutral distorted-wave strengths are not usually
reliable for line intensities and used them to estimate metastable populations
for ionization. The CHIANTI header does not identify the origin of anomalous
Si I row 15→14. Its observed energy order, wavelength magnitude, gf, and
A=132.8 s^-1 remain inconsistent under a simple index reversal. No Si I A
was reassigned. The final publisher HTML was read, but its PDF was not
available for a `pdftotext` check in this environment.

## Production build and focused native validation

The implementation was rebuilt through the actual production target with
`make DEBUG=1 -j4`. The completion controller validated the executable and
linked objects, then tested all six Burgess--Tully types against independent
source strengths, level populations, detailed balance, and cooling. The
maximum relative strength and cooling differences were 4.44e-16 and 6.66e-16,
respectively. Eleven malformed floor records were rejected explicitly; old
version-2 files without the record and legacy files remained readable. These
focused native tests used `-no-ftz`; production compiler defaults were not
changed. No spatial plasma calculation was run for this validation.

Measured results and frozen source/object hashes are in
[the production validation](validation_2026-10-10/bt_floor_production_validation/result.json)
and [build record](validation_2026-10-10/bt_floor_production_validation/build_result.json).
This establishes the implemented cold-strength rule, not the accuracy of
source atomic calculations or completion of the 64-stage production expansion.
