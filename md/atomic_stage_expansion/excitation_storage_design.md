# Excitation-atom storage implementation and validation

## Implementation status

The user authorized implementation of this storage change as a necessary detail
of the original ion-stage expansion, with autonomous continuation unless a
material change to the original plan is required. The dynamic atom arrays,
assumed-shape consumers, explicit version-2 source-node reader, and deterministic
source writer are implemented. The Intel debug production executable was rebuilt
successfully after the coordinated changes. Existing Chebyshev files remain
readable. No new dependency, user parameter, or cooling-grid dimension was added.

The five focused writer tests passed. S I's five-level source-node atom passed
compiled validation at 801 temperatures from 1,000 to 100,000 K and four electron
densities. Its maximum relative strength difference was 6.66e-16, all 150 sampled
source zeros remained exactly zero, and its maximum relative cooling difference
was 6.55e-14. These checks exercise the actual compiled populations, cooling
solver, and line-emissivity routine. Results and full profiles are preserved in
`validation_2026-10-10/source_node_native/SI/`. This is atom-level validation;
it does not establish a spatial photoionization benchmark or production
installation of the expanded roster.

The `atomic-source-native` controller completed all six declared source/native
steps, including Fe VII's 170-state atom. Ne V's separate corrected-source
controller also completed. Their original products and linked-object hashes
remain preserved. S I's five-state atom and S III's full 53-state atom passed
actual cooling registration and were installed. The physical density-ratio
correction and installation evidence are recorded in
[density-ratio cooling correction](density_ratio_cooling_correction.md).

The `atomic-registered-remaining` controller completed Ne V, Fe VI, Ar V,
Mg IV, and Fe VII. Each step built every one of the
64,561 table cells, then checks direct solves at table nodes and interior
locations. These checks establish numerical table behavior, not convergence
of a truncated physical atom outside 1,000--100,000 K.

Extended source diagnostics show that Ar V, Mg IV, Fe VII, and Ne V's earlier
truncated atoms do not meet the 0.3% target throughout 100--1,000,000 K. The
full-bound controllers checked atoms with 197, 157, 189, and 304 states,
respectively. Both Intel and GNU full profiles passed; all four GNU registered
tables passed. The last duplicate Intel Ne V table was deliberately stopped
after the equivalent GNU table passed. Historical objects and products are
preserved. See [extended validation](excitation_temperature_range_validation.md). Each profile or registration step has a
1,800-second limit and an 8-GiB address-space limit. Commands, source/build
hashes, process identifiers, and acceptance results are recorded under
`validation_2026-10-10/full_bound_native/`. Both numerical acceptance and
main-session physical review are required before installing these atoms.

The remainder below records the original design rationale and acceptance
requirements. Descriptions of fixed capacities refer to the implementation
before this change; they are retained as the measured reason for the redesign.

## Measured need and original linked path

The source counts and convergence profiles below were remeasured after the
nonnegative interpolation correction with the exact bound collision-connected
source atoms. The selected levels are the first passing tested candidates,
not a proof of minimality between sampled counts or native-runtime validation.
All comparisons use 801 temperatures from 1,000 to 100,000 K and four electron
densities, 1e-6, 1e2, 1e6, and 1e10 cm^-3.

| Ion | Selected levels | Radiative rows | Collision rows | Maximum sampled source-cutoff difference (%) |
|---|---:|---:|---:|---:|
| Fe IV | 37 | 391 | 645 | 0.000000 |
| Fe VI | 96 | 1,417 | 4,560 | 0.000000 |
| Ar V | 110 | 1,065 | 1,717 | 0.237096 |
| Mg IV | 110 | 1,656 | 5,995 | 0.270676 |
| Fe VII | 170 | 14,164 | 14,365 | 0.195270 |

Fe IV's full 37-state source has 645 collision pairs, exceeding the current
600-row limit even though its levels fit. The largest tested smaller atom
within that limit, 35 states, differs from the full source by 35.9707% at the
highest sampled density. Omitting two states is therefore not an accepted
replacement. Fe IV's missing production cooling was subsequently corrected after the
storage change and full native acceptance; its original failure products
remain preserved.

The source-cutoff target is 0.3%, reserving error for native numerical
representation. These are source-matrix diagnostics, not errors of a
serialized atom evaluated by compiled MoCHII. Complete profiles and arrays
are retained under `validation_2026-10-10/excitation_capacity_positive_source/`
and `validation_2026-10-10/existing_feiv_scope/`. Earlier
`excitation_capacity_study/` products are preserved as historical results.
The main session visually inspected the complete cooling profiles and the
cutoff-difference profiles for all four enlarged atoms at all four densities.

`src/nlevel_mod.f90:58-79` defines `MAXLEV=40`, `MAXTR=600`, `MAXHC=10`,
fixed member arrays in the public `nlevel_atom_type`, and 16 Chebyshev
coefficients for each collision row. `nlevel_load` rejects `NLEV`, `NRAD`, or
`NUPS` above those limits. `nlevel_cooling_prepT`,
`nlevel_cooling_solve`, `nlevel_populations`, `nlevel_emissivities`, and
`nlevel_hcoll_cooling` have 40-by-40 matrix and 40-element population
shapes; the emissivity output has length 600. `level_populations` currently
uses dense Gaussian elimination with partial pivoting.

The linked consumers are material to the change:

| Consumer | Current behavior requiring coordinated change |
|---|---|
| `nlevel_cooling_mod::nlevel_cooling_add` | Loads the same atom; builds 401 temperature by 161 density cooling cells using two fixed 40-by-40 matrices, one matrix assembly per temperature and one dense solve per density. `cooling_mod` and `species_mod` consume its registered table. |
| `lines_mod::lines_write` | Loads each atom and uses `Lline(600)`, `emis(600)`, and `wl(600)` before writing line rows for each transition and optional emissivity images. Its H I collision component also calls the same emitter. |
| `species_mod` | Loads atoms with H-impact rates for thermal corrections; form-9 charge-exchange product populations allocate 40-element populations and two 40-by-40 matrices. |
| Tests and fixed-state driver | `tests/ion_stages/check_excitation_atom.f90`, `check_bound_excitation_atom.f90`, and `check_niv_staging.f90` expose fixed matrix shapes; `tests/cloudy_c25/fixed_state_nlevel.f90` and `tests/hcoll/check_hcoll.f90` also allocate fixed emissivity vectors. |
| Build | `src/Makefile.dep` links `nlevel_mod` into `nlevel_cooling_mod`, `species_mod`, and `lines_mod`; `Makefile` has no LAPACK/BLAS link dependency. A type/routine-shape change requires rebuilding those module users and their test drivers. |

Changing only `MAXLEV`, `MAXTR`, or the reader would leave these other arrays
wrong. With 170 levels, each dense double-precision rate matrix occupies
about 226 KiB rather than 12.5 KiB at 40 levels. Dense elimination scales
approximately as the cube of level count: a 170-level solve has about 77
times the arithmetic work of a 40-level solve. Because the cooling table has
401 x 161 = 64,561 solves for each registered atom, source convergence alone
does not establish an acceptable production build time. Fe VII has all
170 x 169 / 2 = 14,365 pairwise collision rows, so sparse collision storage
does not by itself make its equilibrium matrix sparse. These are size and
complexity estimates, not measured runtime or memory peaks.

## Adopted representation

Allocate the atom's level, radiative, and collision arrays at their actual
`NLEV`, `NRAD`, and `NUPS` sizes after strict header validation. Represent
one collision row with its two level indices, Burgess-Tully type, `C`, `dE`,
source interval, interpolation choice, and source-node offset/count. Store
all scaled abscissae and ordinates in two contiguous atom-owned arrays with
row offsets; a varying number of nodes then costs only the nodes present.
The source values and their original logarithmic status are preserved: type
6 stores its logarithmic ordinates and descaling remains `10**y`, while
types 1-5 retain their stated transformations. Preserve original node order,
source identifier/hash, and the valid temperature interval as metadata in
the file header. The loader must reject duplicate or unsorted nodes, invalid
indices, nonfinite rates, and malformed counts before allocating large arrays.

Evaluate the scaled source curve at runtime with exactly the interpolation
decision made during preparation. For nonlogarithmic nodes that are all
nonnegative, use the original cubic only when all its interior extrema are
nonnegative; otherwise use PCHIP on those nodes. Type 6 retains cubic
interpolation of logarithmic ordinates. Negative source nodes are not
silently replaced: reject an atom whose requested temperature range reaches
an invalid negative-strength region. Clamp the transformed coordinate only
at the source's documented endpoints, then apply the Burgess-Tully descaling
and verify the resulting effective strength is finite and nonnegative.
This matches the present Python source rule in
`chianti_cooling.scaled_collision_strength_interpolator` and `upsilon`.
`nlevel_mod::ups_eval` currently clamps a possibly negative fitted result
with `max(ups,0)`; that behavior must not mask a source or interpolation
defect in the new representation.

Stout Cl VII and Fe VII require a separate source representation: their
physical input is a temperature/Upsilon table, evaluated by PCHIP in
`log(T)` and `log(Upsilon)`. Applying PCHIP directly to their type-2 storage
coordinate would change the adopted source curve. Store the original
13 temperature nodes and strengths with an explicit `LOG_PCHIP_T` model,
retain constant Upsilon below the 9,800 K first source node over the adopted
1,000--9,800 K approximation interval, and use actual T in Maxwellian rates.
Do not claim resonance accuracy below 9,800 K or extrapolate above the final
source node without a stated source decision. The numerical validation
interval and the source's physical validity interval are distinct metadata.

Introduce an explicit new format version in the `nlevel_*.txt` header with
node records and a source model (`BT_SOURCE` or `TEMPERATURE_TABLE`), and
an interpolation enum (`CUBIC`, `PCHIP`, or `LOG_PCHIP_T`). Keep the current
`NLEV/NRAD/NUPS` Chebyshev grammar as a legacy version, read by a distinct
validated branch when the new marker is absent. Both branches yield the same
in-memory physical transition
interface after loading. Do not infer a version from row length or silently
reinterpret an old coefficient row. The generator should reproduce the new
source-node format deterministically from the identified CHIANTI/Stout files
and record source hashes. Existing files remain readable; migrating an atom
requires explicit source and numerical gates, and must not regenerate
unrelated products.

Use allocatable rank-two matrices and rank-one populations sized by `nlev`
inside the n-level solver, cooling registration, form-9 charge-exchange
population path, and H-impact calculation. Make emissivity arrays size
`nrad`, including the line writer and fixed-state drivers. Keep the current
dense equilibrium solve initially, with explicit singular-pivot and
finite/nonnegative-population checks; benchmark its actual 170-level build
time, line time for each leaf, and peak memory before selecting a different
solver. A sparse representation is optional only if the assembled matrices
are measured sparse for the atoms actually loaded. No new LAPACK dependency
is justified by the current build evidence. If dense solve is too slow,
review a solver or cooling-grid redesign against physical accuracy and
runtime before adding a library or silently thinning the atom.

The cooling table storage (`Ssup` and `lamref`) remains allocated separately
for each ion. Its 401-by-161 numerical grid and `low_density`/
`local_ne` user choice do not change in this proposal. The line file and
optional emissivity image retain their present physical meaning and naming;
their transition dimension grows to the loaded atom's actual `nrad`.

## Why a single polynomial is insufficient

`fit_nlevel.py` samples the descaled source strength, fits a single
Chebyshev series, and stores at most 16 coefficients because the current
Fortran type has `coef(16,MAXTR)`. The new positivity-preserving PCHIP source
curve is piecewise polynomial, with potentially discontinuous second derivatives at nodes, that a single 16-term series
cannot uniformly reproduce to the 0.2% transition target. Direct execution
of the current `fit_transition` on Ne V 6->15 found 15 source nodes, type 2,
PCHIP interpolation, and a 4.3723% error with 16 coefficients. A separate
investigation reported about 1.06% as its best single-16 result, but that
search is reproduced in `validation_2026-10-10/nev_collision_fit_diagnostic/nev_pchip_profile.py` and its result files. It is the best tested coordinate, not a proof of a mathematical minimum. The older `NeV_native.log` reports 0.2835% for an
earlier source/fitting state and is stale for the current PCHIP comparison.
Source-node storage avoids claiming that any one
degree suffices. A piecewise Chebyshev representation is a viable alternative
only if segment endpoints, continuity, positivity, source-node agreement,
and worst error over the entire requested temperature interval are checked;
its extra fit and seam logic currently has no demonstrated advantage over
storing the existing nodes directly.

## Acceptance requirements

S I also requires direct source-node storage for a distinct numerical reason.
The unchanged CHIANTI source for transition 4--5 contains consecutive zero
scaled-strength nodes followed by a positive node. Its adopted PCHIP curve
therefore has a finite zero interval followed by positive values. A single
nonzero polynomial cannot reproduce an exact zero interval. The attempted
five-level native file was rejected before any population solve; it is
preserved in `validation_2026-10-10/si_native_correction/SI/`. Positive source
values should be checked relatively; source zeros must be checked explicitly
in absolute units without adding a rate floor. Tabulated zeros are source
data, not proof that the physical collision strength vanishes exactly.
That rejection described the initial polynomial candidate. S I's source-node
atom and cooling fit were subsequently installed after scoped native and
registered-table validation, as recorded in the
[current implementation status](implementation_status.md). A low-density
fit alone still does not establish acceptance of a finite-density atom.

The implementation should first pass focused parser and solver tests for both
file versions, malformed input, one physical transition of every Burgess-
Tully type, source interpolation choice, and H-impact coupling. A compiled
atom must reproduce the source effective collision strengths over the full
temperature interval, including node neighborhoods and interpolation
extrema, with finite nonnegative rates. Then compare compiled populations,
cooling, and line emissivities with independent source calculations over the
complete 801-temperature by four-density profiles for Fe VI, Ar V, Mg IV,
Fe VII, and Ne V. Inspect the full curves and relative-difference plots,
including strong lines and any low-density or high-density changes. Check
registered cooling versus direct compiled solves and the relevant fixed-state
driver. Finally measure table-build wall time and memory on the actual linked
binary. Source-convergence results alone, a successful build, or a golden
output do not pass these gates.

The authorized interface keeps user parameters and the meaning of existing
output quantities unchanged. New source-node files and larger emitted
transition lists are the visible data changes. Material additions to terminal
ionization or Auger physics remain distinct physical work and are not decided
by the storage implementation.
