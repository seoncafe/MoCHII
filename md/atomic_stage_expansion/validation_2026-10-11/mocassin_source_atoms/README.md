# Native validation of MOCASSIN source excitation atoms

54 candidate atoms passed the actual native reader, effective-strength
evaluator, stationary populations, line-emission consumer, and cooling
consumer. The candidate files preserve supplied levels and positive A and
Upsilon nodes. They are not installed in production. Mg I (one collision
temperature) and Ca VII (four nonpositive photon-energy pairs) require
separate physical review. The corrected Fe II model has 142 levels and
1,144 positive spontaneous probabilities; Fe VI retains all 19 levels.

The controller completed in tmux, with a finite 1,800-second limit and
OMP/BLAS restricted to one thread. Root verified all frozen source and
artifact hashes and recorded the numerical output hashes. See
[root_numerical_review.json](root_numerical_review.json).

| Quantity | Measured result |
| --- | --- |
| Candidate atoms exercised | 54 |
| Temperature/density cases | 5408 |
| Maximum Upsilon node relative difference | 3.10862446895e-15 |
| Maximum cooling relative difference | 1.36629596525e-11 |
| Maximum level-population absolute difference | 3.25083293617e-11 |

The reference is an independent dense stationary level matrix using the
source values and the adopted positive log-log PCHIP curve. It is not a
MOCASSIN executable run. All native population curves are saved in
`native/<ion>/profile.txt`. Complete independent curves are saved as NPZ
files with [54 vector comparison pages](figures/excitation_profiles.pdf).
Root inspected every page, including all Fe II population curves, and
verified the figure-controller inputs and all 55 output hashes; see
[root_profile_review.json](root_profile_review.json). The population display
uses a 1e-30 threshold only in the figure, with original values retained in
saved arrays. This review is not final production acceptance.

| Input condition | Native calculation | Independent level calculation |
| --- | --- | --- |
| Source levels and statistical weights | All supplied levels | Same supplied identities and energies |
| Spontaneous probabilities | All positive supplied A; corrected Fe II storage direction | Same emitting-to-lower transitions |
| Electron distribution | Maxwellian | Maxwellian |
| Temperature profile | 17 logarithmic nodes plus all source nodes from 1,000 to min(100,000, last source node) K | Identical |
| Electron density | 1e-6, 100, 1e6, 1e10 cm^-3 | Identical |
| Colliders | Electrons only | Electrons only |
| Cold continuation | Constant first Upsilon; actual temperature in detailed balance | Identical approximation |
| Radiation transfer and photoexcitation | Not included in this level calculation | Not included |
| Spatial domain, potential, spectrum, boundary conditions | Not applicable to this fixed-ion excitation calculation | Not applicable |

## Difference from the original v2 interpolation

The invoked v2 `emission_mod::equilibrium` uses a natural cubic spline of
Upsilon versus log10(T), then uses its value directly for both Maxwellian
rate directions. The candidate representation instead interpolates
log(Upsilon) versus log(T) with PCHIP, preserving positivity. This is an
intentional numerical correction, not an assertion that both codes give
identical results. In a separate 401-temperature sample, the original spline
expression produced 178 negative strength samples for N I and one for
Ne V. These are unphysical negative reaction rates at the expression level;
no full v2 plasma failure is inferred from this diagnostic.

Fe II has no sampled negative cubic values, but the maximum relative
difference of individual positive collision strengths reaches 1594.88179182
where the cubic value is small. This is not a total cooling difference.
The measured interpolation details are saved in
[interpolation_comparison.json](../../phase1_77_mocassin_v2/candidate_atoms/interpolation_comparison.json).
Absolute atomic accuracy, original publication correspondence, broad
physical temperature validity, and the remaining production integration
are unverified.

## Fe II source completeness limitation

The supplied v2 Fe II header identifies an allowed-transition model for
sextets and quartets. It explicitly states that 443 allowed quartet
transitions are set to zero because the listed A source omits their terms;
443 is a quoted header count, not a new measurement. The parsed raw A model
has 35 radiatively closed levels, while the current installed 13-level
MoCHII Fe II atom contains 40 positive forbidden probabilities. Its ground
fine-structure pair has A(2 -> 1)=0.00213 s^-1, absent from the supplied v2
allowed-only table.

Correcting the packed A direction does not restore these missing radiative
processes. The 142-level source candidate is a finite allowed-line model,
not a complete physical Fe II cooling model. Physical integration must
distinguish this supplied model from a model with identified forbidden
transitions and additional allowed transitions. It must not claim full
physical completeness from the level count or numerical matrix agreement.
