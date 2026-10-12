# Additional high-excitation cooling validation

Date: October 10, 2026. These are staged collisional line-cooling fits;
production files and the active ion roster were not replaced.

## Corrections confirmed against the implementation

CHIANTI Ne IV has downward radiative transitions whose upper-level index is
smaller than the lower-level index. Indices identify levels; they do not
establish energy ordering. The stationary cascade calculation now preserves
the supplied upper/lower identities and checks the actual energy direction
where the energy-resolved closed-state balance requires it. It does not flip
radiative rows. A synthetic two-end-state atom with reversed index ordering
checks cascade probabilities independently.

C V, N VI and S VII have ground excitation thresholds high enough that cold
Boltzmann factors underflow double precision. The source is now evaluated
in log space when checking these one-end-state excitation sums. No physical
rate floor is introduced. Signed exponential fits remove the physical threshold
and resolve the remaining temperature dependence using excitation-temperature
offsets. Signed terms represent a positive total mathematical fit; they are
not separate negative emission rates. Full-range fit tolerance remains 1%.
The h c/k_B constant is derived from the same h c and k_B used for excitation
energies, eliminating a preexisting rounding inconsistency.

## Measured acceptance

Every fit was compared to its own stationary low-density source at 801
logarithmically spaced temperatures from 1,000 to 100,000 K. The logarithmic
error remains defined where an ordinary floating-point cooling value is zero.
The same ions, energy vector, collision data, radiative transitions, optically
thin escape, Maxwellian electrons and n_e -> 0 population assumption are used
on both sides. No finite-density, photoexcitation, or nebular comparison is
claimed. The fit domain does not certify the underlying atomic accuracy.

| Ion | Maximum relative fit error (%) | Ordinary source zero points |
|---|---:|---:|
| CV | 0.8690102 | 277 |
| NVI | 0.7239260 | 337 |
| SVII | 0.4071741 | 179 |
| NeIV | 0.0585380 | 0 |

Native linked debug `cooling_load` / `cooling_eval` checks passed all 3,204
samples. The largest relative difference in the well-resolved normal range
is 2.6645352591003757e-15. Near the normal floating-point boundary, each
coefficient multiplication can underflow in the Intel calculation. Acceptance
there uses an absolute roundoff bound of n_terms times the smallest normal
floating-point value, plus the relative tolerance; this does not alter a
coefficient or impose a physical cutoff. Native zero counts differ from NumPy
because their subnormal-number handling differs. The logarithmic source check
provides the physical fit assessment for this tail.

Evidence: [native result](validation_2026-10-10/high_excitation_cooling_offsets/native_underflow_bound/native_result.json)
and [completed fitting controller](validation_2026-10-10/high_excitation_cooling_offsets/controller/state.json).
The initial failed C V fit (2.6154%, above the unchanged 1% requirement) and
native formatting/underflow diagnostic products are preserved. Three-digit
exponents are now printed explicitly by the native validation driver.
Eight focused stationary/underflow tests pass, in addition to the 13 capture
and generation tests. All four complete temperature curves and residual plots
were visually inspected; no clipping of labels or overlapping plot elements
was found. Very cold source values are mathematically positive but not visible
in an ordinary linear floating-point logarithmic cooling plot; residuals still
cover the entire checked interval.

Vector figures:
[ C V](validation_2026-10-10/high_excitation_cooling_offsets/products/CV/CV.pdf),
[N VI](validation_2026-10-10/high_excitation_cooling_offsets/products/NVI/NVI.pdf),
[S VII](validation_2026-10-10/high_excitation_cooling_offsets/products/SVII/SVII.pdf),
[Ne IV](validation_2026-10-10/high_excitation_cooling_offsets/products/NeIV/NeIV.pdf).

## Original energy-source issue and subsequent resolution

The original generation required coherent O V and Cl VI level energies.
Subsequent verified energy selection and native tests resolved that issue;
see [ICFT energy and atom validation](icft_energy_atom_validation.md).
The following counts describe the originally inspected source differences. The CHIANTI O V atom
contains 166 retained bound levels, while the corresponding Stout atom has
238; direct index copying already disagrees in statistical weights starting
at index 87. Cl VI has 171 retained CHIANTI levels versus 283 Stout levels;
matching statistical weights alone does not establish all level identities.
O V theoretical energies contain two rounded zero collision gaps. Mixed
observed/theoretical energies reverse five O V and one Cl VI radiative rows.
Neither source rows nor external energies are silently reordered or repaired.
Final 2014 Be-like and Mg-like journal articles were obtained centrally and
their structure/scattering methods were read; see [source provenance](reference_sources.md).


## Current interpolation correction

The C V and N VI fits above predate the selective PCHIP correction and are
historical products. They exceeded the 1% source-fit target when evaluated
against the corrected source. Only these affected fits were recomputed.
The current 801-point errors are 0.855934% for C V and 0.796310% for N VI;
actual serialized native-reader checks passed. The main session inspected
the complete curves and residuals in
[updated C V](validation_2026-10-10/positivity_correction_cooling_refit/CV/CV.pdf)
and [updated N VI](validation_2026-10-10/positivity_correction_cooling_refit/NVI/NVI.pdf).
Their original ten-state atoms contain affected collision pairs and are
being revalidated with the source-node representation. Updated cooling fits
alone do not validate those original finite-density atoms or activate the
expanded production roster.
