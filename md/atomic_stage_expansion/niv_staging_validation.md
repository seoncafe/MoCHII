# N IV energy selection and staging validation

## Physical model

The final publisher article supplied by the user was verified in
[the source report](niv_final_paper_verification.md). The implemented shared
selector uses Cloudy c25.00 Stout `n_4.nrg` energies while retaining the
original CHIANTI Einstein A values and electron collision identities. The
same energy gaps enter photon energies, collision excitation Boltzmann
factors, and detailed balance. The selector verifies the Stout hash, 238
matching statistical weights, 3,272 radiative pairs, 9,180 collision pairs,
and their positive energy gaps before selecting a smaller model.

Observed energies describe levels 1--29, with BSR energies above. The
external Stout header's claim of 20 updated levels is inaccurate; the
original external file remains preserved. The observed energies formerly
mixed into high CHIANTI levels inverted some transition directions. The
correction keeps each Einstein A attached to its original emitting level.

This remains an intentional approximation to the published scattering
calculation. Near-threshold resonance positions are uncertain, especially
at low temperature. The paper does not establish target-state convergence
above level 87, and original printed BSR indices cannot be directly equated
to the selected 238-state indices above that level. Neither a small numerical
fit error nor a model-cutoff comparison certifies physical accuracy.

## Model comparison and selected size

| Input | 29, 87, and 136-state comparisons |
|---|---|
| Ion | N IV |
| Energies and statistical weights | Same checked source and indices |
| Radiation | Same spontaneous A values within each cutoff |
| Collisions | Same CHIANTI effective strengths within each cutoff |
| Population | Stationary electron collisions and spontaneous radiation |
| Density | Analytic limit n_e -> 0; independent matrix at 1e-12 cm^-3 |
| Temperature | 1,000--100,000 K, 801 points for full curves |
| Changed input | State cutoff and transitions that connect retained states |
| Geometry, gravity, boundary spectrum | Not applicable to isolated atomic coefficients |
| Excluded processes | Photoexcitation, excited-state ionization, electron-impact ionization, transported radiation |

The maximum relative cooling deficits against the supplied 136-state model
are 0.1375974% for 29 states and 0.0102918% for 87 states. Both maxima occur
at 100,000 K. The complete curves were plotted and visually inspected in
[the comparison figure](validation_2026-10-10/niv_staging/products/niv_level_cooling.pdf).
The independent 136-state finite-density matrix agrees with the analytic
limit at 1,000, 10,000, and 100,000 K to floating-point precision, with
finite, nonnegative, normalized populations.

The selected staged atom contains 29 levels, 170 radiative transitions,
and 406 electron-collision transitions. It fits the existing native
`MAXLEV=40`, `MAXTR=600`, and 16-coefficient capacities. No native storage
or file format was changed. The 29-state choice is supported for this
cooling comparison over the numerical interval, not certified for every
individual high-level diagnostic line or every density.

## Fit correction and native checks

The rounded low-density cooling fit has maximum relative error 0.4191986%
against the stationary 136-state source over 801 temperatures. Its numerical
target is 1%. The original collision-strength fitting method failed its
0.2% target: a low-level transition reached 0.76867%. The first controller
stopped before native execution, and its failed products remain preserved.

The fitter now attempts a Chebyshev minimax approximation when its original
least-squares representation misses the target, using at most the existing
16 coefficients. This changes numerical approximation rather than physical
collision strengths. No tolerance was relaxed. The original accepted fits
retain their former computation path. A focused test verifies the three
low-level strength peaks on 1,601 temperatures.

All 406 serialized collision fits were independently checked on 1,601
temperatures. The maximum relative error is 0.1999802%, with no transition
above the 0.2% target. The recovery controller reused the already checked
cooling fit and curves, preserving their hashes, and regenerated only the
affected n-level product.

The actual Fortran reader and n-level solver passed 2,403 cases: 801
temperatures at electron densities 1e-6, 100, and 1e6 cm^-3. All 29
populations were finite, nonnegative, and normalized. The cooling-fit
reader/evaluator differs from Python by at most 1.00295e-11 relative.
The native 29-state low-density cooling differs from the original spline
model by at most 0.1358368%; the independent finite-density matrix differs
by at most 0.1346944% at the three temperatures and three densities tested.
These differences include the collision-strength approximation.

The products remain in staging. No expanded production roster, density
interpolation table, transported nebular model, or matched code comparison
has been activated or calculated. Broader tests would not resolve the
unselected higher-charge reaction probabilities or missing Auger physics.

## Evidence and execution state

`validation_2026-10-10/niv_staging/controller/state.json` preserves the
initial collision-fit failure. `niv_staging_minimax/recovery.json` records
the cause, correction, and copied hashes. Its controller has completed;
`products/native_result.json` and `products/dense_collision_audit.json`
record the native and serialized-fit checks. The separate
`niv_energy_checks/finite_density_level_comparison.json` records full
populations, formerly tied gaps, and selected high-level line contributions.
No calculation controller from this sequence remains running.
