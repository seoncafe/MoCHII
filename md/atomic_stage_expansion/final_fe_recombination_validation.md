# Final-journal Fe recombination implementation

## Physical choice and implementation

The user supplied Zhang et al. (2026), ApJS 285:50, DOI
10.3847/1538-4365/ae80c0. Direct inspection of the final Tables 6/7 and
Eqs. (4)--(6) verifies the selected ground-level S1 DR and RR coefficients
for recombining Fe3+ through Fe6+. These replace inappropriate older
low-temperature DR fits; an absent Badnell record is not a zero rate.

`make_element_data.py` selects the published coefficients for transition
indices 3--6. Index q describes capture by Fe^(q+) into Fe^((q-1)+).
The source domain is 20*q^2 through 1e7*q^2 K. The fits assume an initial
ground level and zero density; they do not supply metastable-parent capture
or DR suppression. Evaluation outside the domain remains an uncertified
analytic extrapolation. The compact DR fit is accepted with its measured
source discrepancy; no new fit to the author data was silently substituted.

Fe V generation now succeeds with the existing reviewed KF96 Fe4 total
charge-transfer estimate and KF99 mean heat. Fe VI/VII still fail the
source preflight because their higher-charge transfer prescription remains
unselected. No expanded production roster is activated. (Corrected
2026-10-11, audit Section 5 item 7: stale. The Langevin rows for Fe5+ and
Fe6+ + H exist (`data/atomic/element_fe.txt`, used with
`par%langevin_charge_exchange`, default `.false.`), and stages I-VII of Fe
are production since 2026-10-11; see
[the restoration results](../ATOMIC_PHYSICS_RESTORATION_RESULTS_2026-10-11.md),
Phase 3.)

## Scoped measured checks

Seven focused Python tests passed, covering ADF09 ground-parent identities,
final PDF identity, positive rates at the source-domain endpoints, invalid
temperature rejection, visible compact-fit discrepancies, source generation,
and refusal to write unresolved added boundaries.

The persistent `mochii-final-fe` controller completed its generation and
native validation steps. The actual linked `species_fractions` calculation
read a seven-stage coefficient fixture. Known photoionization rates and
adjacent population ratios recover RR+DR through the private `alpha_rec`
call path. All 1,204 values (301 for each q=3--6, spanning each source
domain) match the published equations with maximum relative difference
1.0436096431476471e-14. Direct and cached electron-charge sums also agree.
CI and charge transfer are disabled only in this isolated coefficient test;
the fixture does not constitute a physical seven-stage Fe plasma.

| Comparison input | Compact fit and author ADF09 |
|---|---|
| Element and initial parent | Same Fe ion and parent index 1, checked ground level |
| DR calculation | S1 in both sources |
| Density treatment | Zero-density coefficient; no suppression |
| Temperature | Author nodes within 20*q^2--1e7*q^2 K |
| Unavoidable difference | Published compact fit versus author tabulated totals |
| Radiation, geometry, gravity, boundaries | Not applicable to a coefficient comparison; no nebular calculation |

Measured maximum absolute relative DR discrepancies on valid author nodes
are 4.7176%, 4.7126%, 5.1192%, and 3.0684% for q=3,4,5,6, respectively.
Their temperatures are 9,000, 320, 500, and 1,800 K. Consequently, the
published 2% accuracy statement is not certified by these measurements.
The full curves and all valid author nodes were plotted and visually
inspected in [the vector figure](validation_2026-10-10/final_fe_recombination/products/verified_fe_recombination.pdf).
No clipping hides the discrepancy or the temperature range.

## Existing Fe IV correction

The production Fe file's transition 3 previously used an older
high-temperature DR fit. At 10,000 K it gave DR=5.7883e-23 cm^3/s,
where the verified ground-level S1 fit gives 1.4762e-11 cm^3/s. The RR
component is now 2.6502e-12 cm^3/s. This is a physical correction to
missing low-temperature capture resonances, not a byte-identity choice.

Only the RR/DR records of this existing transition and their source comments
were replaced. All other records and the four-stage Fe layout were preserved.
The previous file is retained in the validation products. Total RR+DR changes
by factors 2.74389, 5.24451, and 0.357503 at 1,000, 10,000, and 100,000 K.
Older Fe IV recombination-dependent fractions, thermal results, and quoted
benchmarks are stale for reuse until recalculated. Old products are preserved.

The actual 43-stage production files were then read by the linked
three-temperature charge and mean-heat test at 5,000, 10,000, and 30,000 K.
Direct/cached electron-charge sums, reaction sums, and kinetic heating checks
passed. The largest absolute electron-charge difference was 2.7756e-17.
This verifies affected consumers without claiming a transported nebular result.
No full regression suite or full plasma profiles were calculated: stage
activation, missing cooling, terminal CI, and Auger physics remain downstream.

The captured-electron cooling consumer was checked separately against analytic
Maxwellian RR energy moments and resonance-energy-weighted DR moments. All
301 temperatures from 1,000 through 100,000 K passed, with maximum relative
difference 3.3417408840108465e-9. This directly exercises
`metal_ionization_cooling`, including its derivative of the new total rate.

A 54-stage, 43-transition controlled roster including Fe V also passed the
linked direct/cached charge and mean-heat checks at three temperatures.
Fe4's adopted rate differed from the KF96 expression by at most 3.41477e-6
relative; its adopted mean kinetic heat is 10 eV. This roster contains only
the reactions with selected sources. It is a rate fixture without complete
line cooling and is not activated as production data.

Evidence is in
`validation_2026-10-10/final_fe_recombination/controller/state.json`,
`products/fortran_result.json`, `products/adf09_audit.json`,
`products/production_correction.json`, and `products/production_thermal_closure.txt`.
The controller has completed and is no longer running.
