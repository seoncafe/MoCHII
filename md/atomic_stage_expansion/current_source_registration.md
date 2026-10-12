# Current source-node cooling registration

The corrected absolute-cooling module passed actual registration and public
application checks for staged C V, N VI, and Ar VI atoms. Each calculation
built 64,561 cells; this does not activate these ions in production.

| Ion | Levels | Maximum interior table difference (%) | Maximum application difference (%) | Table build wall time (s) |
| --- | ---: | ---: | ---: | ---: |
| CV | 10 | 0.115760 | 0.113258 | 0.101282 |
| NVI | 10 | 0.002581 | 0.001781 | 0.096269 |
| ArVI | 45 | 0.307473 | 0.310620 | 3.032099 |

The finite source-convergence acceptance covers 1,000--100,000 K and electron
densities 1e-6, 100, 1e6, and 1e10 cm^-3. The numerical registry extends from
100 to 1,000,000 K; successful interpolation throughout this grid does not
establish physical source accuracy or truncated-atom convergence outside the
accepted interval. In particular, 40 Ar VI levels miss the original 0.3%
source-convergence target; its 45-level replacement has a maximum measured
source cooling difference of 0.225088% within the accepted interval.

The main session inspected complete source and native Ar VI profiles and the
failed 40-level difference curves. The Ar VI coronal fit was preserved because
the collision-interpolation correction did not justify regenerating it.

The controller `atomic-added-current` completed its build and all three ion
checks. It uses an isolated GNU O2 build with bounds checks and invalid,
zero-division, and overflow traps. Commands, source/input/executable hashes,
resources, full native check outputs, and numerical acceptance are preserved
in [controller state](validation_2026-10-10/current_added_registration/controller/state.json).
The public check's legacy relative-error metric excludes deeply underflowing
samples; dedicated positive-subnormal tests provide the separate verification
of that numerical behavior. No spatial temperature, density, velocity, or
composition profiles were calculated for this registration check.
