# Explicit He I metastable population coupling

`par%hei_metastable_coupling` now defaults to `.false.` under the user's
instruction to compute additional physics only when explicitly requested.
This changes an existing input default; no input field or public API was
introduced. Explicit `.true.` retains the implemented population and energy
calculation. `metal_freefree` remains `.true.` by the user's exception.

## Source-path comparison

The inspected MOCASSIN v2 production Makefile links the source objects
containing `update_mod` and `emission_mod`; it does not link an external
helium-population physics library. `update_mod` forms ground ion ratios from
photoionization, recombination, and ground-state charge exchange.
`emission_mod.setDiffusePDF` includes triplet recombination emission using
`alpha2tS` and the He II ion density. These are present v2 processes and remain
available in MoCHII's ordinary cascade calculation.

The inspected v2 source has no invoked counterpart of MoCHII's local 2^3S
population feedback: excited-state photoionization, electron-impact ionization,
Penning ionization of hydrogen, triplet transfers in gas heating, and the
population-dependent modification of the diffuse branching fractions.
The distinction is population feedback, not the mere existence of triplet
emission. In MoCHII the existing flag controls the consumers in
`ion_balance_mod.solve_ion_cell`, `thermal_mod.net_rate`, `diffuse_mod`, and
`ion_score_mod`. The diagnostic-only `hei_metastable` flag is separate.

Ground-state He-H charge exchange is present in v2's `chex(2,1,:)` and
`chex(2,2,:)` coefficients and invoked ion-ratio denominators and reverse
terms. MoCHII retains `helium_charge_exchange=.true.`; its corrected physical
coefficients are not replaced with the reference's older fits. The separate
charge-transfer heating option remains disabled by default.

## Actual validation

The full debug production target rebuilt because the default belongs to the
shared parameter type. The actual implicit namelist reader passed seven
cases: omitted coupling, explicitly enabled, explicitly disabled,
diagnostic-only, explicit FUV, coupling with explicitly disabled FUV, and no
helium. For every disabled-coupling case, the invoked ion solver with a leaf
argument exactly matched the solver without that argument. This checks that
it does not accidentally use an excited-state estimator when disabled.

With both coupling and diagnostic omitted, the implicit FUV extension is now
off and the unused FUV floor remains 6 eV. An explicitly requested metastable
population still selects FUV and 4.78 eV if those input fields are omitted.
Explicit FUV selections remain honored. The explicit-no-FUV test verifies
input preservation only: it does not certify the completeness of an
excited-state photoionization integral with that restricted photon band.

The existing explicitly enabled coupled-physics check also passed: collision
nodes and detailed balance, recombination feeding, the trace-population
approximation, gas/radiation/ionization energy accounting, diffuse branch
weights, and solved H/He balances with and without ground-state charge
exchange. Measured worst energy closure was 3.68e-16; diffuse low-density
limit error was 2.49e-13; final H/He balance residual was 6.37e-11.
These are numerical checks, not a new published-data accuracy certification.

[Root review](root_review.json) verifies frozen sources and all controller
artifacts. [Coupled-physics output](native/coupled_physics.log) retains actual
measurements. The controller completed both steps and is no longer running.
No spatial model or full regression suite was run: the full build covers the
shared declaration, while the tests follow its input, population, energy,
and diffuse consumers. Earlier runs that omitted the now-disabled option
represent the older default and are not reference outputs for this choice.

Production stage expansion, shared-ion atom selection, and shell-rate choices
remain incomplete. This default correction does not complete the full goal.

The root README and physics manual now state the changed default. The manual
was compiled twice; PDF pages 26 and 27 were visually inspected for the edited
paragraph, equations, margins, and the following text. No overlap or clipping
was found. Existing bookmark and underfull-box warnings outside the edited
paragraph remain; a complete new manual content audit is not claimed.
