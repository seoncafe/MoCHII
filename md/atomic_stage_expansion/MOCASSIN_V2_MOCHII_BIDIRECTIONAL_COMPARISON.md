# MOCASSIN v2 and MoCHII: differences and implementation status

> **Dated record of the October 9-11 tree, corrected 2026-10-11.** The current
> state is [ATOMIC_PHYSICS_RESTORATION_RESULTS_2026-10-11.md](../ATOMIC_PHYSICS_RESTORATION_RESULTS_2026-10-11.md):
> 77 stages / 66 transitions in production, `phfit2` subshell absorption
> restored, default `metal_excitation_model = 'chianti'`, HEAD's physics on
> by default, and the removed controls (`metal_inner_shell_absorption`,
> `additional_line_cooling`, `high_charge_exchange`, `fe_levels_full`)
> refused with their replacement named. Wrong statements below are corrected
> in place as [the audit](../ATOMIC_STAGE_EXPANSION_AUDIT_2026-10-11.md),
> Section 5, lists; "State note" marks a statement that was right then and is
> no longer. In this document "production" means the default code path and
> data of that tree, not a run-capable binary (audit Section 5 item 4).

Updated: October 11, 2026 (KST).

## Approved changes and current production status

The user approved the paired atomic-data replacements, coherent excitation
model selector, explicit inner-shell option, qualified constant upper
collision strengths, preserved Fe II forbidden radiation, and eventual
77-stage activation on October 11, 2026. Earlier approval-dependent statements
below describe historical checkpoints and no longer represent pending approval.

Six atom/fit pairs (C III, N III, O I, Ne III, S IV, and Si II) have now been
installed together. All twelve originals and hashes are preserved in
[the installation record](validation_2026-10-11/installed_shared_atom_pairs/result.json).
The selected reference cooling listing contains 49 ions and 763 terms;
the extended listing contains 38 ions and 377 terms. Their coefficient
values match the primary files. The historical root listing is not the
production metal selection. Primary cooling products were not regenerated.

Production `metal_inner_shell_absorption` now defaults to `.false.`. The
invoked metal photoionization API selects the physical outer shell; requesting
an active inner shell stops with an explicit diagnostic until vacancy charge
and electron-energy treatment is qualified. This corrects the former implicit
adjacent-stage treatment of every shell. The full production build passed.
(Corrected 2026-10-11, audit Section 5 item 2: it was not only a correction
of the product assignment. Above `einn` the outer-shell default kept a single
PH1 row and deleted valid absorption, e.g. Fe I 70 eV 8.12 -> 0.082 Mb,
Mg II 2p 70 eV 5.96 -> 0.067 Mb, Ca II 3p 45 eV 1.31 -> 0.067 Mb, Cl I-IV
-38 to -45% at all energies, with a discontinuity at `einn`; v2 keeps every
shell in its opacity (`ionization_mod.f90:409-438`), so this default was
below v2 as well as HEAD (audit E2). The explicit `.true.` request could not
recover it: it stopped at the first absorbing packet. Audit Section 5 item 4:
the "full production build" was a `DEBUG=1` build (`-O0 -check all -fpe0`)
that trapped on the NaN default `density_powerlaw_index` before any physics
and could not run the benchmark inputs. State note: the `phfit2` sum is
restored and the option and the routine `outer_shell_photoionization` are
removed.)
Actual packet, opacity, rate-score, heat-score, input-reader, and paired-data
consumer checks completed under the corrected persistent controller. The
initial fixture failure and completed continuation are recorded below.

The coherent source selector and approved hot endpoint models now pass the
full build and native validation: 55 reference atoms (571 levels), initially 61 extended
atoms (1,946 levels), 3,216 direct profile samples, 80 registered cooling
samples, and six rejection cases. See
[the production selection report](validation_2026-10-11/approved_excitation_models/README.md).
The 142-level Fe II atom with its 40 retained forbidden probabilities and
matching 29-term cooling fit has been installed. Its three actual installed
reader samples have maximum fit/emission difference 0.0267510%.

The earlier shell controller failed because a diagnostic assumed that
20 eV was below every inner edge. The production guard correctly rejected
Fe I's active 14.7 eV inner edge. The corrected controller completed 3,822
shell samples and 18 installed-pair reader samples; all source hashes were
verified. See its [root review](validation_2026-10-11/approved_shell_and_pairs_checked/root_review.json).

The prepared 77-stage rate chains now also pass 396 actual charge-transfer,
electron-accounting, population, and mean-heat samples with the selected
reference excitation model. The common H I reader path was preserved after
checking the shared cooling registration caller; a final full build and
reader-only diagnostics passed without rebuilding its numerical table.
See [the 77-stage connection](validation_2026-10-11/staged77_selected_model/README.md)
and [the hydrogen reader check](validation_2026-10-11/approved_excitation_reader_guards_checked/README.md).

The resumed selected-model inventory found that `extended` omitted the
approved Mg I atom. The identical file has now been included there, giving
62 extended atoms and 1,951 levels. Its initial 61-atom validation is
historical; the scoped correction uses the actual input and production
reader without regenerating numerical products. See
[the Mg I correction](validation_2026-10-11/extended_mgi_coverage/README.md).

Ca VII correction and 77-stage activation remain unfinished. Ca VII requires
correction of the confirmed energy and radiation inconsistencies; approval
does not replace physical validation. A newer conditional audit identified
111 pairs allowing one component and 199 pairs allowing two components
under its stated E1/M1/E2/M2 assignment. It cannot establish a complete,
physically accurate atom. The implementation goal is recorded as blocked after
three consecutive audits of this physical source condition, with the earlier
approval resolved. `metal_freefree=.true.` remains.
The old Cloudy c13.05/c17.03 CHIANTI files contain all 310 Ca VII radiation
rows, exactly matching the supplied source. Modern CHIANTI retains 86 after
a recorded branching filter. Recovering the rows does not resolve their
mixed configuration provenance or separate multipole energy corrections.
The checked official source routes did not provide the required original
25-configuration components. See [the source assessment](validation_2026-10-11/cavii_completion_research/README.md).
No matched spatial comparison with MOCASSIN v2 or complete thermal-equilibrium
calculation has been performed. See the [current result](../ATOMIC_IONIZATION_STAGE_EXPANSION_RESULTS_2026-10-11.md).

## Historical checkpoints preceding approval

The 401-node production PCHIP correction passed the full build, 797 compiled
interpolation samples, 1584 thermal samples, 1584 emitted-energy samples, and
24 zero-density checks. Root verified hashes and inspected all six figure
pages. No shared production atom was replaced. The proposed change now pairs
six atoms with six matching cooling fits, because replacing atoms alone would
mix different collision data in thermal and emitted-energy calculations.
The 792 derived matching-source comparisons have maximum discrepancy
0.457773%; the actual combined invocation now also passed 1584 thermal and 1584
emitted-energy samples. The affected O I option passed 18 checks.
Production paired replacement still requires the recorded approval.
See [the twelve-file review](SHARED_ATOM_EXPANSION_REVIEW_2026-10-11.md).

The Mg I constant-strength candidate completed 404 new hot checks and reused
404 verified cold samples. Its 1000-100000 K approximation remains unaccepted
and uninstalled; numerical agreement does not establish atomic accuracy.
See [the physical choice](validation_2026-10-11/mgi_constant_source/README.md).
`metal_freefree=.true.` remains the explicit user exception. Production rate
chains still contain 43 metal stages and 32 transitions; the target is 77 and
66. The global objective is incomplete.

The latest [Auger mean-energy source review](validation_2026-10-11/auger_mean_energy/README.md)
corrects an overly restrictive implementation condition: the current linear
full thermalization of metal electrons needs a consistent unconditional mean,
not a joint final-charge/energy distribution. The published KM mean correction
is negative for Ca II M1 (-5.022 eV); charge thresholds and energy conventions
remain unresolved. This source review does not activate inner-shell rates or
heating and does not change the requested `metal_freefree=.true.` default.

# Qualified Fe II matching cooling representation

The 142-level Fe II candidate with 40 preserved forbidden probabilities now
has a matching 29-term low-density fit. Actual Fortran evaluation passed 801
samples; its maximum source difference is 0.0299785841%. The full coefficient
and error profiles were visually inspected. Derived finite-density predictions
are explicitly distinguished from native thermal invocations in
[the representation report](validation_2026-10-11/feii_radiative_cooling/README.md).

The first finite-density attempt failed before table construction because of
a diagnostic argument-count error. Its corrected continuation completed 40 actual matching-fit application checks,
with a measured 48461-cell table and maximum direct-cooling difference
0.1173885499%. Twenty direct values agreed with the previous checked
compilation to 2.88658e-15; all four comparison pages were visually inspected.
See
[the finite-density report](validation_2026-10-11/feii_density_cooling_checked/README.md).
No production replacement or physical-default choice was made.
`metal_freefree=.true.` remains the requested exception.

# Qualified sulfur charge-transfer dependence on the excitation model

Actual native charge-transfer consumers passed 808 samples with the supplied
5-level and installed 70-level S II product atoms. Existing 4S, 2D, and 2P
ranges match both atoms, so these particular files need no term remapping.
With identical reaction tables, the maximum reverse-rate difference is
44.64932194%, and the maximum mean reverse-energy difference is 0.2184102926 eV.
These source-model differences are not a measured total cooling change.
The selected atom must therefore also reach the charge-transfer population
calculation. All four full-profile pages were visually inspected. See
[the completed comparison](validation_2026-10-11/sulfur_charge_exchange_models/README.md).

The existing rate equations were already correct. A stale comment claiming
that the represented terms sum to one was corrected; their measured minimum
sum is 0.9904478973 in the 70-level model. The subset must not be renormalized.
No public selector or atomic-data replacement was made. The global objective
remains incomplete, and `metal_freefree=.true.` remains enabled.

# Ca VII final publication received; radiation-model inconsistency confirmed

The supplied Landi and Bhatia (2003) final article has been read directly.
The previous article-access dependency is resolved. A 55-pair comparison
shows that Table 3 source A values match the five-configuration column,
whereas Tables 4 and 5 match the twenty-five-configuration column. The paper
adopts the larger radiative calculation and describes observed-energy
corrections. The four nonpositive source photon gaps and mixed radiation
model must be corrected together; an all-theoretical energy substitution
alone has not been accepted. See [the verified source review](validation_2026-10-11/cavii_published_review_checked/README.md).
No production data or public input changed. `metal_freefree=.true.` remains.

## Objective and evidence

Complete the revised first-phase 77-metal-stage roster (66 adjacent
transitions) and the explicit levels, electron collision strengths, and line
cooling represented by the supplied MOCASSIN v2 atom models. Correct physical
and numerical errors instead of reproducing them. Calculate additional MoCHII
physics only when the corresponding input option is explicitly enabled,
except `metal_freefree`, which remains enabled by the latest explicit user
instruction.
The comparison refers to the local `mocassin.2.02.73.2` source tree, not v3.

The [measured atomic comparison](phase1_77_mocassin_v2/README.md) records
individual ions and source hashes. The [invoked physical paths](mocassin_v2_physics_scope.md)
provide source locations; their historical 64-stage completion statements are
superseded by this objective. Table presence, routine existence, compilation,
and production invocation are distinct evidence. No matched nebular run has
been performed. This document does not establish equivalence of dust physics,
geometry, radiative transfer, or every output of the two complete applications.

## Present in v2 data or intended calculations, absent or incomplete in MoCHII

| Item | Verified gap | Required implementation |
| --- | --- | --- |
| Metal stage roster | Current production files declare 43 stages and 32 adjacent transitions; the target is 77 and 66. (Corrected 2026-10-11, audit Section 5 item 3: with 43 stages the terminal-band check refused **every input containing metals** at the default `eion_max` = 100 eV, since the highest carried stages had IPs of 47.22-64.49 eV; 155 of 175 tracked metal inputs and the Lexington benchmarks stopped, with no bypass; audit E1.) | Complete physical rate data and install the expanded chains; capacity alone does not activate them. State note: done 2026-10-11 (77 stages / 66 transitions; a band above the lowest terminal IP, Cl VII 114.2 eV, is refused). |
| Installed line atoms and fits | The selected reference directory contains 55 atoms with 571 levels and 49 primary cooling fits. Mg I has the approved five-level constant-strength approximation; Ca VII is absent. An available atom can calculate absolute cooling without a fit. | Resolve Ca VII physical data; activate upper stages and validate their consumers. Mg I retains an explicit model interval of 1000-100000 K; this interval does not certify atomic accuracy. All seven previously smaller shared-ion atoms now have the approved source coverage; Ca VII still requires correction. |
| Explicit level populations and line emission | The selected reference atoms now contain C III 10, N III 10, O I 7, Ne III 9, S IV 12, Fe II 142, and Si II 12 levels. Their matching thermal and radiation models were installed together. | Completed paired installation and the scoped native checks; full spatial validation remains outstanding. |
| Production oxygen activation | The selected reference model has 10 levels for each of O IV and O V. Explicit `extended` selection provides the installed 147-level O IV and 166-level O V models. Active oxygen still ends at O III. | Activate them with a complete oxygen rate chain; retain explicit source-accuracy qualifications. |
| Transition storage | The former 64-entry limit could not contain 66 transitions. | Changed to `MAX_METAL_ELEMENTS*(MAX_ION_STAGES-1)=66`; full build and staged native consumer checks passed. Production activation remains separate. |
| Bare terminal carbon | The former terminal-edge check demanded a nonzero ionization potential for fully stripped C VII. | Corrected to recognize that a bare ion cannot undergo further photoionization. Bound C VI must still reject an excessive band. |

The low-density cooling fits can include excitation and spontaneous cascades
outside the installed explicit atom. A smaller `NLEV` therefore does not prove
that all higher-level contributions are absent from the low-density total.
Assess that product's actual source path separately.

## Present more broadly in MoCHII than in the inspected v2 calculation

| Item | Source-path distinction | Input control and current status |
| --- | --- | --- |
| Additional ion line atoms | Installed Ne I, S I, Mg III, Fe III, Fe IV, Si I, and Si V atoms have no corresponding available line file in the inspected v2 mapping. | `par%additional_line_cooling=.false.` by default. `.true.` enables their thermal line contributions and line diagnostics. Same selection is used by setup and emission. (Corrected 2026-10-11, audit Section 5 item 1: `.true.` alone enabled none of them. These atoms existed only in `extended/`, so it also needed `metal_excitation_model='extended'`; with the default `'mocassin_v2'` the files were absent, and Fe III cooling and [Fe III] 4658 were gone in every Fe run, audit E4.) State note: the parameter is refused; every stage with an atom in the selected model cools and emits, and the default `'chianti'` model carries all seven. |
| Metal electron-impact ionization | v2 sets metal `collIon=0`; MoCHII has Voronov coefficients and element/stage corrections. | `par%metal_collisional_ionization=.false.` by default. `.true.` enables direct and cached population rates and the matching ionization-energy loss. State note: default `.true.` since 2026-10-11 (HEAD computed it always). |
| High-charge transfer approximation | The adopted q >= 5 unit-transfer Langevin upper limit is not invoked by the inspected v2 balance. | `par%high_charge_exchange=.false.` by default. Explicit `.true.` enables population, hydrogen charge, and mean-heat consumers. Both choices passed staged native validation. (Corrected 2026-10-11, audit Section 5 item 6: in this tree it had no effect, because no production element reached q >= 5 with 43 stages; only the staged files exercised it, audit E8.) State note: renamed `langevin_charge_exchange` (default `.false.`, 22 reactions at q = 5, 6, now reached by the 77-stage roster); the old name is refused. |
| Larger explicit atoms for shared ions | Examples include the installed O II/O III and prepared O IV/O V atoms. This is a broader representation of the same excitation process. | `par%metal_excitation_model='mocassin_v2'` is the default coherent selection. Explicit `'extended'` selects the larger models. Thermal cooling, level populations, line emission, and charge-transfer product populations read the same selected directory. The installed selection passed scoped native validation; Ca VII and 77-stage activation remain incomplete. State note: the default is now `'chianti'` (full CHIANTI 11.0.2, 70 atoms, Ca VII from `ca_7`); `'mocassin_v2'` is kept for parity; `'extended'` is refused. |
| Inner-shell contributions to metal photoionization rates | v2 opacity includes shells, but its cell rate and heat paths select the outer shell. MoCHII now selects the outer shell consistently in opacity, ionization rates, and heat by default. | `par%metal_inner_shell_absorption=.false.` by default. An explicit `.true.` request stops if an unsupported inner shell actually absorbs; a qualified vacancy-charge and electron-energy calculation remains unfinished. The control does not yet provide active inner-shell physics. (Corrected 2026-10-11, audit Section 5 item 2: the default was not a correction of the opacity. It deleted valid absorption above `einn` that v2 keeps in its opacity and HEAD summed (`phfit2`), and the `.true.` request stopped the run, flagging even shells inside the VFKY96 fit such as O I 2s at 28.5 eV; audit E2.) State note: the `phfit2` subshell sum is restored for absorption, rates and heat, and the parameter is refused; the Auger holes keep the flagged one-stage product until the stage matrix exists. |
| Neutral-hydrogen impact excitation | MoCHII supplies H-impact coefficients for C I, C II, and O I; the inspected v2 metal equilibrium forms electron collision rates. | `par%neutral_hydrogen_excitation=.false.` by default. `.true.` enables the joint-collider thermal increment and emitted-line calculation. Installed O I thermal/emitted-energy checks passed. State note: default `.true.` since 2026-10-11. |
| Charge-transfer kinetic energy | The inspected `update_mod::thermBalance` has no charge-transfer kinetic-energy contribution. | `par%charge_exchange_heating=.false.` by default. Explicit `.true.` includes metal-H and He-H terms without removing their ionization reactions. Native thermal selection checks passed. State note: default `.true.` since 2026-10-11. |
| Metal free-free cooling | The inspected v2 thermal balance includes H and He free-free losses only. | `par%metal_freefree=.true.` remains the explicit user-requested exception to the default-off rule. Explicit `.false.` disables its thermal contribution. Native thermal checks passed. |
| He I metastable population feedback | v2 includes triplet recombination emission, but lacks the invoked local 2^3S feedback on populations, gas energy, Penning reactions, and diffuse branches. | Existing `par%hei_metastable_coupling` now defaults to `.false.`; explicit `.true.` retains its coupled calculation. Actual implicit-input and coupled-physics tests passed. Ground-state He-H charge exchange remains shared physics with its default on. State note: default `.true.` again since 2026-10-11; `add_fuv` (from 4.78 eV) is on when the input leaves it out and the 2^3S population (with helium), a metal with a neutral edge below 13.598 eV, or dust is present. |

An additional cooling process is not necessarily negligible. Its contribution
depends on temperature, density, abundances, depletion, radiation, and the ionic
populations. No small-effect claim is made here. The input defaults implement
an explicit scope choice; they do not certify physical completeness for hot,
neutral, dense, or hard-spectrum plasmas.

The new controls are independent. Disabling line emission for an additional
ion does not remove that ion's charge, opacity, recombination, or other shared
processes. Disabling metal collisional ionization removes both its population
source and its associated ionization-energy sink; it retains recombination
kinetic-energy cooling and shared photoheating.

```fortran
&parameters
  par%metal_collisional_ionization = .true.,
  par%additional_line_cooling = .true.
  par%neutral_hydrogen_excitation = .true.
/
```

The parameters are read through the existing `/parameters/` namelist. Existing
inputs that omit them use `.false.`. Previously enabled extra-ion line cooling
and metal collisional ionization consequently change default behavior.

State note (2026-10-11): this example no longer runs. `additional_line_cooling`
is refused with its replacement named; the other two items default to
`.true.`.

## Errors in the preserved external reference

The inspected v2 Fe II packed radiation reader transposes the intended rate
orientation relative to its equilibrium consumer. Fe VI (19 declared levels)
and Ca VII (27) exceed the default 17-level output array and cause an explicit
stop when invoked. Ca VII also has four positive A values with nonpositive
energy differences in its supplied level ordering. The v2 hydrogen
collisional-ionization ratio lacks the electron-density factor, and two
Ca-labeled charge-transfer assignments overwrite Ar indices. These source
findings must not be copied into MoCHII. Preserve the external files for direct
comparison, implement physically oriented rates and energies, and identify
which corrected models have passed numerical and physical review.

(Corrected 2026-10-11, audit Section 5 item 9: the list above omitted six
v2 defects found by the audit (Sections 1 and 4).
- The collisional line cooling sums only levels 1-10 of each atom: the loop
  over `forbiddenLines(elem,ion,j,k)` runs `j, k = 1, 10`
  (`update_mod.f90:1101-1111`, read for this note).
- Fe II never enters the cooling: its atom is solved into
  `forbiddenLinesLarge` (`update_mod.f90:1279-1293`), which that sum does
  not read (read for this note).
- There is no dielectronic recombination for He.
- Elements without DR data are given the oxygen DR.
- The He+ recombination cooling uses the hydrogen coefficients.
- The free-free opacity loop fills only bin 1: for `i > 1`, `FFOpacity(i)`
  is set only to zero when `FFOpacity(i-1) < 1e-35` and otherwise keeps its
  initial zero, so only `i = 1` is computed (`ionization_mod.f90:369-392`,
  read for this note).
The He DR, oxygen-DR substitution and He+ recombination-cooling items are
taken from the audit; their v2 source lines were not reread here. None of
the six is copied into MoCHII: its line cooling uses every level of the
selected atom, Fe II included; its recombination is Badnell (Zhang et al.
2026 for Fe3+ and Fe V-VII) with each ion's own coefficient; its free-free
is evaluated at every frequency.)

Reading an Auger table is not an invoked Auger calculation. Likewise,
recombination-pumping code with zero/default coefficients does not prove
nonzero excited-level feeding. Neither is automatically a missing v2 process
that MoCHII must emulate.

## Validation and remaining work

The first implementation validation uses the actual namelist reader, linked
species and n-level modules, direct and cached stationary populations,
ionization-energy loss, and terminal photon-band entry point. It checks all
154 line-scope decisions (77 ions with each option state) against the measured
v2 file inventory. Synthetic carbon coefficients test numerical consumers,
not the accuracy of carbon atomic data. Physical-source validation remains
mandatory for production expansion.

Results and controller records are stored in
[the optional-physics validation directory](validation_2026-10-11/optional_physics/).
The initial build linked successfully but an existing verifier incorrectly
required objects for every source file, including a module outside the
production target. The verifier now queries the actual Makefile object list.
The first numerical fixture also placed its synthetic CI coefficient in the
wrong column; the failed products are preserved, and a corrected fixture uses
a separate directory. A completed controller is numerical evidence only.

Remaining work is to correct and install the Ca VII atom, activate the
prepared 77-stage rate chains, and validate their full production thermal
and radiation consumers. The seven shared-ion atomic replacements and
coherent excitation selector are installed. Default-off controls are
implemented, but an enabled inner-shell request remains explicitly rejected
at an active unsupported edge. A matched spatial comparison with MOCASSIN
v2 has not been performed. This objective is not complete.

## Additional implementation and verification on October 11

`neutral_hydrogen_excitation` now requires explicit input. Species setup omits
its retained thermal atoms when disabled, `metal_cooling_H` returns before
population/rate work, the thermal consumer skips the term, and metal line
output omits the neutral-H collider. The underlying scientific library routines
remain callable for explicit source diagnostics. The default change can be
important in neutral gas; enable the option when that physics is required.

The completed [controller validation](validation_2026-10-11/neutral_hydrogen_option/root_review.json)
used actual installed Ne I and O I files, 36 temperature/density cases,
absent/present namelist options, and the thermal and photon-energy consumers.
The O I joint-collider energy increment agrees with the corresponding emitted
energy difference. The earlier 154 ion-selection and direct/cached ionization
checks also passed against the rebuilt modules. No full production line-output
file or spatial comparison was generated.

The rate generator now includes Mg4+ + H -> Mg3+ + H+ using the final KF96
Table 1 coefficients `(6.36, 0.55, 3.86, -5.19)`, in units of
`10^-9 cm^3 s^-1`, over 1,000--30,000 K. KF99 Table 1 supplies the separate
8.60 eV mean kinetic heat. Both final tables were extracted and visually
inspected, and their methods retain the ground-reactant, low-density
radiative-relaxation assumptions. Endpoint rate continuation is an
approximation outside the fit interval; the mean heat is not a resolved
product-state distribution. The inspected v2 already has this Mg4+ rate.
The published-rate check and charge-transfer-table checks passed.

The adopted q>=5 unit-transfer Langevin upper limit is extended to the twelve
additional high-charge directions required by the 77-stage roster. It remains
an upper-limit approximation without a certified temperature-accuracy range,
not a measured rate. The independent capture and kinetic-energy tests passed.
[All eleven staged rate files](phase1_77_mocassin_v2/staged_rates/manifest.json)
now contain 77 stages and 66 adjacent transitions. They have not replaced
production atomic files. The S II atom is included unchanged as the dependency
of the sulfur product-term charge-transfer calculation.

The [actual 77-stage native check](validation_2026-10-11/staged77_charge_transfer_complete_inputs/README.md)
passed the 100 eV photon-band entry point and 198 reaction/temperature samples
covering direct/cached populations, charge conservation, and mean-heat
accounting. Root verified source and output hashes. An initial attempt lacked
the required S II atom and its failed record is preserved. These rate checks
do not establish line-atom completeness, source accuracy, or default-scope
parity of every reaction. Production still has the earlier 43-stage roster.

## October 11: default-off high-charge approximation

The user confirmed that additional physical options must be disabled by
default. The four implemented controls (`metal_collisional_ionization`,
`additional_line_cooling`, `neutral_hydrogen_excitation`, and
`high_charge_exchange`) all default to `.false.`. This statement does not
mean the unfinished controls for larger shared-ion atoms, inner-shell rates,
and other energy terms have already been implemented or audited.
(Corrected 2026-10-11, audit Section 5 items 5 and 6: the controls for
larger shared-ion atoms (`metal_excitation_model`) and for shell selection
(`metal_inner_shell_absorption`) did exist in this tree; and
`high_charge_exchange` had no effect in production, since no element reached
q >= 5 with 43 stages. State note: `metal_collisional_ionization` and
`neutral_hydrogen_excitation` default to `.true.` since the restoration;
`additional_line_cooling` and `high_charge_exchange` are refused, the latter
renamed `langevin_charge_exchange`, default `.false.`.)

The [high-charge validation](validation_2026-10-11/high_charge_option/README.md)
passed both actual namelist choices with 198 samples each, including 66
high-charge samples and preserved low-charge coefficients. The controller
completed; root verified its frozen sources and artifacts. Production remains
at 43 metal stages and 32 adjacent transitions.

## October 11: explicit metal free-free exception

The latest user instruction keeps `metal_freefree=.true.` by default.
`charge_exchange_heating` instead defaults to `.false.`. These choices are
independent of charge transfer in the ionization balance. The
[scoped thermal validation](validation_2026-10-11/thermal_defaults_complete_inputs/README.md)
passed 108 native cases and independently checked the net-energy differences.
This supersedes any general statement that all additional physical terms
are disabled by default. Remaining unimplemented scope controls are still
listed as incomplete.
State note (2026-10-11): `charge_exchange_heating` defaults to `.true.`
again (HEAD's default).

## October 11: native finite-level source validation

54 coherent supplied v2 atoms now have candidate native files and completed
numerical consumer checks, including corrected Fe II and complete 19-level
Fe VI. Mg I and Ca VII remain excluded pending physical source review.
This changes preparation status, not the installed 43-stage production
roster or missing-installed-model counts quoted above. The source-table
interpolation is intentionally corrected to remain positive. Original
N I/Ne V cubic expressions produce negative strengths in sampled intervals;
individual Fe II strength differences can also be large. Neither exact
v2 numerical equality nor negligible cooling differences are claimed. See
[the measured results and matched excitation inputs](validation_2026-10-11/mocassin_source_atoms/README.md).

## October 11: complete curve review and Fe II qualification

All 54 comparison pages were visually inspected and their source/product
hashes verified. The independent strict low-density calculation differs
from the native density-floor column by at most 2.44430287891e-7 relative
over the sampled curves; this is a finite-density comparison, not a
physical collision-strength accuracy bound.

The original Fe II file contains allowed transitions only. Its raw A table
leaves 35 levels radiatively closed and does not include the ground
fine-structure probabilities present in the installed MoCHII 13-level atom.
Correcting its packed reader direction does not complete forbidden-line
physics or the omitted allowed transitions documented in its own header.
Model selection must state this restriction explicitly. The
[source-profile report](validation_2026-10-11/mocassin_source_atoms/README.md)
records the distinction. Cooling fits remain candidate products until
compiled-reader, temperature-range, and production-integration checks pass.

## October 11: cooling coefficient preparation checkpoint

52 supplied-model cooling fits passed the actual compiled reader and
evaluator on 801 temperatures each. They remain candidate products, with
required temperature bounds and production integration unfinished. Ne VI
and Fe V require a better numerical representation; a small partial fit
error does not establish acceptable conditioning. The
[compiled-reader report](validation_2026-10-11/mocassin_cooling_native_exact_endpoints/README.md)
records the measured errors and preserved failed attempts. No installed
coverage counts or active stage counts changed.


## Absolute cooling integration and Fe II diagnostic

Native atoms without exponential fits are now connected to the production
cooling consumer. The focused actual linked Ne VI/Fe V tests passed 528 samples
and eight zero-density cases. Maximum matrix comparison error was 0.124251%.
The complete two-page curves were visually inspected; no expanded production
atom installation or spatial calculation is claimed. See
[the report](validation_2026-10-11/atomic_cooling_without_fit/README.md).

The assumption that every additional MoCHII cooling term is negligible fails
for the measured Fe II excitation diagnostic. Adding the 40 installed
forbidden A values to otherwise identical 142-level source data increases
Fe II cooling by 5.33244 at 10,000 K and electron density 100 cm^-3. This is
not total nebular cooling. The matched-input table, full temperature curves,
and limitations are in [the report](validation_2026-10-11/feii_forbidden_effect/README.md).
Default selection requires the user's physical choice; current installed
Fe II was not changed. `metal_freefree=.true.` remains the explicit exception.


## Mg I cold source preparation

The single-temperature MOCASSIN Mg I data now have a restricted five-level
native candidate under the adopted cold constant-strength approximation.
The actual reader, level solution, line emission, and cooling passed 404
samples over 1,000--10,000 K; a 10,000.1 K evaluation correctly failed.
Maximum native/independent cooling difference was 2.18603e-13. Full cooling
and five-level population curves at all four densities were visually reviewed.
See [the matched-input report](validation_2026-10-11/mgi_cold_source/README.md).
This candidate is not installed, has no verified high-temperature accuracy,
and does not complete 1,000--100,000 K Mg I coverage. The previous 54-candidate
inventory remains unchanged; this restricted candidate is separate.


## Ca VII energy-basis investigation

The measured source/current-CHIANTI comparison identifies mixed observed and
theoretical energy ordering as the numerical cause of four negative gaps.
All-theoretical energies give positive gaps for all 310 original positive A
values, but that physical replacement is not yet accepted. The current
CHIANTI table has only 86 of those pairs, with identical stored A values.
Simply using it would omit 224 original pairs. See
[the data and article-access report](validation_2026-10-11/cavii_energy_basis/README.md).
Final published methods and transition definitions remain unverified; no
external source or production Ca VII atom was changed.


## Explicit He I metastable coupling

The existing `hei_metastable_coupling` default is now `.false.`. The production
build and seven actual implicit-input cases passed, as did the explicitly
enabled population, energy, and diffuse checks. Default inputs also stop
implicitly enabling the FUV extension for this omitted population. Ordinary
He I recombination cascades and ground-state He-H charge exchange remain
available; this selection concerns the additional excited-population feedback.
See [the source-path and validation report](validation_2026-10-11/hei_metastable_default/README.md).
Earlier runs using an omitted option are historical outputs of the older
default, not verification of the present default. No spatial products were
regenerated.
State note (2026-10-11): `hei_metastable_coupling` defaults to `.true.`
again, and `add_fuv` (from 4.78 eV) follows the 2^3S population (with
helium), a metal with a neutral edge below 13.598 eV, or dust.


## Missing-atom installation checkpoint

23 supplied finite-level atoms are now installed with exact validated bytes
and no overwrites. Cl V now contributes absolute cooling within the existing
production chlorine chain even without a fit; its 264 scoped samples and
four zero-density checks passed, with maximum matrix difference 0.00600158%.
All installed paths were opened by the actual compiled reader. The new
measured missing-atom count is two, while 25 exponential fits remain absent
and seven explicit installed models remain smaller. Production still has
43 stages and 32 transitions. See
[the installation and matched-input report](validation_2026-10-11/missing_atom_installation/README.md).
The earlier candidate-only statements describe previous checkpoints and are
superseded for these 23 files; they do not imply 77-stage activation.


## Consistent shell-selection component

The private outer-shell photoionization and heat calculation passed 1,764
coefficient samples across the current and staged rate chains. It is invoked
below the existing PHINNER boundary with the former values preserved; an
outer-only default has not yet been selected above that boundary. No new
public input was added. See
[the code, matched coefficient inputs, and validation](validation_2026-10-11/outer_shell_photoionization/README.md).

A further v2 audit confirms that opacity loops over subshells while cell rate
and heat integrals select only the stored outer-shell window. The outer
4s window also depends on 3d/3p lower pointers and the photon-bin grid.
The corrected MoCHII model must retain valid outer tails and select the same
absorption scope in opacity, rates, and heat. The earlier proposal to keep
full shared opacity while removing all inner-vacancy rates is superseded by
this consistency requirement. Full inner-hole products and energies remain
unfinished.
(Corrected 2026-10-11, audit Section 5 item 2: the "valid outer tails"
requirement was not met by the outer-only default installed afterwards; it
removed the valid absorption above `einn` (audit E2). State note: the
`phfit2` subshell sum is restored and `outer_shell_photoionization` removed;
the channel-6 ground continua use the outer-shell part of the cross section
through `species_outer_shell_cross_section`.)

## October 11: original v2 temperature endpoint behavior

The original v2 archive confirms endpoint-constant collision strengths in
`interpolation_mod::splint`, as invoked for three or more temperature nodes.
The single-node branch is constant; the two-node branch is linear in log10(T)
without its own endpoint clamp. This distinguishes v2 numerical use from the
current candidate upper-temperature rejection. The four supplied tables ending
below 100000 K are C I, S III, Ar VI, and Mg I. No hot continuation has been
accepted or installed by this source inspection. See
[the original-source review](validation_2026-10-11/v2_collision_temperature_endpoints/README.md).

The registered-cooling diagnostic was corrected to stay within finite source
bounds and to stop reporting an unmeasured fixed table size. Eighty actual
C I/C III samples passed. Production physics and input defaults did not change.
See [the diagnostic correction](validation_2026-10-11/registered_temperature_scope/README.md).

## October 11: qualified hot endpoint candidates

The original v2 endpoint assumption is now prepared for C I, S III, and Ar VI
in isolated modules through 100000 K. Direct/native matrices passed 2412
samples; 1212 original-range profile comparisons are identical, and 120
registered cooling checks have maximum discrepancy 0.159297%. Nine bound
probes passed. Appending a synthetic node was rejected as a representation
choice because it changes previously accepted in-range PCHIP strengths.
Instead the private evaluator retains original source polynomials and clamps
only upper evaluations. The corrected continuation completed and all three
figure pages were inspected. No public input, production atom, hot physical
acceptance, or atomic-accuracy claim follows from this preparation. See
[the full review](validation_2026-10-11/upper_temperature_models_checked/README.md).
