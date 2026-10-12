# Current goal status: blocked after three consecutive source audits

The October 11 documentation update is complete separately from the physics
goal. See [the implementation results](../ATOMIC_IONIZATION_STAGE_EXPANSION_RESULTS_2026-10-11.md)
and [the documentation review](validation_2026-10-11/documentation_update/README.md).

The goal status was set to `blocked` after three consecutive resumed turns
verified the same unresolved Ca VII physical source condition. No new candidate
article or accepted corrected atom has appeared. Production remains 43 metal
stages and 32 adjacent transitions, with the original target of 77 and 66
preserved. This is not completion or a user-requested pause. No calculation
controller is running. Continuation requires complete qualified radiation
components/energies or an explicit material approximation decision.
[The source audit](validation_2026-10-11/cavii_completion_research/resumed_blocked_audit.json)
records the verified state. Earlier active-status checkpoints below are historical.

# Current state after the October 11 approval

The prior approval and article-access dependencies are resolved. The approved
six paired replacements, Fe II paired replacement, coherent model selector,
and finite hot collision-strength model have been installed. The default-off
inner-shell guard and selected-model consumers passed their scoped native
checks. The scoped controllers are terminal completed:
`mochii_approved_shell_pairs_checked_20261011` and
`mochii_approved_excitation_20261011`, and
`mochii_77_selected_model_20261011`. The final H I common-reader correction
and parser checks also passed the final production build without regenerating
the H I numerical table. There is no running calculation or
configured model-event delivery service at this checkpoint.

The prepared 77-stage chains now also pass 396 actual charge-transfer samples
with the selected reference model. Next: resolve Ca VII radiation and energy
inconsistencies together, then activate
and validate the complete rate/line-cooling roster. The global objective is
incomplete. Source uncertainty must not be presented as missing user approval.
See [the current report](validation_2026-10-11/approved_excitation_models/README.md) and
[the comparison](MOCASSIN_V2_MOCHII_BIDIRECTIONAL_COMPARISON.md).


After goal resumption, an additional source investigation completed: the ADS
API query succeeded, two later Ca VII radiation papers were identified, and
the official NORAD dipole-only dataset was preserved in references. Five
local Cloudy Stout versions have 26 levels, 15 distinct radiation pairs, and
no collision records, so they cannot replace the full source model. The later
papers' final PDFs remain inaccessible through the checked publisher and ADS
routes; their methods and suitability are unverified. See
[the additional evidence](validation_2026-10-11/cavii_completion_research/README.md#additional-source-investigation-after-the-goal-resumed).
No calculation controller is running. The previous goal turn only reported
status; this resumed turn added source evidence and corrected current-state
documentation. Ca VII and production 77-stage activation remain unfinished.

A subsequent resumed turn completed a read-only selected-stage audit and
corrected a confirmed omission: `extended` now retains the identical approved
Mg I atom. Actual input selection and production readers loaded 62 atoms with
1951 levels. No numerical products were regenerated. Both selections now
lack only Ca VII among supplied source models; production declarations remain
43 metal stages and 32 adjacent transitions. The requested later Ca VII
articles have not been added to references at this checkpoint. See
[the correction](validation_2026-10-11/extended_mgi_coverage/README.md) and
[the current coverage](validation_2026-10-11/selected_stage_coverage/README.md).
No calculation controller is running. This turn made implementation progress;
the full goal remains active and incomplete.

A further source diagnostic checked all 310 Ca VII radiation pairs under
conditional fixed-strength E1/M1/E2/M2 rescaling. Selection rules leave 111
single-component and 199 mixed-component pairs. The conditional frequency
factors range from 0.00153375758918 to 3.56296319229 across individual pairs;
these are not total-cooling errors or atomic-accuracy bounds. No corrected
atom was generated or accepted, and no calculation is running. The later
article PDFs remain absent from references. This turn added physical source
evidence; accepting Ca VII still depends on component data or a justified
new approximation. See
[the conditional assessment](validation_2026-10-11/cavii_completion_research/README.md#conditional-radiation-frequency-correction-assessment).

The next resumed audit verified that no later Ca VII candidate article was
added to references and no accepted Ca VII atom exists in either selected
directory. Production remains 43 stages and 32 transitions. There is no live
calculation to wait for. The unresolved condition is physical source data or
an explicit material approximation, not missing approval for the established
implementation. This is the first consecutive blocked audit after the recent
progress turns; the goal remains active. Recording this checkpoint is not
implementation progress. See
[the resumed audit](validation_2026-10-11/cavii_completion_research/resumed_blocked_audit.json).

The sections below are historical evidence; their old missing-article and
pending-approval statements are superseded by this current checkpoint.

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

# Historical checkpoint before receipt of the Ca VII final article

The last three consecutive goal turns made no implementation progress because
the same approval and source-access dependencies remained unresolved. The
latest automatic continuation is not an approval response. Independent
numerical preparation has completed; do not repeat completed calculations to
substitute for those decisions.

Direct inspection still measures 43 production metal stages and 32 adjacent
transitions. The staged target is 77 stages and 66 transitions. The proposed
`metal_excitation_model` and `metal_inner_shell_absorption` inputs are absent.
`metal_freefree=.true.` remains the explicit user exception. Process inspection
found no running calculation controller or Fe II/S II diagnostic executable.
The Fe II finite-density, sulfur charge-transfer model, and upper-temperature
candidate controllers are terminal completed.

Resume the dependent implementation after the user's response to the concrete
five-item proposal: paired replacement of six atom/fit sets; coherent reference
and extended excitation selection; default-off inner-shell absorption with
unsupported vacancy regimes rejected; qualified hot-endpoint continuation for
Mg I, C I, S III, and Ar VI; and retention of the 40 existing Fe II forbidden
probabilities as a physical default contribution. See the
[paired replacement review](SHARED_ATOM_EXPANSION_REVIEW_2026-10-11.md),
[excitation-model design](SHARED_EXCITATION_MODEL_SELECTION_REVIEW_2026-10-11.md),
and [shell-selection design](OUTER_SHELL_INPUT_DESIGN_REVIEW_2026-10-11.md).

Automatic approval review rejected the earlier atomic-data overwrite for
insufficient explicit approval. Do not bypass that rejection. Public-input and
data changes also require prior confirmation under the global AGENTS.md.

The local references directory contains 37 files and does not contain the
requested Landi and Bhatia (2003), ApJ 589, 1075-1084, DOI 10.1086/374686 final
article. Its methods are required for accepting corrected Ca VII energies and
radiation probabilities; see [the source-access record](validation_2026-10-11/cavii_energy_basis/README.md).
The goal is incomplete. Its blocked status records these dependencies, not
successful completion or an intentional pause.

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

# Qualified upper-temperature model preparation

C I, S III, and Ar VI upper-endpoint candidates passed 2412 combined direct
samples, 1212 exact in-range profile comparisons, 120 registered cooling
samples, and nine guard probes. All three full-profile pages were visually
inspected. Maximum registered/direct cooling difference was 0.159297%.
The original failed controller is preserved; its corrected continuation
completed. No production data or input changed, and hot atomic accuracy is
unverified. See [the candidate review](validation_2026-10-11/upper_temperature_models_checked/README.md).

Appending a synthetic hot node would change original in-range strengths by
up to 2.829553% (C I) and 1.271571% (S III). The candidate instead retains the
original source polynomials and clamps only upper-endpoint evaluations.
This preparation does not adopt the hot physical model or complete the global
77-stage production objective. Previous approval and physical choices remain
unanswered.

# Independent source-model audit and diagnostic correction

A fresh inventory measured twelve installed shared-ion atoms larger than the
supplied v2 model, ten in active rate chains. Their general input selection
remains unimplemented. The proposed coherent source-model selection also
includes form-9 charge-transfer product-term populations, not only thermal
cooling and lines. See [the design review](SHARED_EXCITATION_MODEL_SELECTION_REVIEW_2026-10-11.md).

The original v2 archive confirms constant endpoint collision strengths in its
multi-node spline path. Current source candidates' upper-temperature rejection
differs for C I, S III, Ar VI, and Mg I. No hot approximation or public input
was adopted. The registered-cooling diagnostic now respects finite source
bounds and stops claiming an unmeasured fixed table size; 80 C I/C III samples
passed. See [the correction](validation_2026-10-11/registered_temperature_scope/README.md).
All calculation controllers from this sequence are terminal completed. The
previously requested production replacement and Mg I hot-model decisions are
still unanswered; they are not inferred from automatic continuation.

# Current calculation state

The corrected PCHIP, Mg I constant-strength candidate, and matching-source
candidate controllers are terminal completed. All recorded sources and
artifacts were verified by root. No production atomic replacement is running. The isolated paired
thermal/emission controller `mochii_shared_matched_native_20261011` completed
1584 thermal and 1584 emitted-energy samples, 24 zero-density checks, and
six comparison pages. The O I hydrogen-excitation controller completed
18 samples. Root verified both results and visually inspected all six pages.
All these controllers are terminal; there is no calculation running in this
sequence. See [the actual combined review](validation_2026-10-11/shared_atom_matched_native/README.md).
No model completion-event watcher is configured. The 801-node
continuation remains failed after its source-change check, and the first
matching-source preparation remains failed after selecting an obsolete partial
fit record. Their products are preserved, not silently restarted or relabeled.

The next dependent actions are explicit twelve-file paired replacement
approval, Mg I hot-model acceptance, and the previously requested Fe II and
shell-selector choices. Independent preparation and correction may continue.
See [the paired review](SHARED_ATOM_EXPANSION_REVIEW_2026-10-11.md) and
[Mg I report](validation_2026-10-11/mgi_constant_source/README.md).

# Latest independent source work

The read-only KM Eq. (7) audit completed for 14 source subshells. It corrected
an overly restrictive joint-energy requirement for the current linear metal
heat consumer and rejected unconditional use of the published mean correction
because Ca II M1 gives -5.022 eV. No calculation controller is running for
this short source review. The public shell selector and Fe II physical-default
choice still await the previously requested answers; they were not inferred
from the `metal_freefree=.true.` instruction. Ca VII publisher metadata now
confirms pages 1075-1084; final methods remain inaccessible.
See [the source review](validation_2026-10-11/auger_mean_energy/README.md).

# Public shell-selector design awaiting confirmation

The tested outer-shell component is ready for a consistent opacity/rate/heat
selector. The proposed `metal_inner_shell_absorption=.false.` default and
rejection of unsupported vacancy regimes are recorded in
[the reviewable input design](OUTER_SHELL_INPUT_DESIGN_REVIEW_2026-10-11.md).
The global instruction requires confirmation before changing public inputs.
No such parameter was added. Independent source work remains possible;
this is not a claim that the entire goal is blocked or completed.

# Latest outer-shell component checkpoint

A private consistent outer-shell cross-section and heat component passed
1,764 coefficient checks in current and staged rate chains. Production
below-PHINNER arithmetic is retained; no above-PHINNER default or public input
selector has been changed. The full build and diagnostic controller completed;
root verified sources and outputs. The existing inner-vacancy product/energy
limitations remain uncorrected in the above-edge production calculation.
See [the report](validation_2026-10-11/outer_shell_photoionization/README.md).
The full objective remains incomplete.

# Latest installed-atom checkpoint

23 absent supplied finite-level atoms are installed without overwrites. All
new paths opened successfully; the active Cl V no-fit cooling term passed
264 samples and four zero-density checks. Root verified hashes and visually
inspected the full curves. The controller completed both recorded steps and
is no longer running. Source accuracy is qualified; upper-stage invocation
remains unfinished. Production still has 43 metal stages and 32 transitions.
See [the report](validation_2026-10-11/missing_atom_installation/README.md).

# Latest He I default checkpoint

The existing `hei_metastable_coupling` now defaults to `.false.`. The full
production debug build, seven actual implicit-input cases, and existing
explicitly enabled coupled-physics checks completed and passed. Root verified
source and artifact hashes. No controller in that sequence remains live.
See [the report](validation_2026-10-11/hei_metastable_default/README.md).
This continues the requested opt-in selection, not full goal completion.

# Latest verified continuation checkpoint

The no-fit atomic cooling, Fe II forbidden-transition, and restricted Mg I
cold-excitation controllers completed all recorded commands and checks. Their output hashes, numerical profiles,
and figures have been reviewed. No calculation from those sequences is live.
See [no-fit cooling](validation_2026-10-11/atomic_cooling_without_fit/README.md),
[Fe II](validation_2026-10-11/feii_forbidden_effect/README.md), and
[Mg I](validation_2026-10-11/mgi_cold_source/README.md).
Production activation remains unfinished. Fe II default selection is awaiting
user clarification; independent source and integration work can continue.

# Current continuation status

The main session is continuing the scientific work directly. The separate
model watcher remains stopped after its execution tools failed; automatic
model reasoning or a completion notification in this chat is not claimed.
Finite calculation controllers execute only their recorded commands and checks.

Production still declares 43 metal stages and 32 adjacent transitions. The
current 77-stage, 66-transition first-phase roster remains staged, and the
overall work is not complete. The former 64-stage target is superseded. Source-node storage and dynamic consumers are implemented.
Installed source/native-verified corrections include eleven low-density fits,
S I, full S III, Ca IV ground fine structure, full Fe IV, and full C I (42
levels), C II (53), N I (50), N II (58), and S II (70). These five correct unjustified
unsuppressed residual cooling and missing high-density bound excitation.

The absolute registered-cooling correction is installed and its actual build
and positive-subnormal/zero checks passed. All earlier remaining-source
registration checks completed. All four full-bound Intel/GNU source profiles
passed; all four GNU full-bound tables passed. The duplicate Intel Ne V table
was deliberately interrupted after the equivalent GNU calculation passed;
its saved state is an interruption, not a completed step. The current
C V, N VI, and Ar VI registration controller also completed all its steps.

Independent work continues on remaining low-ion density dependence, final
source temperature intervals and adopted cold Upsilon continuation, missing
cooling sources, terminal reverse rates, and nonadjacent Auger energy
accounting. A stationary nonadjacent solver has been tested independently but
is not linked or invoked by production. Production roster activation and
matched full nebular comparisons remain unfinished.

See [implementation status](implementation_status.md),
[compact-atom correction](existing_compact_atom_correction.md),
[subnormal correction](subnormal_excitation_cooling_correction.md), and
[extended source validation](excitation_temperature_range_validation.md).
No overall completion or expanded production acceptance is claimed.

C IV (24 states) and Si V (169) cooling atoms and fits are installed, with
35 fits among the 43 active metal stages after corrected Ne I installation. The actual production BT-floor
build and focused native validation completed. The existing O III optical
A values are corrected to final published totals and passed native emission
and cooling checks. Source-limited registry grids are installed and passed actual rebuilt Intel
validation, including invalid-input and subnormal/zero checks. Ar II ground
fine-structure cooling is installed. Investigation continues on Cl V and
coherent complete Ca I excitation sources.

The corrected Ne I 31/144 model is installed after all seven outer
controller steps and all three private product steps completed. Root checked
43 outer and 46 product frozen hashes, full cooling and all-level population
curves, and the fit curve. Native/source cooling differs by at most 2.22e-15;
fit error is 0.0317686%, table application error 0.0041481%. Installed files
preserve exact validated bytes and passed a production debug open check.
The earlier failed rounding-identity pilot and the physically stale 31/142
model remain preserved. Neither Ne I controller is running now. Mg I remains
unaccepted. Additional Si I, Fe I, and Ar I source investigations completed without
altering the approved public roster or data layout. Their unresolved inputs
and the remaining ionization/thermal requirements are recorded in
[the remaining physical requirements](remaining_physical_requirements_2026-10-10.md). See
[the Ne I/Mg I audit](source_validity_ne_i_mg_i.md).

The private oxygen radiative controller `oxygen_radiative_v2_20261010`
completed both steps (controller PID 4861 at launch). Its
[durable state](validation_2026-10-10/oxygen_valence_convergence/controller_state_v2/state.json)
records source and executable hashes, process identifiers, exit statuses, and
artifacts. Candidate construction passed before the native successor started.
The native step computed 401 temperatures from 1,000 to 100,000 K, four
densities, and three radiative-rate models for each of O II and O III, with
finite memory, CPU, and wall-time limits. This originally private physical-model
diagnostic was followed by the reviewed oxygen installation described below;
it is not a nebular comparison.
An earlier 4 GiB controller failed during OpenBLAS startup; its failure record
is preserved separately. The current controller limits BLAS to one thread and
uses an 8 GiB limit. No uncertain partial calculation is treated as completed.

The separate tmux controller `oxygen_independent_20261010` (PID 6048 at
launch) automatically advanced after that native controller completed.
It completed original collision-node preservation, effective strengths,
independent Python populations and cooling, level balance, and collision-to-
radiation energy closure. Its
[state](validation_2026-10-10/oxygen_valence_convergence/successor_state/state.json)
and
[commands](validation_2026-10-10/oxygen_valence_convergence/successor_plan.json)
persist the authorized successor. Neither controller installs an atom or
claims physical acceptance from a zero process exit status.

Main-session review verified both completed oxygen controllers' frozen
source and artifact hashes and inspected the full profile PDFs. Independent
native cooling errors were `2.9643e-14` for O II and `2.0428e-14` for O III.
The conditional unidentified same-parity rate tests are not atomic-accuracy
bounds. Generator integration and full-atom installation subsequently completed.
The isolated Intel optimized controller `oxygen_registered_20261010` completed
both full registered tables. Main review checked all 21 build-manifest hashes
and both results; maximum application errors are 0.197768% and 0.143563%.
Low-density coefficients were preserved. No oxygen calculation is pending.

## Preserved Mg III execution history (superseded by installation)

Mg III private native and independent profiles completed over 801 temperatures
and four densities. A final artifact-name mismatch left the original recovery
controller marked failed; its numerical comparison passed and a separate
verification-only audit checked existing source/product hashes. No numerical
rerun was performed. At that point, collision/radiation energy closure and guarded source
integration remained to be checked. The later corrected installation below
supersedes that pending status.

The subsequent official NIST identity check found a misassigned Mg III
observed energy at level 167 (J=0 source label versus J=1 observed level).
The source selection now requires a guarded correction to that state's own
theoretical energy, 634430 cm^-1. The old native profiles are preserved as
physically stale conditional diagnostics. A new finite workflow was then
prepared for the corrected energy and its downstream A, populations, and
cooling; no production fit or installation is authorized before validation.

The corrected Mg III sequence subsequently ran in persistent tmux session
`mochii_mgiii_corrected_20261010`, controller PID 25986 at launch. Its
[state](validation_2026-10-10/mgiii_corrected_level_identity/controller/state.json)
freezes 40 source/input hashes. Pilot validation completed before the native
full167 successor. The queued successors are independent full-atom balance,
comparison, collision/radiation ledger, finite result manifest, and
synchronous private atom/fit/native-coefficient/registered-table validation.
The physical ID167 correction and guarded E1 model were reviewed before
launch. No installation command is in this sequence. Main-session review of
completed results and complete profile plots precedes installation.

Both corrected Mg III controllers completed all guarded steps. Main review
verified 40 and 44 frozen source/input hashes before subsequent header-only
source edits, and inspected the complete source and fit PDFs. The corrected
full167 atom and fit are installed with unchanged numerical records and
explicit ID167 source-error annotations. The current production-module
driver opened the installed atom and matched the validated four-density
10,000 K profile exactly. That installation brought production to 34 electron-excitation fits.
At that checkpoint, the later Ne I installation brought it to 35 fits
among 43 metal stages and 32 adjacent transitions. Both sequences completed;
no overall completion is claimed.


## End of the completed Ne I and source-identity sequence

At the end of the Ne I checkpoint, no calculation controller remained running
from the Ne I, Mg III, or oxygen sequences. The preliminary Si I, Fe I,
and Ar I audits at that time did not yet support another installation. The
later original NIST Si I recovery and installation supersede that checkpoint. Root rechecked the relevant raw rows and
source counts. At that checkpoint, production had 43 stages, 32 adjacent transitions,
and 35 ordinary metal cooling fits. The remaining source inputs and material
physical closure requirements are listed in
[remaining physical requirements](remaining_physical_requirements_2026-10-10.md).
This is a completed sequence within an incomplete overall expansion;
no automatic model-resume watcher or chat completion delivery is configured.

## Latest scope decision and short diagnostics

The user retained 64 metal stages and 53 adjacent transitions. The 165-stage
roster is excluded from this task. The original-node terminal reverse-rate
diagnostic completed locally, and root verified its nine recorded input hashes.
Its one-added-state inequality applies only to the restricted adjacent,
electron-only stationary model; sparse common source nodes do not cover the
whole temperature interval. Root also checked the Auger diagnostic's three
input hashes, probability-column conservation, and ground-state event-energy
bounds. These completed diagnostics do not activate either model in production.

The Si I/Mg I and Fe I/Ca I source investigations have completed. Some
individual atomic data and final-publication methods remain unavailable.
There is no new long calculation controller or model-resume watcher running.
The outstanding physical conditions remain in the
[requirements report](remaining_physical_requirements_2026-10-10.md).

## Completed Si I numerical validation and installation

The persistent eight-step Si I controller completed successfully. Its
[durable state](validation_2026-10-10/si_i_radiative_completion/controller/state.json)
records every completed command, child PID, exit status, artifact hash, and
63 frozen source/input hashes. Root verified these hashes and full PDFs,
installed the coherent atom and fit, and exercised the actual installed
file with the current native driver. See the
[installation review](validation_2026-10-10/si_i_cooling_products/INSTALLATION_REVIEW.md).
No controller is running from this completed Si I sequence; its launch
session was `mochii_sii_complete_20261010` on socket `mochii-workflows`, PID
65855. No automatic model-resume or UI completion notification is configured.

Production remains at 43 stages and 32 transitions, now with 36 ordinary
metal electron-excitation fits. Seven excitation sources and the Auger
joint charge/energy and terminal-stage validity requirements remain unresolved.
The overall 64-stage expansion is not complete.

## Completed O IV full-bound conditional diagnostics

The five-step O IV controller completed at 2026-10-10 14:05:08 UTC. Its tmux session was `mochii_oiv_full_20261010` on socket `mochii-workflows`, controller PID 6069. Root verified all 53 frozen source hashes, all five completed-step artifact hashes, the saved native product hashes, and all seven PDF pages. See the [durable state](validation_2026-10-10/o_iv_full_completion/controller/state.json) and [root review](validation_2026-10-10/o_iv_full_completion/root_completion_review.json). This controller is no longer running.

The 147-level atom passed original-node serialization, source/native populations and cooling at 401 temperatures and four densities, stationary energy balance, an 801-temperature cooling fit, the actual GNU coefficient reader, and the Intel registered cooling calculation. Maximum source/native cooling discrepancy is 1.39888101103e-14; normalized energy discrepancy is 4.37526002907e-14; fit error is 0.287515276%; registered application error is 0.277913875%. The published 3,200 K cold boundary uses the user's adopted constant-Upsilon approximation below that temperature, with Maxwellian factors evaluated at the actual gas temperature.

These are conditional numerical diagnostics: `physical_model_accepted: false`, `production_installed: false`. Root independently reproduced the original radiative-source audit. Most of the higher gf/A inconsistencies are also present when CHIANTI A agrees with the original published ADF04 rate; they do not prove that A is wrong and do not justify wavelength-based rescaling. The source correspondence investigation continues, and a separate O V full-bound sequence is being prepared. No automatic model-resume or UI completion event delivery is configured. Production still has 43 metal stages and 36 ordinary metal cooling fits; the 64-stage expansion is incomplete.

## Completed O V full-bound conditional diagnostics

The O V controller completed all five steps at 2026-10-10 14:16:31 UTC. Its session was `mochii_ov_full_20261010`, socket `mochii-workflows`, controller PID 11329. Root verified all 58 frozen input hashes, all step artifacts and saved native products, and all seven PDF pages. See [durable state](validation_2026-10-10/o_v_full_completion/controller/state.json) and [root review](validation_2026-10-10/o_v_full_completion/root_completion_review.json). No calculation from this sequence remains running.

The 166-level source/native cooling discrepancy is 2.02060590482e-14; normalized energy residual is 6.63342718617e-14; maximum fit error is 0.139092942%; registered application error is 0.099054140%. O V remains a conditional, uninstalled source model. Its numerical completion is distinct from physical adoption and overall 64-stage completion.

Root subsequently adopted and installed the qualified O IV source-data cooling model after verifying the actual A/energy consumer path, preserving the previous conditional diagnostics and executed generator sources. The [O IV installation review](validation_2026-10-10/o_iv_full_completion/installation/README.md) records the installed reader checks and source limitations. Production still has 43 metal stages/32 transitions; adding the O IV data files does not activate the fourth oxygen stage. No background calculation controller or automatic model-resume/UI completion notification is currently running.

## Latest installed oxygen data

Both O IV and O V full-bound source-data models are now installed after qualified root adoption and actual installed native/fit-reader checks. The original conditional controller results remain preserved. Source atomic accuracy is still unverified. O V installation supersedes the preceding conditional/uninstalled checkpoint; see [the installed model](validation_2026-10-10/o_v_full_completion/installation/README.md). No controller is running. The production ELEMENT roster remains 43 metals/32 transitions; the 64-stage objective is not complete.
