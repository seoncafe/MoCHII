# Current production status (October 11, 2026)

The approved atomic expansion remains incomplete. Production declares 43 metal
stages and 32 adjacent transitions; 77 stages and 66 transitions are prepared
but not active. The default `mocassin_v2` selection contains 55 atoms (571
levels) and 49 cooling fits (763 terms); `extended` contains 62 atoms (1,951
levels) and 38 fits (377 terms). Both include the approved five-level Mg I
constant-strength model over 1,000--100,000 K. Six shared atom/cooling pairs
and the 142-level Fe II atom with 1,184 radiative transitions, including 40
preserved forbidden transitions, are installed. Only Ca VII lacks an accepted
reference atom. Its 27-level source has four nonpositive photon gaps and
mixed five- and twenty-five-configuration radiation probabilities; numerical
conversion alone cannot establish physical accuracy.

The earlier approval-dependent statements in this cumulative record describe
their dated checkpoints. The approval was received and applied where the
physical model was qualified. `metal_freefree=.true.` remains the explicit
default exception. See the [current result](../ATOMIC_IONIZATION_STAGE_EXPANSION_RESULTS_2026-10-11.md)
and [Ca VII audit](validation_2026-10-11/cavii_published_review_checked/README.md).
No full spatial comparison or full thermal-equilibrium run has been completed.

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
The current 177-level O III optical energy ratio was also remeasured directly
in twelve samples: 2.9839683327599866--2.983968332759988. See
[the current measurement](validation_2026-10-11/current_oiii_optical_ratio/README.md).
All calculation controllers from this sequence are terminal completed. The
previously requested production replacement and Mg I hot-model decisions are
still unanswered; they are not inferred from automatic continuation.

# Latest completed numerical review and pending paired data replacement

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

# Latest Auger source judgment

The final KM article's mean-energy correction was evaluated for 14 subshell
rows mapping to nine historical PH1 channels. Ca II M1 gives -5.022 eV,
which is inadmissible as mean kinetic energy. The earlier requirement for a
joint charge/energy distribution before any heating calculation was corrected:
a consistent mean suffices for the current linear full thermalization of
metal electrons. Charge onsets and energy conventions remain unresolved.
No production heat/rate code changed. `metal_freefree=.true.` remains the
user-requested default exception. See [the source review](validation_2026-10-11/auger_mean_energy/README.md).

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

# Atomic Ionization Stage Expansion: Implementation Status

Date: October 10, 2026 (KST).

## Confirmed scope

The user's revised first-phase goal is **77 metal ion stages and 66 adjacent
transitions**, with explicit levels, collision strengths, and line cooling
covering MOCASSIN v2's supplied atom models. Stages I--VII of the existing
eleven metal elements are included; H and He are counted separately.
The earlier 64-stage target is superseded. Historical 64-stage checks below
remain evidence of completed preparation, not current completion.
The 165-stage complete-charge roster remains outside this first-phase scope.
See the [measured coverage comparison](phase1_77_mocassin_v2/README.md).
The new roster does not establish that terminal electron-impact losses vanish
or that additional inner-shell physics is complete.

## Current completion assessment

October 11 update: the [bidirectional comparison](MOCASSIN_V2_MOCHII_BIDIRECTIONAL_COMPARISON.md)
records the user's revised default scope and incomplete requirements. Metal
electron-impact ionization and line atoms absent from the v2 file mapping now
require explicit `metal_collisional_ionization` and `additional_line_cooling`
input options, respectively. Transition capacity is 66, and the terminal
photon check recognizes bare C VII. A full debug build and scoped native input,
154 line-selection, direct/cached population, energy-loss, and bare/bound
carbon checks passed. Root verified all controller source and result hashes.
These checks do not establish missing-atom completion or a matched nebular run.
The 36-fit count below describes available files in the historical 43-stage
roster, not the number selected under the revised default options. Further
controls for larger shared-ion atoms and other physical differences remain
incomplete. No calculation controller is currently running.

Latest checkpoint: full O IV (147 levels) and O V (166 levels) source-data cooling models are installed and have passed actual installed native and coefficient-reader checks. These imports explicitly retain unverified atomic accuracy. Direct comparison with the supplied v2 files finds 25 missing installed atom/fit pairs and seven smaller explicit atoms; the former 2-of-21 oxygen checkpoint belongs to the superseded 64-stage roster. The production ELEMENT roster still declares 43 metal stages and 32 adjacent transitions, so the approved 77-stage expansion is incomplete. The combined coefficient listing has 39 species including H I and 279 terms. No controller is currently running.

Si I's corrected 35-level/250-radiative-pair atom and seven-term cooling fit
are now installed after all eight validation steps completed. Root checked
63 frozen input hashes and visually inspected the complete profiles. The
maximum fit error is 0.0673568%, table application error is 0.136933%, and
source/native cooling disagreement is 1.08358e-13. The actual installed atom
was opened by the current native driver at 10,000 K and four densities and
exactly matched the validated populations and cooling. Its atomic accuracy
remains explicitly unverified; high-density conditional length/velocity
sensitivity reaches 18.8906%. See the [installation review](validation_2026-10-10/si_i_cooling_products/INSTALLATION_REVIEW.md).

The overall implementation is **not complete**. Direct reading of production
`ELEMENT` records still gives 43 metal stages and 32 adjacent transitions;
the earlier 64-stage, 53-transition roster remains staged, and the revised
77-stage, 66-transition roster is not yet integrated. The active
43-stage roster now has 36 ordinary electron-excitation cooling fits.
C IV, Si V, Ar II, corrected full Mg III, and corrected Ne I atoms/fits are installed and
native-validated. Full O II/O III atoms and the O III optical A correction
are also installed.

The original Ne I 31/142 source-node profiles are preserved as physically
stale after two missing spontaneous decays from level 15 were confirmed.
The corrected private 31-level atom has 144 radiative rows, 465 collision
rows, and 6,975 original nodes. A source-identity test checked both restored
Stout rates, unchanged original rows, and a unique radiatively closed ground
state. Its inferred common BT floor is 1000.36 K; this is not a published
physical validity bound. The corrected source and product controllers completed all native
profiles, energy balances, coefficient-reader checks, and registered-table
checks. Root reviewed all frozen hashes and full curves, then installed
exact validated atom and fit bytes. Maximum source/native cooling error is
2.22e-15, fit error 0.0317686%, and table application error 0.0041481%. A
first pilot failed on an exact comparison of rounded level energies; its
state is preserved and the recovery measured the 0.001 cm^-1 rounding
difference before applying the 0.002 cm^-1 identity tolerance.
Mg I's published 500--20,000 K BSR matrices do not yet supply an accepted
100,000 K model. See the [Ne I/Mg I source audit](source_validity_ne_i_mg_i.md).

The full O II/O III source atoms now have coherent electric-dipole probabilities
with their selected observed/theoretical energies. The generator preserves E1
line strengths with the physical frequency cubed and the actual emitting
statistical weight. Noninverted same-parity rates remain explicitly approximate
source rates. The 196-level O II and 177-level O III atoms are installed;
their numerical records exactly match independently validated candidates.
At 401 temperatures and four densities, native cooling differed from independent
balance by at most `2.9643e-14` and `2.0428e-14`, respectively. The installed
low-density cooling coefficients were preserved because their coronal cooling
was directly verified unchanged. Full 401-by-161 cooling tables passed registration; maximum application errors
were 0.197768% and 0.143563%, respectively;
see
[oxygen radiative consistency](oxygen_valence_convergence.md).
The initial Mg III radiation audit identified 16 inverted electric-dipole
rows. The subsequent NIST identity correction reduces these to 14 and the
full physical-frequency correction is installed and native-validated; see the
[Mg III transition diagnostic](validation_2026-10-10/mgiii_radiative_consistency/README.md).

The dynamic n-level storage and source-node representation are implemented
as an authorized detail of the original expansion. The linked Intel debug
executable was rebuilt successfully. Existing Chebyshev files remain readable.
The new reader stores original collision nodes, preserves source zeros, and
uses the adopted cubic/PCHIP interpolation without clipping negative strengths.
The public consumers allocate their actual level and transition dimensions.

S I, S III, Fe IV, Fe VI, Ar V, and Mg IV passed compiled source/population/
cooling checks over 801 temperatures and four densities. Fe VII and Ne V also passed compiled source/profile checks. The source-node
controller completed all six declared steps, and Ne V's separate corrected-source
controller completed both steps. All five declared registered cooling checks also completed in
the `atomic-registered-remaining` controller. See
[source-node storage and validation](excitation_storage_design.md) and
[controller progress](validation_2026-10-10/source_node_native/controller/progress.md).
The changed hydrogen-impact consumer also passed its published-rate,
Boltzmann, and thermal/emissivity checks for C I, C II, and O I.

Ar II and Ca IV's published ground fine-structure pairs now have checked
Stout readers and source-node atoms. Three source tests and independent
compiled 801-temperature/four-density checks passed. Their collision tables
remain valid only up to their final source nodes; Ar II's table ends at
158,000 K. The source-limited registry now ends at this physical source
ceiling and rejects unsupported higher temperatures. Ground-pair validation does
not establish cooling from higher excited states. Ca IV was installed after
registered-table checks; Ar II is now also installed after GNU and actual
rebuilt Intel source-limited registry checks.

The full 100--1,000,000 K source diagnostics exposed failures of the earlier
truncated Ar V, Mg IV, Fe VII, and Ne V atoms outside their original accepted
interval. Full bound atoms passed all four compiled source checks with both
Intel and GNU builds. All four GNU registered tables passed. The last duplicate
Intel Ne V table was deliberately interrupted after its GNU counterpart passed;
it is recorded as an interruption, not a passed step. These full-bound atoms
remain staged. See [extended validation](excitation_temperature_range_validation.md).

Eleven production low-density cooling corrections were installed after
source/native verification. The earlier corrected five-level C I atom has
since been replaced by its 42-state source atom, alongside complete source
models for C II, N I, N II, and S II. S I cooling and its five-level source-node atom, and S III's full 53-level
source-node atom, were installed after registered-cooling checks. Missing
cooling data, terminal ionization losses, nonadjacent Auger ionization and its
energy accounting, production activation, and matched nebular comparisons
remain incomplete. These statuses do not imply that the expanded production
roster is active.

Ca IV's two-level ground fine-structure cooling and full-source Fe IV's
37-level cooling/atom were also installed after actual source, serialized
native-reader, and registered-table checks. These correct existing missing
cooling within the 43-stage production roster. Ar II is installed following
the validated source-limited registry extension. See
[extended source convergence](excitation_temperature_range_validation.md).

S I's missing low-density cooling has a staged, existing-format constrained
signed-exponential solution. Its earlier single-polynomial native representation was rejected because
source zeros followed by positive values cannot be reproduced by one nonzero
polynomial. Its new source-node representation now passes compiled validation.
S I cooling and its source-node atom are now installed; original failed and
staged products remain preserved. The [storage implementation](excitation_storage_design.md)
records direct source-node evaluation and native checks. The Fe VII final
journal article is now available in the requested references directory;
all 189 printed energies and 57 printed radiative rates were checked as
recorded in [source verification](fevii_final_paper_verification.md).

Automatic completion detection worked, but the first actual automatic model
continuation could not use execution tools (`code-mode host closed its
stdout`). Its successful model exit does not certify a successful scientific
continuation. No status update or source change was made by that execution.
The model watcher was stopped after confirming that execution tools were
unavailable in the continuation session. The main session is continuing
directly; no automatic model continuation is claimed. See
[completion detection record](validation_2026-10-10/completion_watch/README.md).

## Continued implementation and new physical decisions

The user requested autonomous continuation through the approved plan unless a
significant issue was found. The next charge-transfer work is complete:
original N4 and four Cl/Ca total rates are implemented, and a 53-stage staging
roster was validated. Four existing Cl/Ca omissions were corrected in the
43-stage production roster. See
[further charge-transfer validation](kf96_total_recombination_validation.md).

The final 1997 paper identifies the higher-charge Langevin expression as a
capture upper limit. The user supplied the final 2026 Fe and 2017 N IV
journal PDFs, and both have been read directly. Fe ground-level RR/DR fits
are implemented; 1,204 linked rate checks passed. The existing Fe IV
low-temperature recombination omission was corrected, with production
charge/heat consumers checked. See [Fe validation](final_fe_recombination_validation.md).
N IV now uses one checked energy vector in cooling and n-level generation.
Its staged 29-level atom was revalidated against the current source at 801
temperatures and four densities after the interpolation correction. The
maximum full-source cutoff difference is 0.289351%, and native source cooling
agrees within 2.67e-15. Earlier three-density checks remain historical; see
[current source results](validation_2026-10-10/niv_current_source/native/level29/result.json). A 54-stage controlled
rate roster including Fe V passed charge and heat checks. The user subsequently adopted the ten higher-charge Langevin upper-limit
reactions, and a 64-stage, 53-transition rate roster passed linked native checks;
see [Langevin adoption](langevin_adoption_validation.md). The separate
high-charge mean heat uses the physically consistent exit-channel
2.86(q-1) eV; the inconsistent literal KF99 prose is explicitly identified. No expanded production activation or nebular comparison has
occurred. New papers and author data are stored in `~/MoCHII/references/`;
see [reference storage](reference_sources.md).

## Scope and current state

The initial implementation prepares MoCHII for the approved 100 eV target:
64 metal stages and 53 adjacent transitions across the existing 11 elements.
The user approved seven-stage output axes and the old-to-expanded restart
behavior. **Production atomic files still declare 43 metal stages and 32
transitions. The 21 proposed stages are not yet activated.**

The [source inventory](source_inventory.md) and its
[JSON records](source_inventory.json) identify actual local sources, hashes,
rate coverage, thresholds, and unresolved physical inputs. Direct execution
found PH1 and Voronov CI records for all 21 additions, PH2 for 19, local
Badnell M=1 RR and DR for 18 each, and CHIANTI line/collision files for 20.
Cl VII has Stout data whose collision-temperature grid begins at 9,800 K.
These are source measurements, not validated physical models.

## Completed implementation

| Change | Behavior and reason |
|---|---|
| Shared stage storage | `define::MAX_ION_STAGES=7` and `MAX_METAL_ELEMENTS=11` supply the species and line arrays. The packet capacity remains 64, sufficient for the 53-transition target. |
| Physical ionization energies | `ETH`, PHOTO cutoffs, and occupied-valence PHSHELL cutoffs use the identified NIST ASD snapshot. Inner PH1 edges and Voronov CI parameters retain their independent source values. All eleven current element files were updated consistently. |
| Terminal photon-band check | The main transport setup refuses a band above an active terminal ion's physical ionization potential. The current 43-stage roster therefore rejects 100 eV; the proposed 64-stage terminal edges all exceed 100 eV. This check does not establish collisional or Auger completeness. |
| Charge-transfer fit ranges | Thirteen existing KF96 higher-boundary fits use individually resolved form-8 tables; N's first-boundary fits use 100--50000 K. Private recombination-table capacity is 1024, with finite positive increasing-temperature checks; heat tables retain capacity 40. |
| Nine added charge-transfer boundaries | Selective generation combines KF96 published total rates with KF99 Table 1 mean kinetic heat. A 52-stage, 41-transition staging roster passed compiled rate and heat checks. Ground-state reactants and low-density radiative relaxation are assumed; production remains at 43 stages. |
| Element reader | Rejects invalid stage/transition bounds, nonconsecutive transitions, missing required rate records, duplicates, unknown records, nonfinite principal coefficients, invalid thresholds and RR parameters, and a mismatched element identity. DR counts are checked before indexing. |
| Adjacent ionization equilibrium | Direct and cached calculations call `adjacent_ion_fractions`. Logarithmic populations prevent ratio/product overflow. Exact zero-rate links are treated without a rate floor; multiple independent stationary populations cause an explicit nonunique-equilibrium error. This routine does not implement Auger transitions. |
| Cooling registration | Internal capacity is 78 (H I plus all supported metal stages); storage of the actual cooling tables remains allocated only when needed. |
| Fe level selection | `_full` files are selected only for Fe II/III in cooling and line calculations. Later Fe stages use ordinary filenames. |
| Missing line cooling | Setup prints every stage lacking a cooling coefficient file. The message exposes missing physics; it does not supply the omitted cooling. |
| Free-free cooling | Gaunt factors cover charges 1-6. The preserved Hummer implementation is used only within its stated range; outside it, the existing van Hoof table is integrated. The Python table is regenerated from direct Fortran execution. |
| HYBRID CI coefficients | All 53 target metal transitions, plus three H/He transitions, match the corresponding Cloudy c25.00 constant-ratio entries. Raw Voronov behavior and existing lower-stage factors are preserved. This remains a constant-ratio model, not direct Dere data interpolation. |
| Selected generation | Element, cooling, and n-level generators have main entry points, explicit selection, and output-directory options. Files are replaced atomically; element files are rendered before any selected output is replaced. Added boundaries with unresolved charge exchange/product energies fail before writing. Fe ground-level q=3--6 uses reviewed Zhang 2026 RR/S1 DR; the ten higher-charge Langevin limits are selected and staged; added cooling still blocks activation. |
| Expanded restart | Common old rates are retained; added rates start at zero. Missing transition rates refuse `gas_niter=0` and require convergence before output. A fresh transport pass precedes each solve in the existing main loop. |
| Optional reaction outputs | Each requested charge-exchange output is initialized and filled independently, including when the other optional output is absent. |
| Build consistency | A compiler-settings prerequisite recompiles objects when the compiler or flags change, preventing native debug links from retaining incompatible IPO objects. |

## Existing defects corrected

The Gaunt correction also affects existing charges at low temperature: the
former Hummer extrapolation produced negative free-free cooling coefficients
outside its stated range. Previous low-temperature results using those values
are stale until recalculated. Their small or large numerical impact does not
justify keeping an invalid coefficient. Existing numerical products are
preserved; no nebular model has been recalculated in this work.

An old restart could contain more elements than the new run. Previously,
limiting read blocks to the new element count could omit a later active
element's rates after earlier disabled elements filled that count. The reader
now retains the supported old roster before matching active element names.
The synthetic restart test places a disabled C block before the active Fe
block to check this correction.

The CHIANTI `.elvlc` parser now follows the official reader's first 87
fixed-width characters. The previous last-three-fields implementation still
misread N IV records with whitespace-separated extra metadata. The source
record for level 3 now reads g=3, E_obs=67,272 and E_th=67,458 cm^-1.
Optional comma-separated and whitespace-separated columns are both tested.

The previous cooling generator assumed that all ions occupy the ground
state at zero electron density. C III and Si III each have two radiatively
stable states in the supplied model. Their collision-determined stationary
populations and radiative cascades are now included, with independent
finite-density balance checks. Both production cooling fits were corrected;
old coefficients are preserved in the validation directory. The aggregate coefficient table reflects the installed primary products.
After the S I, Ca IV, and Fe IV additions it contains 31 ions and 182 terms; all primary coefficients
are retained at their serialized precision. Dependent older model results are stale until recalculated.
No nebular model has been recalculated.

A signed exponential expansion fits total cooling when a positive-coefficient
sum cannot reach the required accuracy. Individual fit coefficients are not
emission rates. The serialized total is positive and within 1% on 801 points
between 1,000 and 100,000 K. Additional positivity checks between 1 and 10^9 K
do not establish extrapolation accuracy. The unchanged Fortran file format
and evaluator accept these coefficients. See [cooling validation](cooling_validation.md).

The original mixed-energy N IV generation was blocked by 14 radiative
ordering reversals. This earlier defect is resolved by the checked Stout energy
vector used consistently in cooling and the accepted staged 29-level atom.
Original external files and failed products remain preserved.

## Validation

Validation commands and checks run under the persistent
`~/.codex/scripts/calculation_workflow.py` controller in `tmux`, with finite
timeouts, source hashes, command/PID records, expected artifacts, and JSON
acceptance checks. Failed build and fixture attempts are preserved, with
explicit recovery notes; they are not overwritten or silently restarted.

The completed native debug verification is recorded in
[full check results](validation_2026-10-09/final_preparation_corrected/products/full_result.json)
and the follow-up
[optional-output check results](validation_2026-10-09/optional_reaction_output/products/extended_result.json).
The optimized/IPO build and balance checks are recorded in
`validation_2026-10-09/optimized_matching_flags/`; the completed rate and
negative-input checks are in
[optimized extended results](validation_2026-10-09/optimized_precision_check/products/extended_result.json).
Legacy Gaunt quadrature uses default-REAL intermediates: the observed largest
debug/optimized reference difference was 3.612e-8. Only the Gaunt comparison
uses a 1e-7 relative tolerance; other rate comparisons retain 1e-10.

- The complete executable builds with bounds and floating-point checks.
- 25 controlled adjacent-stage checks cover analytic balance, normalization,
  extreme rates, one-way limits, disconnected systems, and invalid rates.
- A synthetic seven-stage Fe fixture checks storage, finite normalized
  fractions, direct/cached charge closure, all six packet transitions, and
  stage/transition FITS output readback. It is not physical Fe atomic data.
- Old three-transition rates are preserved in a six-transition restart;
  the three new rates are zero, reconvergence is required, and a fixed-state
  incomplete restart is rejected.
- Existing charge-exchange conservation, cell convergence, detailed balance,
  H/He charge exchange, and heating checks pass with the actual production
  element files. After the October 10 rate corrections, the measured largest
  reaction conservation deviations are 3.565e-16 (ionizing) and 3.169e-16
  (recombining).
- All 246 Gaunt nodes are positive. Charges five and six at 100 and 1,000 K
  agree with a separately coded 4,096-point integration of the existing table
  to a maximum relative difference of 5.3283e-7. This checks quadrature and
  table usage, not the accuracy of the source atomic calculation.
- Python rate references agree with the compiled Fortran, including all six
  Gaunt charges. All 56 selected HYBRID factors match the actual Cloudy table.
- Twelve malformed element files are rejected with the intended diagnostics.
- The initial nine generation tests cover file safety, both optional energy-column
  formats, a synthetic atom with two stable states, actual C III/Si III
  against independent finite-density statistical equilibrium, the C II
  single-state limit, and rejection of inconsistent N IV energies.
- Corrected C III and Si III fits have maximum dense-grid errors of 0.1279%
  and 0.4388%. The compiled optimized Fortran evaluator agrees with Python
  to 1.09e-13 and 3.05e-13 respectively across 801 temperatures.
- The S V staged fit passes with 0.278725% maximum error and 1.84408e-13
  Fortran/Python difference. Its density and level convergence remain open.
- All 20 added CHIANTI energy files are read in the source audit; reversed
  radiation energies and nonpositive collision gaps are listed for review.
- Generation headers use the actual KST date; no unrelated products were
  regenerated solely to correct old header dates.

No full regression suite, matched Cloudy/MOCASSIN model, nebular profile
comparison, or production stage-expanded run has been performed. This scope
tests the affected internal storage, rates, parsing, cooling, and I/O; it does
not establish physical completeness of the proposed atomic roster.

The Hummer validity limits and van Hoof grid limits were checked against the
preserved implementation/data headers. Those earlier validation steps did not newly verify the Hummer and van Hoof
journal articles: ADS name resolution and publisher access failed. This is
distinct from the subsequently verified Fe and N IV final PDFs. No new
Hummer or van Hoof published accuracy claim is made.

## Reproducible entry points

The October 10 implementation and measured checks are detailed in
[physical threshold validation](physical_threshold_validation.md). The complete
debug executable, 77 compiled potential entries, all 32 production photon
thresholds and energy budgets, 26,049 compiled charge-transfer samples, and
25 focused Python tests pass. Three malformed rate tables are rejected.
Completed controller records are in `validation_2026-10-10/`; failed attempts
and their recovery notes are preserved. Current calculations are identified
in [remaining cooling workflow](remaining_cooling_workflow.md).

The separate [charge-transfer audit](charge_exchange_expansion.md) identifies
the original nine added-boundary KF96 total rates with accepted KF99 mean kinetic heat
for selective staging. The [mean heat validation](kf99_mean_heat_validation.md)
records the final publisher source, approximation limits, and the current
28 passing Python checks. Compiled evaluation covered 18,027 rate samples
(maximum relative interpolation error `5.1541e-6`) and 27 heat samples
with zero discrepancy from the adopted published means. Charge and reaction
accounting passed at three prescribed temperatures. This does not establish
state-resolved branching, reverse rates, or a nebular model comparison.
[Fe recombination](fe_recombination_expansion.md) includes a 301-temperature
source-fit comparison and a visually inspected vector PDF; that historical comparison preceded selection of the final Zhang RR/DR fits
and is superseded by the measured Fe source and implementation reports. The
[energy consistency audit](atomic_energy_consistency.md) shows that switching
all levels to theoretical energies would not repair the CHIANTI source
inconsistencies. None of these investigations activates the missing stages.

These commands are implemented. Use an explicit staging directory; do not
replace production files until their sources and numerical fits are accepted.

```bash
python3 tools/fitting/inventory_stage_expansion.py
python3 tools/fitting/make_element_data.py --element c --output-dir /tmp/mochii_atomic_selected
python3 tools/fitting/fit_cooling.py --ion CI --output-dir /tmp/mochii_atomic_selected --no-plot
python3 tools/fitting/fit_nlevel.py --ion CI --output-dir /tmp/mochii_atomic_selected
python3 tests/ion_stages/run_checks.py build --work-dir /tmp/mochii_stage_checks
python3 tests/ion_stages/run_checks.py full --work-dir /tmp/mochii_stage_checks
```

`--stages ELEMENT=N`, `fit_cooling.py --ion-spec ELEMENT:STAGE`, and
`fit_nlevel.py --ion-spec ELEMENT:STAGE:LEVELS` support explicit selection.
The element generator intentionally refuses unresolved added boundaries.
The selected n-level count is a candidate model size, not a convergence result.

## Next physical work

1. The ten higher-charge reactions now use the user-adopted Langevin
   upper-limit approximation and have passed the 64-stage rate checks. N4 uses the original published
   fit, and the four existing Cl/Ca omissions have been corrected. Fe4 total
   rate and mean heat are available with its radiative-transfer caveat. Refine product branching, reverse reactions, and cascade-photon
   treatment where needed beyond the accepted ground-state, low-density
   mean heat approximation. The generator prevents silent acceptance of
   unresolved added boundaries.
2. Fe final-publisher verification, linked ground-level RR/DR and captured-electron cooling checks are complete.
   Preserve the measured compact-fit limitations; complete added-stage cooling
   before Fe VI/VII activation.
3. N IV energy consistency and staged 29-state cooling/density checks are complete.
   C V, N VI, S VII and Ne IV staged fits now pass dense source and native
   checks; see [cooling validation](high_excitation_cooling_validation.md).
   Continue convergence checks for the other selected atoms. The user adopted
   constant Cl VII Upsilon below 9800 K; its Stout conversion and numerical
   verification are staged. Address current-stage cooling gaps as well.
4. Complete the collisional terminal-stage assessment and nonadjacent Auger ionization/energy
   accounting, then activate reviewed production data and perform matched
   fixed-state and nebular comparisons with full profile inspection.

The `mochii-remaining-cooling` controller completed all 12 selected fits
and their native checks at 09:09:57 KST on October 10. No calculation remains
running in that sequence. All 12 figures were subsequently inspected; remaining physical limitations are
recorded with each source. See [workflow and limits](remaining_cooling_workflow.md).
No nebular model calculation or automatic model-resume watcher is running.

## Subsequent excitation-source completion

O V and Cl VI energy consistency, staged cooling and 29-level atoms now pass
complete 801-temperature, three-density comparisons with full source atoms.
See [ICFT validation](icft_energy_atom_validation.md). All 21 added ions now have staged low-density cooling-fit/native checks,
including Cl VII with the adopted cold approximation; this count does
not establish finite-density convergence for every other added atom.
Cl VII uses the adopted cold constant-Upsilon approximation described in the
[source-domain assessment](clvii_temperature_and_core_states.md). Its low-density
fit now passes 801-temperature checks (maximum source-fit error
0.5728601711129544%, native fit-evaluation difference 4.440892098500626e-16).
The 32-bound-level finite-density validation passed all 2,403 temperature/density
samples (maximum source/native cooling difference 0.1348220130716005%).
All 496 serialized strength curves pass 1,601-temperature checks at
0.1999403980527781% maximum relative error; cold strengths are constant.
Cooling and all excited-state population profiles were visually inspected.
The final Cl VII controller completed; none remains running.
All three ICFT controllers completed; none remains running. Production remains
at 43 metal stages. The next independent terminal-CI source assessment is
preliminary and is not a runtime implementation or completed physical audit.

## Latest scoped numerical verification

The main session inspected the complete current C V, N IV, and N VI cooling
fits and the N IV full-source temperature/density profiles. Current rounded
low-density fit errors are 0.855934%, 0.469710%, and 0.796310%, respectively,
on the checked 1,000--100,000 K grid. These products remain staged; their
fit accuracy does not activate the expanded roster. N IV's 29-state source
atom also passed actual registration of 64,561 cells and public cooling
application checks. The greatest checked interpolation/application difference
is 0.163454%.

A separately confirmed existing numerical defect causes positive subnormal
registered cooling to be diminished by dividing twice by a `tiny()` floor.
The installed correction stores absolute cooling and preserves exact zeros
separately. Its focused subnormal/zero tests and S I, Fe IV, and Ca IV
registration checks passed. The actual production build and four tests linked
to its rebuilt objects also passed; these tests explicitly enable gradual
underflow and do not change production compiler defaults. The separate
existing partial-atom residual approximation failed seven of nine corrected
combined-cooling checks (maximum C II discrepancy 27.7687%) and remains
under investigation.
See [subnormal cooling correction](subnormal_excitation_cooling_correction.md).

The isolated nonadjacent stationary-ion entry point passed main-session
Intel debug/no-FTZ and GNU bounded/FPE-checked builds. Analytic Auger
fractions, reduction to adjacent equilibrium, unique closed classes,
invalid inputs, unresolved ill-conditioned networks, and representable
1e-310 and 1e-320 fractions were checked. The module is not included in
the production object list and no production caller invokes it. This is
a tested numerical component, not implemented Auger ionization or heating.
See [main validation](validation_2026-10-10/nonadjacent_stationary_main_review/main/result.json).

The source-node writer now accepts exactly degenerate electron-collision
pairs with deterministic indices while retaining their original positive BT
scaling energy. No artificial physical level splitting was introduced. Six
writer checks passed. An actual compiled two-level atom at five temperatures
and four densities recovered populations 3/8 and 5/8 exactly and zero cooling,
as required by equal physical energies and statistical weights. Positive
radiative transitions still require a nonzero physical gap. See
[native degeneracy check](validation_2026-10-10/degenerate_collision_native/result.json).

C V, N VI, and the source-converged 45-level Ar VI atom passed current-module
registered-table and public-application checks. The maximum measured
application differences were 0.113258%, 0.001781%, and 0.310620%, respectively.
Their isolated finite controller completed all four steps; products remain
staged. See [current registration](current_source_registration.md).

Full source-node C I (42 levels), C II (53), and N I (50) atoms now replace
the former five-level models after native source/population, public cooling
registration, and affected H-impact checks passed. Their density-ratio cooling
removes the physically unjustified unsuppressed fit residual. The earlier
nine-ion failed gate remains historical. Its seven failing compact atoms have
now been replaced after focused source/native and registered-table checks;
the complete historical nine-ion command was not rerun. Low-density fit files and the 43-stage production
roster were preserved. See [compact-atom correction](existing_compact_atom_correction.md).

Full N II (58 levels) and S II (70) source-node atoms are also installed
after full-source/native and registered-cooling checks. Their maximum
checked public-application differences are 0.138828% and 0.136925%.
O II and O III were subsequently replaced by the independently validated
196- and 177-level atoms. Their radiative probabilities use physical
photon energies with preserved E1 strengths, not symmetric pair ordering.

### Production validation of the source temperature floor

The actual production target rebuilt successfully, and a driver linked to its
`nlevel_mod.o` passed the six-type source-strength, population, detailed-balance,
and cooling checks. Eleven invalid `BT_TMIN_K` records were rejected. See
[measured results](validation_2026-10-10/bt_floor_production_validation/result.json).
The source approximation holds Upsilon constant below its stated physical
temperature; Maxwellian factors retain the actual temperature. The production
metal roster remains 43 stages and 32 adjacent transitions.

### C IV and Si V cooling installed

Source-node atoms with 24 and 169 bound levels, respectively, and annotated
low-density fits are installed in `data/atomic`. The maximum rounded fit
errors over 801 temperatures from 1,000 to 100,000 K are 0.344146% and
0.929616%. The adopted constant-Upsilon continuation starts below 3,200 K
and 5,000 K; its physical accuracy is not certified by numerical fitting.
Source/native full profiles and registered cooling tables passed before
installation; installed bytes were verified and opened by native readers.
No numerical products were regenerated for installation alone.

The active 43-stage roster now has 32 ordinary electron-excitation fits;
[the current coverage record](validation_2026-10-10/production_cooling_coverage_current.json)
distinguishes absent files from all other physical terms.

### Existing O III optical radiation corrected

The two optical A values now use final Storey and Zeippen (2000) Table 3
totals: 6.791e-3 and 2.046e-2 s^-1. The actual native 5007/4959 energy-flux
ratio changes from 2.89105976 to 2.98396833. Native emission and cooling were
checked at 801 temperatures and four densities, agreeing with independent
Python cooling to 1.49e-14 relative. The largest finite-density cooling
change was 0.0290678%; the coronal cooling is unchanged to numerical
precision, so its low-density fit was preserved. This is a correction of an
existing source-attribution/physical-expression defect, separate from stage
expansion. See [native measurements](validation_2026-10-10/oxygen_radiative_provenance/native_main/result.json).

### Source-limited registered temperature grids and Ar II

The cooling registry ends at the smallest last source temperature among
`TEMPERATURE_TABLE` collision rows, preserving its original nodes and adding
the exact final endpoint. BT/legacy grids retain their original 401 rows.
The final interval, including a ceiling inside the original final grid cell,
and explicit nonfinite/zero-input rejection passed GNU and actual rebuilt
Intel validation. Four subnormal/exact-zero fixtures passed again because
they exercise the modified common interpolation path. Public interfaces,
source nodes, and user inputs were unchanged.

Ar II's published two-state ground fine-structure cooling is now installed.
Its source ends at 158,000 K; unsupported higher-temperature use is rejected
rather than silently extrapolated. The existing 0.202152% fit over
1,000--100,000 K was preserved. Production now has 33 ordinary electron
cooling fits among 43 metal stages, while the 64-stage roster remains staged.
See [source-limited validation](validation_2026-10-10/source_temperature_registry/README.md).

### Mg III full-bound private validation

The corrected 167-level atom preserves higher-level E1 strength with the
physical frequency cubed and emitting statistical weight. All 801 temperatures
and four densities passed native/Python comparison, with maximum cooling
error 1.13243e-14 and population error 2.10942e-15. The 27-level truncation
loses up to 17.2406% of full-atom cooling and is rejected. The original
controller records a failed final artifact-name check despite a successful
numerical comparison; a verification-only audit checked the preserved results
without recomputation. These pre-correction profiles are historical conditional diagnostics; the
subsequent NIST identity correction, full validation, and installed products
are recorded below. Unresolved multipole identification and source accuracy
remain explicit. See
[the physical assessment](mgiii_radiative_consistency.md).

### Mg III observed-energy assignment defect

Direct retrieval of the official NIST Mg III energy table found that source
level 167 (`2s2 2p5 7p 1S`, J=0) carries an observed energy of
628148.375 cm^-1 belonging to a `7d` J=1 level (NIST 628149 cm^-1).
The Stout NIST-derived list independently contains the latter J=1 identity.
The mismatched observation must not set the J=0 photon or collision energy.
The correction uses its original theoretical energy, 634430 cm^-1, while
preserving every external input byte and the other level identities. The
previous full167 profiles remain a conditional diagnostic and are physically
stale for production acceptance. The guarded selection correction and new scoped full167 validation
subsequently passed; the installed correction is recorded below. This corrects an established source defect
within the approved excitation scope and does not change a public interface
or the stage roster.

Two preexisting plot-label typos (`C III]` and `Si III]`) encountered in
`fit_cooling.py` were corrected to `C III` and `Si III`. Only labels changed;
no coefficients or existing figures were regenerated for that correction.

### Corrected full Mg III cooling installed

The NIST-guarded 167-state source atom and low-density fit are installed.
Source ID167 uses its own theoretical 634430 cm^-1 energy; all other selected
energies and original collision nodes are unchanged. The source preserves
higher E1 line strengths, retains first27 A values, and explicitly approximates
112 unidentified multipole rates. Native/Python cooling agreed to 2.57572e-14
relative over 801 temperatures and four densities. Gross collision/radiation
energy residual was at most 4.13663e-14. The rejected first27 truncation
differs by up to 17.1870% from the corrected full source. The low-density
fit error is 0.210705%; its actual native coefficient reader differs by
1.13554e-11 on all 801 points. The full 64,561-cell table passed; maximum
application error was 0.000534832%. Main review verified both controllers'
frozen hashes and inspected the complete cooling and fit PDFs before copying
unchanged numerical records. The actual production-module driver opened the
installed atom at 10,000 K and four densities, reproducing the validated
profile exactly. No full profile or table was rerun for installation.

At the Mg III installation checkpoint, production had 34 ordinary electron-excitation fits among 43 stages/32
adjacent transitions. The collected text table now contains 35 primary files
including H I; collection performed no fitting or atomic-data regeneration.
Ne I validation followed and is complete in the next section. The overall expansion remains incomplete
overall; the 64-stage roster is not activated. See
[installation and main review](validation_2026-10-10/mgiii_corrected_level_identity/main_review.json).

### Corrected Ne I cooling installed

The complete 31-state Ne I atom restores only the two identified omitted
spontaneous channels, preserving all 142 original A values and all 6,975
collision nodes. Its production atom and fit bytes exactly match validated
private products. Scoped full profiles and source strengths, statistical
and energy balances, the GNU coefficient reader, and the Intel density-table
registration passed. The production debug reader opened the installed file
and reproduced the validated 10,000 K/four-density records exactly.

At the Ne I installation checkpoint, production had 35 ordinary fits and
the aggregate listed 36 species including H I and 233 terms. The subsequent
Si I installation raises the current measured counts to 36 ordinary metal
fits and 37 aggregate species/240 terms. Seven fits remain absent: Ar I,
Ca I, Ca III, Cl I, Cl V, Fe I, and Mg I. Production still has 43 stages
and 32 adjacent transitions. Expanded production activation and the overall
project remain incomplete.

## Resumed 64-stage implementation: O IV source review

The user has resumed the 64-metal-stage objective. O IV full-atom numerical
validation is the next prepared calculation, but root first detected a
source-frequency inconsistency requiring radiative provenance review. The
current selector returns 147 states, 2,319 radiative rows, and 10,731 collision
pairs. For 343 higher E1 pairs, stored gf and A do not imply a frequency
within 1% of any available observed/theoretical/wavelength frequency. No
correction is inferred from this diagnostic alone. See the
[direct source review](validation_2026-10-10/o_iv_main_radiative_review/README.md).
Production remains at 43 metal stages and 36 ordinary cooling fits.

## Completed O IV conditional full-atom verification

The full O IV numerical sequence completed all five steps. Root verified 53 frozen input hashes, five step artifacts, saved native products, and seven full-profile PDF pages. The selected atom has 147 levels, 2,319 radiative pairs, 10,731 collision pairs, and 160,965 original nodes. Source/native cooling agrees to 1.39888101103e-14 across 401 temperatures and four densities; normalized energy balance error is 4.37526002907e-14. The 20-term fit reaches a maximum 0.287515276% error on 801 temperatures, and actual registered application reaches 0.277913875%. See the [root completion review](validation_2026-10-10/o_iv_full_completion/root_completion_review.json).

No O IV production installation follows from these numerical passes. Root independently reproduced the original ADF04 comparison: 303 selected CHIANTI A values differ from the original by more than 0.1%, while 256 of the 343 higher gf/A discrepancies retain the original A within 0.1%. This establishes that arbitrary wavelength-based A rescaling is unsupported, not that either source is automatically physically wrong. Radiative source correspondence remains under review. The next O V full-bound sequence is being prepared separately; its temperature bounds must be established for that ion from the actual source.

## O V distributed temperature boundary

Root inspected the final Fernandez-Menchero et al. (2014), A&A 566 A104 Methods and the actual Stout `o_5.coll` temperature row. The final paper prints a general 16,000--160,000,000 K interval, while the identified distributed O V strengths have 13 finite temperatures from 5,000 to 50,000,000 K. These are not treated as identical claims. The current generator explicitly records the 5,000 K table boundary and applies the already authorized constant-Upsilon approximation below it; actual gas temperature still enters Maxwellian factors and detailed balance. Root directly verified constant finite strengths at 1,000, 2,000, and 5,000 K for all 13,695 selected collision pairs, and checked all recorded source paths. Original published-source accuracy remains unverified. The O V writer includes the actual Stout collision file in its source hashes and distinguishes a distributed table boundary from a published physical validity limit.

The executed O IV generator files were preserved before this O V change in [the source archive](validation_2026-10-10/o_iv_full_completion/executed_source_archive/record.json). Existing production products were preserved; the change applies to newly generated O V diagnostic atoms.

## O IV installed; O V numerical sequence completed

Subsequent root review confirmed that the invoked O IV n-level cooling path uses level energies and imported A, with no dependence on the inconsistent gf or wavelength metadata. Root independently reproduced the updated original-index source audit and inspected the actual writer/native consumers. No wrong A was established. The source probabilities are explicitly retained as atomic-model estimates, with unverified detailed replacement provenance and absolute accuracy. This supports qualified adoption of the full-bound cooling model without an unsupported wavelength rescaling. Earlier conditional validation records remain preserved. See the [installation review](validation_2026-10-10/o_iv_full_completion/installation/README.md).

The O IV 147-level atom and 20-term fit are installed, with every validated numerical record preserved. The current generator reproduces the atom records. The actual installed native consumer at 10,000 K/four densities exactly matches saved cooling and populations, and the installed GNU coefficient reader passed 801 temperatures. The aggregate listing now contains 38 species including H I, with 260 terms. Production ELEMENT records remain unchanged at 43 metals/32 transitions; only 36 ordinary metal fits belong to this active roster. O IV data installation does not activate a fourth oxygen stage.

O V's five-step controller also completed. Root verified all 58 frozen hashes, all step artifacts and saved native products, and seven full-profile plot pages. Its 166-level source/native cooling agrees to 2.02060590482e-14, normalized energy residual is 6.63342718617e-14, maximum 19-term fit error is 0.139092942%, and registered application error is 0.099054140%. Its model remains conditional and uninstalled pending radiative-source review. See the [root O V review](validation_2026-10-10/o_v_full_completion/root_completion_review.json).

The readiness audit's previous behavior silently replaced its original snapshot. Root corrected this internal report writer to refuse overwriting it, verified preservation, and saved the post-O-IV inventory separately. This report fix changes no physical data or numerical calculation. Two delegated source investigations ended with a workspace-credit exhaustion error after saving their earlier findings; their unfinished temperature/source reviews are not represented as completed. No calculation controller is currently running.

## Latest oxygen completion checkpoint

O V is now also adopted as a qualified source-data cooling model and installed. As with O IV, the actual collision/radiation calculation uses identified level energies and imported positive A values; inconsistent unused gf/wavelength fields alone do not establish incorrect A. Source probabilities remain estimates, absolute accuracy is unverified, and the observed/theoretical selection does not recalculate scattering resonances. The 5,000 K distributed collision boundary and its adopted cold continuation are explicitly distinguished from the general final-paper interval. This qualification is not a new photoexcitation or absorption model. The [installation review](validation_2026-10-10/o_v_full_completion/installation/README.md) records the limits.

Every validated numerical record is preserved. The current generator reproduces the O V atom, installed native consumers at 10,000 K/four densities exactly match saved cooling/populations, and the actual installed GNU coefficient reader passed 801 temperatures with maximum relative discrepancy 6.12732087291e-12. The aggregate listing now has 39 species including H I and 279 terms. The [latest direct inventory](validation_2026-10-10/upper_stage_readiness/post_oxygen_installation.json) finds 2 of 21 added-ion atoms/fits installed, with 19 remaining. This supersedes the earlier O V pending-installation checkpoint. Active production still has 43 metal stages/32 transitions and 36 ordinary metal fits. Final 64-stage activation, seven original-roster source gaps, Auger charge/energy closure, and terminal-stage validity are unfinished. No background controller remains running.

## October 11: neutral-H option and revised rate preparation

The neutral-H collider requires `par%neutral_hydrogen_excitation=.true.`;
its default is `.false.`. Thermal and emitted-line consumers use the same
choice. A rebuilt native validation passed 36 installed Ne I/O I
thermal/emitted-energy cases, as well as the earlier 154 scope decisions and
ionization checks. Root checked all controller source and artifact hashes.
The library's explicit collider routines remain available for diagnostics.

The generator includes the published Mg4+ charge-transfer fit and 8.60 eV
mean heat, checked against the final KF96/KF99 tables. The previously adopted
q>=5 Langevin upper-limit rule now covers all required high-charge boundaries
for the revised roster. Twelve scoped Python tests passed. The 77-stage,
66-transition rate files are staged separately from production and retain the
unchanged S II product-state dependency. Missing line atoms, shared-ion level
selection, and other extra-physics controls remain incomplete.

The actual staged 77-stage/66-transition reader and conservation check also
passed: 198 reaction/temperature samples and the 100 eV terminal photon-band
entry point, including bare C VII. Root verified source/product hashes. The
first attempt lacked the S II dependency and its failed record is preserved.
No production rate file was replaced; line-atom completeness and remaining
extra-reaction default controls are not established by this check.

## October 11: high-charge option completed

`high_charge_exchange` defaults to `.false.`, as requested for additional
physics. Its central rate selection covers direct and cached populations,
hydrogen charge transfer, and mean kinetic energy. Both actual namelist
choices passed 198 staged native samples, including 66 high-charge samples.
Root verified 98 recorded hash entries and both controller artifact hashes.
The controller completed and has no next calculation scheduled. See
[the scoped validation](validation_2026-10-11/high_charge_option/README.md).
The four implemented additional-physics switches default to `.false.`.
Further controls and 77-stage line-atom completion remain unfinished;
production is still 43 metal stages/32 transitions.

## October 11: thermal options and explicit exception

`charge_exchange_heating` defaults to `.false.`; `metal_freefree` defaults
to `.true.` as explicitly requested by the user. Rebuilt debug objects and
108 native thermal cases passed. The four input choices change only the
selected kinetic energy/free-free loss at fixed temperature, with metal-only,
helium-only, and combined reaction consumers checked. Root verified frozen
hashes and both artifacts. No spatial run or unrelated data regeneration was
performed. See [the thermal validation](validation_2026-10-11/thermal_defaults_complete_inputs/README.md).

The MOCASSIN source reader can now return complete collision tables and
spontaneous probabilities for import diagnostics without changing its default
inventory. All 56 available files parsed successfully: only Ca VII has
nonpositive radiative energy gaps (four), and only Mg I has a single
collision temperature. No pair mixes zero and positive collision strengths.
These measured conditions are recorded in
[the source-table import audit](phase1_77_mocassin_v2/source_table_import_audit.json).
No candidate line atom was installed by this source-table audit.

## October 11: MOCASSIN finite-level candidates

54 source atoms are converted and passed actual native level/cooling/line
consumers against an independent dense matrix at 5,408 temperature/density
cases. Root verified frozen hashes and numerical outputs. Fe II retains
142 levels with 1,144 positive A values, correcting the external reader's
orientation error; Fe VI retains 19 levels. Mg I's singleton temperature
and Ca VII's four nonpositive photon-energy pairs remain under physical
review. The positive log-log PCHIP correction differs from v2's natural
cubic: the original expression gives 178 negative sampled strengths for
N I and one for Ne V. This diagnostic does not constitute a complete v2 run.

No candidate was installed, no production rate file was replaced, and
production remains at 43 metal stages/32 adjacent transitions. Figures and
visual profile review, fitted cooling, shared-ion selection, and full
77-stage activation remain incomplete. The source-atom controller completed
and has no successor scheduled. See
[the source-atom validation](validation_2026-10-11/mocassin_source_atoms/README.md).

## October 11: source profiles and compiled cooling fits

All 54 vector profile pages were visually inspected at four electron
densities, including every supplied level. Root checked all figure source
and product hashes. Independent source arrays are preserved. The strict
low-density limit differs from the native 1e-6 cm^-3 column by at most
2.44430287891e-7 relative on the sampled profiles. See
[the full profile review](validation_2026-10-11/mocassin_source_atoms/root_profile_review.json).

52 candidate exponential cooling fits now pass the actual production
reader/evaluator on 41,652 temperatures, with a maximum relative source
difference of 0.0049953117071. The 0.5-percent
criterion is unchanged. Failed preparation series and their products are
preserved; only unfinished fits were recalculated. Two exponential fits
(Ne VI and Fe V), source reviews for Mg I and Ca VII, numerical temperature
range enforcement, shared-ion model selection, and production integration
remain incomplete. No atom, fit, or expanded rate chain was installed by
these checks; production remains 43 metal stages/32 transitions.

The diagnostic final-temperature roundoff problem is corrected. A misleading
shared Python table error now identifies the actual table ceiling rather
than labeling every caller Cl VII. The low-density physical equations were
exposed through a source-independent entry without changing their body;
existing C III, N IV, O IV, and Fe II calls were checked. All calculation
controllers for this checkpoint have reached terminal states; no successor
calculation is currently scheduled. See
[the compiled coefficient validation](validation_2026-10-11/mocassin_cooling_native_exact_endpoints/README.md).


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
