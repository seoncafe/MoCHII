# Interim atomic-physics scope relative to MOCASSIN v2

> **Dated record, corrected 2026-10-11.** The current state is
> [ATOMIC_PHYSICS_RESTORATION_RESULTS_2026-10-11.md](../ATOMIC_PHYSICS_RESTORATION_RESULTS_2026-10-11.md):
> stages I-VII of all 11 metals (77 / 66) in production; metal absorption,
> rates and heat on the `phfit2` subshell sum; default
> `metal_excitation_model = 'chianti'`; metal collisional ionization,
> neutral-H excitation, charge-exchange heating and the He I 2^3S coupling on
> by default. Corrections follow
> [the audit](../ATOMIC_STAGE_EXPANSION_AUDIT_2026-10-11.md), Section 5.

**October 11 status correction:** The historical comparison below describing
shell-summed metal photoionization in MoCHII predates the approved
`metal_inner_shell_absorption=.false.` default. The invoked default now uses
the outer shell; explicitly enabling an active inner shell stops with a
diagnostic because vacancy products and energy accounting are not qualified.
The 77-stage roster remains staged, with 43 metal stages and 32 transitions
active. See the [current result](../ATOMIC_IONIZATION_STAGE_EXPANSION_RESULTS_2026-10-11.md).
(Corrected 2026-10-11, audit Section 5 item 2: the outer-shell default was
not an improvement on the shell-summed treatment described below. It
deleted valid absorption above `einn` (Fe I 70 eV 8.12 -> 0.082 Mb, Mg II 2p
70 eV 5.96 -> 0.067 Mb, Ca II 3p 45 eV 1.31 -> 0.067 Mb, Cl I-IV -38 to
-45%), which v2 keeps in its opacity; audit E2. Audit Section 5 item 3: with
43 stages every metal input was refused at the default 100 eV band; audit
E1. State note: the `phfit2` sum is restored, `metal_inner_shell_absorption`
is refused, and the 77-stage roster is in production. The linked "current
result" is itself superseded by the restoration results above.)

**Status: interim source-path assessment, updated October 11, 2026; not completion of the atomic-stage expansion.** This report compares the local `~/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2` tree with the currently invoked MoCHII paths. It does not assign MOCASSIN 3.0 behavior to v2, infer active species from array capacity, or claim that proposed MoCHII physics is implemented. Dust, grains, and wider radiative-transfer features are outside this atomic audit unless an invoked path below requires them.

## What MOCASSIN v2 actually calculates in the inspected path

`source/set_input_mod.f90:82,176-186` sets `nstages=7` by default and accepts an override no larger than 10. The cell update is called at `source/iteration_mod.f90:808`; its `updateCell` calls `ionBalance` at `source/update_mod.f90:365`. The balance forms adjacent rate ratios through `min(Z,nstages-1)` (`:1457-1499`), their products (`:1505-1517`), and normalized fractions through `min(Z+1,nstages)` (`:1530-1543`). The last represented stage has no rate to a higher unrepresented stage in this calculation. This is balance **capacity** for an enabled element, not proof that every stage has cooling or line data or that a particular input enabled the element. For example, the stored HII20 benchmark sets `nstages 6` (`benchmarks/gas/HII20/input.in`); its actual scope differs from the default.

The v2 opacity path uses individual metal shells: `source/ph_mod.f90:396,789-822` constructs shell cross sections with `phFitEl`, which implements the Verner PH1/PH2 prescription (`:506-580`); `source/ionization_mod.f90:409-438` adds the stored shell entries to continuum opacity. The invoked metal **rate** path is narrower: `source/update_mod.f90:180-218,248-268` obtains `outShell` and integrates the cross section at that shell's pointer. Its thermal path also selects `outShell` (`:1127-1173`) and integrates the excess above the outer ionization threshold (`:1200-1229`). Consequently, shell-resolved opacity does not establish shell-resolved ionization rates or Auger heating in v2.

`source/ph_mod.f90:276-278` invokes `makeCollIonData` and `makeAugerData`; the latter reads yield arrays at `source/hydro_mod.f90:36-54`. The invoked `ionBalance` does not use these Auger arrays, and `source/update_mod.f90:1485` sets `collIon=0` for elements above H before constructing its ratios. Hydrogen has a collisional coefficient at `:1369-1379`, but the ratio numerator at `:1487-1499` adds that coefficient to photoionization rates without the electron-density factor, so its dimensional form is inconsistent. The H charge-exchange fit is set to zero outside 6,000–50,000 K at `:1462-1467`. Two Ca-labeled coefficient rows overwrite Ar indices at `:1440-1443`. These are observed defects in preserved external source and limit quantitative use of the corresponding v2 results; they were not changed here.

(Corrected 2026-10-11, audit Section 5 item 9: six further v2 defects belong
in this list. Line cooling sums levels 1-10 only (`update_mod.f90:1101-1111`);
Fe II, solved into `forbiddenLinesLarge` (`:1279-1293`), never enters the
cooling; there is no DR for He; elements without DR data take the oxygen
DR; the He+ recombination cooling uses the hydrogen coefficients; and the
free-free opacity loop fills only bin 1 (`ionization_mod.f90:369-392`). The
three cited locations were read for this note; the other three items are
from the audit (Sections 1 and 4). See also the
[bidirectional comparison](MOCASSIN_V2_MOCHII_BIDIRECTIONAL_COMPARISON.md),
"Errors in the preserved external reference".)

The v2 file-name reader looks through `min(Z+1,10)` at `source/hydro_mod.f90:511-518`, and the later line-data availability flag depends on whether a file opens (`:701-726`). Neither stage capacity nor an open file verifies a scientifically complete line atom. The inspected `source/` invocation establishes intended calculation paths; no v2 executable/source identity check or converged v2 model was run for this report.

## MoCHII scope before the shell-selection correction (historical)

Direct reading of the eleven production `data/atomic/element_*.txt` `ELEMENT` records gives **43 metal stages and 32 adjacent transitions**. The revised 77-stage, 66-transition roster is prepared in staged inputs and is not production input; it supersedes the earlier 64-stage target. `src/define.f90:92-93` permits seven stages for eleven elements, but this storage capacity does not activate stages. `src/main.f90:95` invokes the terminal photon-band check; `src/species_mod.f90:225-249` tests the currently loaded terminal potentials, so the current 43-stage roster rejects a 100 eV upper photon edge (that is, every input containing metals at the default band was refused; audit Section 5 item 3, E1; state note: the 77-stage roster is in production since 2026-10-11 and the check refuses a band above Cl VII 114.2 eV). The proposed roster passes this photon test according to the direct threshold assessment, but the test does not evaluate terminal electron-impact ionization or Auger destinations.

MoCHII's metal opacity and rates use the same `species_cross_section` call at `src/species_mod.f90:417-434,1424-1469`: PH2 outer fits below `einn`, then the sum of occupied PH1 shells above it. The packet scorer integrates those cross sections into rates and shell-edge primary-electron heating (`src/ion_score_mod.f90:220-226`; `src/species_mod.f90:491-506,998-1017`). Thus MoCHII computes a shell-summed metal **photoionization rate** at active inner edges, which is broader than the outer-shell rate used by the inspected v2 cell update. Its current balance still assigns that total rate only to the next stage (`src/species_mod.f90:1475-1497`), however. The adjacent solver is used in both direct and cached metal electron calculations (`:1090-1107`; `src/ion_stage_equilibrium_mod.f90:1-12`). No shell-dependent multiple-electron branch is invoked.

The current metal thermal path adds `metal_heating` and subtracts `metal_ionization_cooling` through `src/thermal_mod.f90:82-101`. `src/species_mod.f90:1266-1283,1382-1416` includes the scored primary photoelectron, adjacent collisional-ionization potential loss, and a recombination kinetic-energy term. Its `species_cross_section` comment and formula at `:1424-1469` identify the missing Auger energy and fluorescence distinction. The secondary-electron estimator scores H I, He I, and He II primary electrons only (`src/ion_score_mod.f90:180-196`; `src/gas_rates_mod.f90:126-155`), so metal Auger electrons are not passed through that partition. Existing metal CI, RR/DR, charge exchange, and cooling calls exist, but their physical scope is tied to the **loaded 43-stage roster**, available coefficients, enabled options, and model temperature. The full 77-stage production cooling and photoionization calculation has not been activated. Staged rate consumers passed their scoped 198-case validation, which does not establish complete production line cooling.

| Process | Inspected MOCASSIN v2 invoked calculation | Current MoCHII invoked calculation | Remaining requirement for the proposed 100 eV roster |
|---|---|---|---|
| Metal shell opacity | Individual shell cross sections added | PH2/PH1 switched cross section used by packet opacity | Check every added stage and shell against source fit |
| Metal photoionization rate | Outer-shell cross section in cell update | Sum of active PH1 shells above inner edge | Resolve shell rate into product charge states |
| Metal ion balance | Adjacent rate ratio and product | Adjacent logarithmic stationary balance | Nonadjacent matrix with validated Auger branches |
| Metal electron-impact ionization | Explicitly zeroed in v2 balance | Voronov and HYBRID-adjusted rates only when `metal_collisional_ionization=.true.`; default off (state note: default on since 2026-10-11) | Validate terminal loss or enforce a justified T-domain truncation |
| Metal photoelectron heat | Excess above outer threshold | Primary excess above each active shell edge | Include Auger electron, fluorescence, excitation, and secondary-energy channels |
| Auger data | Yield table read; no invoked balance use | No yield table read into balance | Map physical shell indices and energy-dependent permissible products |
| Line and cooling data | File existence and active-stage dependent | Loaded n-level and cooling data for available stages | Complete and validate all proposed stages and matched profiles |

The [Auger audit](auger_ionization_energy_audit.md) measured nine energetically reachable +2 channels below 100 eV using local Verner PH1 thresholds and Cloudy's reconstructed Kaastra and Mewe yields (Mg, Fe, and Ca: three each). All nine products fit the **proposed** roster; none is represented as a nonadjacent transition in production. Cloudy's yield file is a data source and its shell rows are not an exact copy of the published Table 2. The article's mean Auger energy does not supply a validated, branch-resolved thermal ledger when combined directly with MoCHII's different edge and potential tables. The [terminal CI assessment](terminal_ionization_assessment.md) evaluated the staged 53 internal HYBRID-adjusted transitions plus one additional transition for each element at 801 temperatures from 1,000 to 100,000 K, without photons or hydrogen charge exchange. It is neither a runtime terminal-rate implementation nor a temperature-domain acceptance test. That diagnostic measured a maximum added Mg V fraction of 2.377e-4 and an added terminal CI energy fraction of 6.152e-3. Terminal HYBRID factors return unity, and source-validity limitations remain explicit in that assessment. The numerical size does not establish permission to discard a physically possible loss.

## How the completion report will distinguish the physics

The final report will separate three categories rather than treating every newer table as a new physical process:

1. **Physics shared with the inspected v2 calculation.** Adjacent metal photoionization and recombination equilibrium, charge exchange where an invoked coefficient exists, and line excitation and radiative cooling where an active atom is available. The invoked v2 recombination path explicitly adds `radRecFit` and `diRec` in `source/update_mod.f90:1632-1649`. Updated energies, RR/DR coefficients, and larger line atoms improve data or representation within these process classes; they are not automatically new physical processes.
2. **Physics calculated more broadly by MoCHII, or adopted approximations awaiting activation.** Metal electron-impact ionization and its ionization-energy cooling are available in MoCHII only with explicit `metal_collisional_ionization=.true.`; the default is off. They are absent from the inspected v2 metal balance. Shell-summed metal photoionization rates also exceed that v2 rate path. The staged high-charge H(1s) charge-transfer prescription assumes unit-probability Langevin capture; it is an adopted upper-limit approximation, not a measured reaction rate. Its mean kinetic release uses the exit-channel expression `2.86*(q-1)` eV. The adopted constant lowest-temperature Upsilon continuation changes the extrapolation assumption only; Maxwellian rates still use the actual temperature. The Cl VII and Fe VII Stout tables start at 9,800 K, and accuracy below that source temperature remains unverified. Production activation and numerical fitting acceptance are separate from physical validity.
3. **Physics still required for the proposed domain.** Nonadjacent Auger products, their electron/photon energy partition, and outgoing terminal collisional ionization are incomplete. They must not appear as completed extensions until their actual runtime paths and appropriate diagnostics pass. The revised 77-stage roster alone does not close the network against these processes.

This atomic report does not claim that MoCHII reproduces the entire MOCASSIN v2 code. Dust physics, three-dimensional radiative transfer, line transfer, and other processes need their own source-path and input-matched comparison before such a claim.

## Validation limits and next decision

The MoCHII stage counts were reread from production element records, and v2 claims above were checked against declarations, call sites, and the invoked calculations rather than inferred from plans alone. The Auger script was run and its nine-channel JSON checked; the terminal CI numbers cited here are identified as **previous direct source-fit diagnostic results**, not remeasured in this report. No matched MOCASSIN/MoCHII input set, full temperature/velocity/density/composition profiles, converged model, or numerical agreement claim exists. Gravitational potential, domain, spectral shape and integral, boundaries, active elements, reactions, and energy terms must be matched and documented before a model comparison.

The remaining physical work is to preserve the physical shell index in atomic input, establish energy-consistent product probabilities and electron/photon energy accounting, implement and validate nonadjacent equilibrium, and resolve terminal electron-impact losses over an explicit temperature domain. The atomic input layout and any public control for the accepted temperature domain are costly-to-reverse choices and require confirmation before changing them. This report marks those items as incomplete; it does not convert the staged 77-stage capacity into a claim of fully computed physics.

## Newly verified excitation details

The inspected MOCASSIN v2 `oiii.dat` already uses Storey and Zeippen (2000)
Table 3 optical A values 6.791e-3 and 2.046e-2 s^-1. MoCHII's existing
six-level O III atom previously retained older values despite the CHIANTI
source attribution. The earlier six-level diagnostic measured 5007/4959 = 2.98396833
instead of the old 2.89105976. A new direct evaluation of the currently
installed 177-level atom measures 2.9839683327599866--2.983968332759988 in
twelve temperature/density combinations. Its actual energy-derived wavelengths
are vacuum wavelengths; 5007/4959 are conventional optical labels. See
[the current native measurement](validation_2026-10-11/current_oiii_optical_ratio/README.md). This is a correction
within spontaneous line emission shared with v2, not an added physical
process or a complete O III atom comparison.

C IV and Si V cooling fits and full 24-state and 169-state source atoms
are now installed in MoCHII. Their cold Upsilon continuation is an explicitly
adopted approximation. Their source-node, population, fit, and registered
cooling checks do not establish spatial agreement with MOCASSIN.

For Mg III and Si I, [the actual v2 availability assessment](mg_iii_si_i_mocassin_scope.md)
distinguishes a filename in `fileNames.dat` from a successfully opened and
invoked line atom. For Cl V and Ca I,
[the source assessment](remaining_chlorine_calcium_excitation_sources.md)
likewise separates source tables from an actual line-cooling calculation.

Ar II's published ground fine-structure cooling is now installed after
source-limited table registration passed actual native checks. This does
not represent the full Ar II atom or certify a global thermal temperature
domain. Direct file counting now gives 39 ordinary excitation fits in `data/atomic/`,
including H I and 38 metal fits. File presence does not mean every fit is
active: stage ranges and explicit physical options still select consumers.
The 77-stage production activation, nonadjacent Auger production path, and matched
nebular comparison remain incomplete.


Corrected full Mg III (167 states) and Ne I (31 states, two missing
spontaneous decays restored) cooling are now installed and native-validated.
These electron-excitation models exceed the local v2 line-data availability:
Mg III is listed but its numerical file is absent, and no Ne I collision file
is present. Larger level models remain the same physical excitation process;
this difference does not establish recombination-feeding or full nebular
agreement. Their adopted cold continuations and source-accuracy limitations
are stated at the actual generator and in the production headers.

## Current defaults and remaining shared-model selection

`metal_freefree=.true.` is the user-requested exception. Extra-ion line cooling,
metal collisional ionization, neutral-H excitation, high-charge Langevin
transfer, charge-transfer kinetic heating, and He I metastable coupling are
disabled by default and require their respective explicit input controls.
Ground-state helium charge transfer remains enabled shared physics. The
public controls for larger shared-ion excitation models and consistent shell
selection are not yet implemented. (Corrected 2026-10-11, audit Section 5
item 5: both controls existed in this tree, `metal_excitation_model`
(`'mocassin_v2'` / `'extended'`) and `metal_inner_shell_absorption`. State
note: `metal_collisional_ionization`, `neutral_hydrogen_excitation`,
`charge_exchange_heating` and `hei_metastable_coupling` default to `.true.`
since the restoration; `additional_line_cooling`,
`metal_inner_shell_absorption` and `high_charge_exchange` are refused, the
last renamed `langevin_charge_exchange` (default `.false.`); the default
excitation model is `'chianti'`.) Twelve installed shared-ion models are
larger than the supplied v2 files; ten belong to active production chains.
See [the measured model-selection review](SHARED_EXCITATION_MODEL_SELECTION_REVIEW_2026-10-11.md).

The original v2 archive holds endpoint collision strengths constant outside
multi-node tables. Current candidate upper-temperature rejection therefore
differs from that v2 assumption. This matters for C I, S III, Ar VI, and
Mg I below 100000 K; no hot approximation is silently adopted. (State note,
2026-10-11: the restoration removed the stop; tabulated rows hold Upsilon at
their end node, as HEAD and v2 do, and setup prints the atoms whose tables
end below `te_max`.) See
[the original source review](validation_2026-10-11/v2_collision_temperature_endpoints/README.md).
The current [bidirectional comparison](MOCASSIN_V2_MOCHII_BIDIRECTIONAL_COMPARISON.md)
records the paired six-ion candidates, their completed actual consumers, and
the pending production replacement. Earlier table sizes and source-path line
numbers in this report are historical locations; routine names and current
invocation evidence control implementation claims.
