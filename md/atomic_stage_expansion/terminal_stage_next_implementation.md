# Terminal-stage closure within the 64-ion target

The user's latest decision retains 64 stages and 53 adjacent transitions.
Additional charge states may be used in private diagnostics, but expanding
the production roster is outside this task.

## Physical decision

The selected 64-ion, 53-transition roster closes **direct outer-shell
photoionization by photons at or below 100 eV**: all eleven terminal
ionization thresholds exceed 100 eV, with the smallest at Ca VI (108.784 eV)
and Mg IV next (109.265 eV). It does not close electron-impact ionization.
At every finite electron temperature the terminal ion can lose an electron;
the current adjacent-stage solver omits that outgoing rate and its reverse
recombination because the terminal transition does not exist in the element
file. Inner-shell absorption can also produce nonadjacent Auger products.
These are separate physical limitations, not failures of the logarithmic
equilibrium solver.

The immediate implementation should keep the 64-ion target staged and make
its validity testable. Do not silently set a terminal rate to zero or insert
its ionization energy as cooling while retaining a closed 64-stage balance.
The latter would remove energy without the corresponding particle and reverse
reaction. Do not expand the public stage roster, output arrays, or user
temperature limits without review of that changed interface and its products.

## Evidence from the current code and diagnostic

The staged files at
`validation_2026-10-10/langevin64_exit_heat/products/atomic` contain the
roster C 5, N 6, O 5, Ne 5, S 7, Ar 7, Mg 4, Fe 7, Si 5, Cl 7, and Ca 6:
64 ions and 53 adjacent transitions. The installed `data/atomic` roster is
smaller and is not the subject of this estimate. `main.f90` calls
`species_terminal_photoionization_check` before setting up the ionization
bands; that check rejects a transported upper photon energy above any active
terminal threshold. It does not inspect electron collisions. In
`species_mod.f90`, the loops forming ionization rates, stage fractions, and
ionization-energy cooling stop at `nstage-1`. The packet path scores each
absorbed inner-shell event into the listed adjacent transition.
`ion_stage_equilibrium_mod.f90` correctly solves the supplied adjacent rates,
including exact zero-rate boundaries, but has no nonadjacent transition matrix.

The existing [terminal diagnostic](validation_2026-10-10/terminal_ionization/result.json)
and [complete profiles](validation_2026-10-10/terminal_ionization/full_profiles.npz)
evaluate 801 temperatures from 1,000 to 100,000 K with no photons, fully
ionized hydrogen, no charge transfer, and adjacent electron processes only.
One additional stage is included separately for each element. The same
internal `dere_hybrid` factors are used on both sides of each comparison; no
factor is tabulated for any of the eleven added boundaries, so the diagnostic
uses the source Voronov rate with factor one. Electron density cancels from
these particular equilibrium ratios. The arrays were read directly to check
the following values:

| Temperature | Mg V fraction with the extra stage | Added Mg terminal CI energy as a fraction of the extended CI energy |
|---:|---:|---:|
| 10,000 K | 2.900e-95 | 1.592e-93 |
| 30,026 K | 9.992e-27 | 9.737e-25 |
| 100,000 K | 2.377e-4 | 6.152e-3 |

Mg has the largest measured terminal effect: its maximum added stage fraction
is 2.377e-4, maximum change to an existing stage fraction is 1.941e-4, and
maximum added fraction of CI ionization-energy loss is 0.6152%, all at the
upper sampled temperature. The next largest added energy fraction is Si at
1.064e-6; O is 4.461e-7. These are fractions of the diagnostic's CI energy
loss, not fractions of total nebular cooling. Numerical zero at the cold end
can result from underflow; it is not proof of a physically zero rate. The
complete 801-point arrays, rather than these four temperatures, must be used
when testing a proposed validity boundary.

The reverse rates in this restricted calculation came from local Badnell
M=1 RR/DR records for ten ions and CHIANTI Fe VIII RR/DR records for Fe.
The Fe VIII DR coefficients agree with Table 3 of the final published
[Schmidt et al. (2008) article](https://doi.org/10.1051/0004-6361:200810834),
as checked in `terminal_ionization_assessment.md`. Its experimental plasma
range starts at an electron energy of 0.2 eV, above the 1,000 K diagnostic
endpoint. Published fit-to-calculation intervals for the other ten reverse records
have since been assessed in
[the reverse-rate provenance audit](terminal_reverse_rate_source_validity.md).
O VI and Ne VI DR fits are used below their published 125,000 K accuracy
interval in this diagnostic; the archived source tables permit sampled
residual checks but do not extend that published guarantee. Transition-specific
Voronov errors, terminal H charge-transfer rates, and energy-consistent
nonadjacent Auger yields remain unresolved. Therefore this diagnostic quantifies a restricted electron-only
comparison; it cannot certify that the full 64-ion plasma model is closed
from 1,000 to 100,000 K.

## Next implementation and acceptance conditions

1. Keep the current 64-ion staged files and 100 eV photon check. Add a
   read-only source diagnostic for all eleven terminal boundaries that
   evaluates the local Voronov rate, identified RR+DR reverse rate, terminal
   population estimate, and omitted ionization-energy term across complete
   temperature profiles. Give each source fit its verified temperature range;
   identify the `dere_hybrid` factor-one continuation as uncalibrated rather
   than a measured correction. Include the actual active hydrogen state and
   radiation conditions when a spatial calculation is compared. Do not use
   the electron-only arrays as a substitute for that comparison.
2. Determine terminal H charge-transfer channels and inner-shell Auger
   branching from identified sources. A capture upper limit is not a measured
   charge-transfer rate, and a missing table row is not evidence for zero.
   Compare each channel's loss with the adjacent electron loss under the same
   conditions. For Auger products, retain the product charge and electron
   energy accounting; adding their rate to the adjacent ionization rate would
   give the wrong stage distribution.
3. Decide the physical operating domain from these complete rates and matched
   state profiles. A bounded 64-ion approximation is acceptable only where
   terminal loss, missing reverse flow, and Auger escape satisfy an explicit
   physical accuracy requirement throughout the domain. Record that range at
   the place where the approximation is used and enforce it before consuming
   the stage fractions. The existing 100 eV test establishes only the direct
   photon condition. The current 1,000–100,000 K source diagnostic alone does
   not establish a thermal validity range or a numerical acceptance threshold.
4. Retain the approved 64-stage production roster. If a proposed operating
   domain fails the truncation assessment, record the failing conditions and
   keep that domain unaccepted. A diagnostic with higher charge states must
   include their particle and energy equations; it does not authorize a
   production roster change. Selecting a narrower public temperature range
   or a new approximation also requires an explicit physical justification
   and review of its effect on inputs and outputs.

The existing `make_element_data.py` refuses an added boundary without both
`CX_HIGH` and `CX_HEAT` data, which prevents an unreviewed stage extension.
The next diagnostic can run without changing production files or user inputs.
Neither its completion nor a low terminal fraction should by itself authorize
activating the staged 64-ion roster.
