# Auger implementation requirements for the 100 eV ion roster

## Physical judgment

The current target is 77 metal stages. The earlier 64-stage audit already
contains the destinations of the nine nonadjacent
photoionization channels that are energetically possible below 100 eV in the
current Verner PH1 and NIST threshold audit. The source material supports a
shell-resolved **charge-channel diagnostic**. It does not yet support an
energy-conserving thermal implementation. Kaastra and Mewe (1993) publish
electron-count probabilities and a mean Auger energy for each initial vacancy,
plus mean fluorescence line yields. For the currently invoked metal heating
calculation, all emitted metal-electron kinetic energy is thermalized: a
consistent unconditional mean is sufficient for the aggregate linear heating
integral. A joint distribution of final charge and electron energy is not a
prerequisite for that restricted calculation. It is needed for conditional
charge/energy outputs; nonlinear secondary deposition requires at least the
emitted-electron energy spectrum. Charge yields must still be physically valid
at each photon energy. KM primary-vacancy energies differ from the PH1 edges,
and their published mean-energy correction is negative for Ca II M1. Keep
unsupported inner-shell physics out of production until those actual energy
and onset problems are resolved. Do not assign unaccounted excitation or
fluorescence energy to heat.

The [direct mean-energy review](validation_2026-10-11/auger_mean_energy/README.md)
supersedes the earlier blanket requirement for joint branch energies. It does
not establish a usable heating table or activate a production calculation.

## Source facts and correction to the earlier audit

The final journal article is J. S. Kaastra and R. Mewe (1993), *Astronomy and
Astrophysics Supplement Series* **97**, 443–482, bibcode
`1993A&AS...97..443K`. The local final-article scan is
`~/MoCHII/references/Kaastra_Mewe_1993_AAS_97_443.pdf` (copied byte-for-byte from the existing local final-journal scan, SHA-256 `4b22f56f4c38a5f5dc24cd553eb49b17a7b40c87364e9754a635577e79e02cc5`). The 40-page scan is raster-based; `pdftotext -layout` exposes its ADS bibcode but not the article text, so the earlier page inspection remains necessary.
Its journal page 446, Section 2.2, describes a branching cascade that removes
an Auger transition when its electron energy is nonpositive. Journal page 447,
Section 3 and Eqs. (5)–(7), compares mean Auger electron energy, mean photon
energy, and successive outer-electron binding energies. It reports energy
discrepancies usually within 10% but up to 70% for some M1 vacancies,
partly because of final-state excitation and subshell correlation. The paper
offers a mean-energy correction; it does not publish the joint branch energies
needed for nonlinear secondary-electron deposition.

The deposited CDS `ReadMe`, `table2`, and `table3` are at
`~/RT_Codes/ExoAtmosphere/references/kaastra_mewe_1993/`. Table 2 defines
the initial ion stage, KM vacancy subshell, primary-vacancy energy, mean Auger
energy, energy discrepancy, and probabilities for one through ten emitted
electrons (including the primary). Table 3 gives mean fluorescence line yields
and approximate photon energies. Its `Delta` field counts Auger electrons
emitted before a line, not the final electron count for that cascade. Neither
table gives an event-by-event electron spectrum or a joint distribution of
final charge, electron energy, fluorescence energy, and residual ion excitation.

The earlier `auger_ionization_energy_audit.md` incorrectly stated that CDS
Table 2 lacked Mg I 2s and Ca I 3s rows. Direct fixed-column reads find
Mg I `(Z,stage,KM shell)=(12,1,2)` and Ca I `(20,1,5)`, both with two-electron
probability 1.0000. KM shells K/L1/L2/L3/M1/M2/M3 map to PH1 1s/2s/2p/2p/
3s/3p/3p, respectively. The corrected
[audit JSON](validation_2026-10-10/auger_audit_cds_checked/auger_channels.json)
checks all nine reachable two-electron channels against the deposited CDS
rows; all nine mappings are complete. The greatest difference from Cloudy's
transcribed probabilities is 0.00001, consistent with the last stored digit.
The original [audit JSON](validation_2026-10-10/auger_audit/auger_channels.json)
is preserved. Cloudy's `mewe_nelectron.dat` header itself warns that its shell
resolution differs from the KM convention; this requires an explicit map, not
an assumption that numeric shell indices match. These are calculated cascade
probabilities, not independently measured yields.

Verner's `~/RT_Codes/Verner_Fortran_sub/phfit2.f:1-23,35-59` uses physical
shell indices 1s=1 through 4s=7. MoCHII's `PHSHELL` rows currently retain
only orbital angular momentum and fit coefficients
(`tools/fitting/make_element_data.py:715-761,905-916` and
`src/species_mod.f90:85-96,1424-1469`), so the original physical shell
index cannot safely be recovered from the row order. The 100 eV audit found
31 active PH1 rows and nine reachable two-electron channels, from Mg, Fe,
and Ca. It found no reachable product beyond the approved 64-ion roster and
no reachable third-electron channel. This is a threshold test, not a
spectrum-weighted rate.

| Source and shell | PH1 edge (eV) | KM edge (eV) | Sequential product threshold (eV) | CDS two-electron probability |
|---|---:|---:|---:|---:|
| Mg I 2s | 94.00 | 92.0 | 22.682 | 1.0000 |
| Mg I 2p | 54.90 | 54.3 / 54.0 | 22.682 | 0.9967 |
| Mg II 2s | 98.81 | 104.0 | 95.179 | 0.9999 |
| Fe I 3p | 66.00 | 60.9 / 59.0 | 24.102 | 0.9999 |
| Fe II 3p | 76.17 | 82.9 / 81.0 | 46.851 | 0.9999 |
| Fe III 3p | 87.05 | 104.9 / 103.0 | 85.564 | 0.9999 |
| Ca I 3s | 48.30 | 46.0 | 17.985 | 1.0000 |
| Ca I 3p | 34.43 | 28.4 / 28.0 | 17.985 | 0.9950 |
| Ca II 3s | 60.37 | 58.0 | 62.785 | 0.9999 |

The table shows two distinct threshold problems. Ca II 3s has a product
threshold 2.415 eV above its PH1 edge, so its two-electron probability cannot
apply immediately at that edge. Fe III 3p has CDS vacancy thresholds above
100 eV while its PH1 cross section starts at 87.05 eV. The source does not
state a photon-energy-dependent yield that bridges these intervals.

There is also a direct energy contradiction if the tabulated mean Auger energy
is combined unchanged with the PH1 primary photoelectron and NIST final-ion
binding energies. For Mg II 2s, `95.179 + 20.1 - 98.81 = 16.469 eV`; for Ca II
3s, `62.785 + 7.1 - 60.37 = 9.515 eV`. These positive deficits occur before
adding any fluorescence or residual excitation. At the channel onset, even
`E - E_PH1` for the primary electron plus those two terms exceeds the absorbed
photon energy. This demonstrates incompatibility of the unmodified averages
and thresholds, not a measured 16.469 or 9.515 eV physical energy error.

## Current invoked path and required changes

`src/species_mod.f90:1445-1469` sums PH1 shell absorption above `PHINNER`,
but its `heat` value retains only the primary electron excess over each shell
edge. `src/ion_score_mod.f90:220-226` calls
`species_rate_path_score`, which scores one total rate and one heat integral
for each adjacent transition (`src/species_mod.f90:491-506`).
`species_gamma_compute` converts them to `egam/eheat` (`:994-1017`). Both
the direct fraction calculation and `species_metal_sums` then call the
adjacent-rate solver; `metal_heating` consumes `eheat` (`:1266-1283`). The
secondary-electron estimator in `src/ion_score_mod.f90:180-196` only scores
H I, He I, and He II primary electrons. Thus correct nonadjacent products
would change the fractions, free electrons, opacity feedback, line and free-free
emission, and the thermal calculation.

After a consistent mean-energy model and applicable charge yields are obtained,
or a justified state-resolved cascade model is obtained, the
implementation needs these coordinated changes:

1. Version the atomic input and retain each PH1 physical shell index, its
   cross-section coefficients, an explicit mapping to the deposited KM
   subshell rows, and the adopted threshold convention. Validate the mapping
   and all active one- and two-electron probabilities. A change to this data
   layout and to public output meaning requires prior confirmation.
2. At each packet energy, resolve the PH1 shell cross sections. Below
   `PHINNER`, retain the distinct PH2 outer fit with its stated range. Above
   `PHINNER`, assign each shell absorption to represented final stages using
   physically applicable energy-dependent branch probabilities. Require
   `E` to exceed the sequential final-ion threshold; do not renormalize
   forbidden empirical branches without a conditional cascade calculation.
   Add shell and branch onset energies to the transport band's edge list.
3. Score photon absorption by source and destination stage, preserving the
   sum of branch rates as the shell absorption rate. Build an element-wise
   stationary rate matrix with the existing adjacent CI, recombination, and
   charge transfer and the new nonadjacent photoionization entries. Solve the
   normalized nonnegative stationary populations with a check for multiple
   closed classes. Use the same fractions in the electron closure, thermal
   terms, packet opacity, and output. The current adjacent solver must be
   recovered when all nonadjacent rates vanish.
4. For a state-resolved calculation, close each absorbed-photon ledger:
   `E = sum(IP) + K_primary + sum(K_Auger) + E_fluorescence + E_excitation`.
   Deposit each electron's kinetic energy through the applicable thermal and
   secondary-ionization/excitation calculation. Account for fluorescence as
   escaping or transported radiation and for residual excitation through an
   identified decay channel. Reject negative residual energy and any event
   that exceeds its photon energy. For the existing full thermalization of
   metal electrons, an unconditional mean instead gives
   `Q_mean(E,s) = E - I_s + <EA_s>`. Require a consistent mean ledger for
   fluorescence, final binding, and residual excitation; this formula alone
   does not prove conservation. Do not copy the same unconditional mean into
   every charge branch. If nonlinear metal secondary deposition is later
   enabled, integrate the deposition law over the electron spectrum rather
   than evaluating it at the mean total kinetic energy.

The charge-path matrix can be built and tested in an isolated diagnostic with
the source probabilities. It cannot be activated as a physically complete
photoionization and thermal path until photon-energy-dependent onsets and
a consistent mean-energy model are supplied; a state-resolved cascade is one
way to establish them. MOCASSIN v2 (`source/ph_mod.f90:278`, `hydro_mod.f90:36-54`,
`update_mod.f90:1459-1543`) reads Auger data but its invoked ion balance
uses adjacent ratios. MOCASSIN 3.0 (`source/update_mod.f90:1833-1855`) adds
yield-weighted terms with a population-ratio approximation. Cloudy c25.00
(`source/ion_solver.cpp:649-668`) uses shell rates and electron-count yields
but clips products at its configured high stage. Its thermal code
(`source/cont_gammas.cpp:382-435`) uses a separate threshold and secondary
approximation. These paths confirm that the source probabilities can be used
for charge accounting; none supplies the missing MoCHII energy distribution.

## Acceptance checks

Check PH1 shell cross sections against the original `phfit2.f` at each shell
edge and through 100 eV. Compare every imported probability and physical
shell key with CDS Table 2, including both KM p subshells; check normalization
at source precision and ensure every positive represented branch has a
documented onset. At Ca II 3s and Fe III 3p, test the discrepant threshold
intervals explicitly. Test matrix column conservation, nonnegative normalized
fractions, the all-adjacent limit, and a known nonzero two-electron rate.
Reject any product outside the active roster instead of redirecting it to the
terminal stage.

For a mean model, test its mean energy ledger and charge onsets without
claiming event-resolved conservation. If branch energy information exists,
also test each one- and two-electron event at onset and at 100 eV, including
fluorescence, excitation, and any secondary deposition. Compare full temperature, velocity,
density, and composition profiles with any external plasma calculation only
after matching the spectrum and its integral, stage limits, chemistry,
boundaries, and enabled processes. The present audit has no matched spectrum
or plasma run. Its proven output is reachability and source-table mapping.
