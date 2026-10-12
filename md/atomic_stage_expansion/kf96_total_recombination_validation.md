# Published Total Charge-transfer Recombination: Further Validation

Date: October 10, 2026 (KST).

## Physical decisions and implementation

N4 now uses the original KF96 Table 1 total fit over 1000--300000 K.
Cloudy's changed fit has unidentified provenance and is not silently preferred
because it is newer. This adds one usable boundary to selective generation.
The validated staging roster has 53 metal stages and 42 adjacent transitions;
it is not activated in production.

Four existing Cl3/4 and Ca3/4 recombination reactions were incorrectly absent
from production. KF96 supplies their total rates and KF99 Table 1 supplies
mean gas kinetic heat. The generator and the two production element files now
include them. Ca uses the printed 1000--30000 K Table 1 intervals. Cl uses
1000--30000 K within the Table 2 calculation span (log10(T/K)=3,3.5,4,4.5),
which is a selected calculation domain rather than a printed fit interval.
The Landau--Zener rate estimates and the ground-reactant, low-density
mean-heat approximation are identified next to the generated records.
Endpoint holding outside these domains asserts no physical accuracy.
No reverse excited-product thermal rate is inferred.

Fe4 total-rate and mean-heat entries also exist in the source helper. Its
Table 2 estimate carries a radiative-charge-transfer caveat, and Fe stage
generation remains blocked by the separate recombination review. It was not
included in the compiled staged-model comparison.

Source equations, coefficients, qualifications, and original papers are in
[the reaction decision](remaining_charge_transfer_decision.md). New papers
and the author Fe DR archive are stored in `~/MoCHII/references/`; project
reference links preserve historical controller paths. No numerical products
were regenerated for this reference relocation.

## Matched comparison conditions

| Condition | Reference | Compiled calculation |
|---|---|---|
| Rates | Nine previously accepted rows plus N4 and four Cl/Ca rows | Same fourteen serialized tables |
| Temperature | 2001 logarithmic samples inside each selected domain | Same samples plus two endpoint-continuation samples |
| Mean heat | Independently transcribed KF99 Table 1 constants | Actual `CXHEAT` reader and heat consumer |
| Plasma checks | Prescribed T=5000,10000,30000 K; n_H=100, n_e=10, n_HI=n_HII=50 cm^-3 | Same fixed inputs; each metal abundance 1e-4 and each photoionization rate 2e-8 s^-1 |
| Radiation transport, gravity, spatial boundaries | Not evaluated | Controlled local rate tests, without a nebular solution |
| Product branching and diffuse cascade absorption | Outside the adopted mean-heat treatment | Not introduced |

## Measured validation

The finite tmux controller `atomic_kf96_totals_20261010` completed Python
checks, staged generation, compiled evaluation, and figure generation.
Its [state and hashes](validation_2026-10-10/kf96_total_recombination/controller/state.json)
and products are retained.

- 31 fitting Python checks passed before the subsequent independent audits.
- 28,042 compiled rate samples passed. The maximum relative interpolation
  error was `5.331534752e-6` (Ca3), below the `1e-5` requirement.
- 42 compiled mean-heat samples matched the adopted values exactly.
- Normalized finite ion fractions, electron accounting, reaction-resolved
  hydrogen transfer, and rate-weighted heat passed at all three temperatures.
- After correction, the actual 43-stage, 32-transition production roster passed
  the same local charge and heat consumers. Four production reactions passed
  8,012 compiled coefficient samples and twelve heat samples. Its maximum
  relative rate error was again `5.331534752e-6`.

The [production result](validation_2026-10-10/kf96_total_recombination/production_correction/result.json)
records hashes. Original Cl/Ca files are preserved beside that result as
`element_cl.txt.before` and `element_ca.txt.before`. Other transition numerical
records in those files were checked against staging before replacement;
only `CXR`, `CXRTAB`, and `CXHEAT` in transitions 3 and 4 changed numerically.
Other production element files were not replaced.

Full curves are in [JSON](validation_2026-10-10/kf96_total_recombination/products/curves.json)
and the visually inspected [vector figure](validation_2026-10-10/kf96_total_recombination/reviewed_plot/added_charge_transfer_heat.pdf).
The original figure said all intervals were published; a plot-only correction
now states selected intervals. The original PDF and numerical products remain
in the controller products. This correction did not recompute rates.

No full regression suite or nebular profiles were calculated. The affected
path is element generation, coefficient reading, ion fractions, hydrogen
charge transfer, and gas kinetic heating. Previous nebular products remain
preserved and are stale for the corrected Cl/Ca reactions.

## Newly established requirements before activation

The final Ferland et al. (1997) paper calls its high-charge Langevin rule a
capture upper limit. The Fe-specific rounded coefficient must not be treated
as a verified species-specific transfer rate for every element. The new
`charge_exchange_capture.py` diagnostic uses the correct reduced mass and
H(1s) polarization potential; it is not connected to production generation.
Its four focused checks passed. The ten proposed high-charge boundaries and
classical barrier diagnostics are recorded in
[capture bounds](validation_2026-10-10/kf96_total_recombination/capture_bounds.json).
No temperature interval or unit reaction probability is falsely certified.

The Fe V--VII source investigation found a 2026 author dataset, but provisional
fit transcriptions disagree with that dataset at some nodes. Three direct
ADF09 parser checks passed; this is not production fit acceptance. See
[Fe source decision](fe_recombination_source_decision.md).

N IV's conflicting observed energies are authentic NIST entries, and its
computed transition identities require an ion-specific energy-model decision.
The n-level generator now refuses the known inverted N IV radiative rows
rather than silently moving Einstein A to a different emitting state. One
focused preservation/rejection test passed. Its existing low-density cooling
guard remains. See [N IV source decision](niv_atomic_source_decision.md).

Production still has 43 stages. Full 64-stage activation is not physically
accepted, and no controller remains running.
