# Published Mean Charge-transfer Heat: Staged Implementation

Date: October 10, 2026 (KST).

## Physical acceptance

The nine added recombination reactions listed in
[the rate audit](charge_exchange_expansion.md) now have identified published
thermal data. Kingdon and Ferland (1999), ApJ 516, L107--L109,
[DOI 10.1086/312008](https://doi.org/10.1086/312008), Table 1 supplies mean gas
kinetic heat. The original final
[publisher PDF](references/Kingdon_Ferland_1999_ApJ_516_L107.pdf) was downloaded
from IOPscience, read with `pdftotext -layout`, and its table was visually
checked. Its SHA-256 is
`2c681fe32532baddee5d19a67e472cefc51f6a0135de01ab82eb9aa182e78def`.

These constant means can represent a ground-reactant, low-density
approximation. They do not establish individual product terms or temperature
dependent branching. Some entries use avoided-crossing estimates rather
than resolved channel calculations. Collisional redistribution at high
density and absorption of diffuse cascade photons require further treatment.
The fit interpolation errors below measure numerical reproduction of that
adopted model, not its physical uncertainty.

KF96 Table 1 supplies the total rates, whereas KF99 Table 1 supplies the heat.
The earlier statement that these nine heat sources remained unidentified is
superseded by the publisher check. The separate
[product-channel](charge_transfer_kf96_channels.md) and
[Cloudy/MAPPINGS](charge_transfer_cloudy_heat.md) investigations document the
source and limitations. MAPPINGS stores the same heat in units of
`k_B * 10000 K`; its agreement is another transcription, not independent
atomic evidence.

## Implemented behavior

`tools/fitting/charge_exchange_expansion.py` now records nine published
temperature intervals and `KF99_ADDED_HEAT`. Its Cloudy comparison checks
rate coefficients, intervals, and heat independently. Production generation
does not depend on Cloudy's heat column: the published constants supply the
new heat records directly.

`make_element_data.py` can selectively write these boundaries using the
existing `CXR 8`, `CXRTAB`, and `CXHEAT` records. The source and physical
approximation appear next to each new heat record. Recombination rates have
their individual published intervals and adaptive log-log interpolation.
Holding an endpoint outside a fit interval remains an explicit continuation
without asserted physical accuracy. No reverse rate is inferred from a
total forward rate or an assumed product excitation.

The selected staged roster is C I--V, N I--IV, O I--V, Ne I--V, S I--V,
Ar I--V, Mg I--IV, and the existing Fe I--IV, Si I--V, Cl I--V, Ca I--V.
It contains **52 stages and 41 transitions**, nine more transitions than
production. This is an atomic-rate test fixture with controlled ionization
rates. Its missing cooling files are reported. It is not a proposed default
roster or a complete 100 eV plasma model. The existing terminal-band check
would still reject 100 eV for several of its elements.

The production files continue to declare **43 stages and 32 transitions**.
Their SHA-256 hashes were checked before and after every validation step;
none were replaced. No production Fortran routine, user parameter, or file
format changed in this step. Unresolved boundaries and Fe recombination
still cause generation errors. A missing heat record remains an error even
when its total rate exists.

## Matched conditions and numerical checks

| Condition | Published-expression reference | Compiled calculation |
|---|---|---|
| Reactions | Nine KF96 recombination fits | Same nine selected transitions |
| Temperature | 2001 logarithmic points within each fit interval | Same points, plus one outside each endpoint |
| Heat model | KF99 constant mean heat | Same serialized mean in `CXHEAT` |
| Rate model | Direct analytic expression | Serialized adaptive table; log-log interpolation |
| Density, spectrum, gravity, boundaries | Not parameters of the rate-coefficient comparison | Not used in the coefficient comparison |
| Product redistribution and cascade absorption | Outside the adopted mean-heat treatment | Not introduced by these records |

- 28 focused Python tests passed. New checks cover all nine heat values,
  physical energy budgets, full published rate intervals, actual selective
  generation, and rejection/preservation when a heat record is removed.
- 18,027 compiled rate samples passed. The maximum relative interpolation
  error is `5.154022608e-6` (N3), below the `1e-5` requirement. The largest
  new table contains 500 nodes, below the private capacity of 1024.
- Nine mean energies were read and checked at 5000, 10000, and 30000 K:
  27 heat samples, with zero difference from the independently transcribed
  published values. The reaction directions and absence of an inferred
  reverse were checked.
- Controlled fraction tests at those temperatures used `n_H=100`, `n_e=10`,
  `n_HI=n_HII=50 cm^-3`, abundance `1e-4` for each metal, and prescribed
  photoionization rates `2e-8 s^-1` for all included transitions. These are
  held inputs, not a self-consistent transported radiation solution.
  Fractions are finite and normalized. Direct and cached electron sums,
  reaction-resolved hydrogen/metal charge transfer, and rate-weighted thermal
  energy agree within the declared numerical tolerances.

The full thermal-coefficient curves are saved as
[JSON](validation_2026-10-10/kf99_added_mean_heat/corrected_products/curves.json)
and a visually inspected
[vector PDF](validation_2026-10-10/kf99_added_mean_heat/corrected_products/added_charge_transfer_heat.pdf).
The plotted quantity is `k(T) * mean_heat`, in `erg cm^3 s^-1`. It contains
no density, temperature, velocity, or composition profile of a nebula.

The finite controller sequences and source hashes are in
`validation_2026-10-10/kf99_added_mean_heat/`. The initial test used the wrong
index for the reaction-output array and stopped; the actual production API
was correct. Its outputs and recovery note are preserved. Corrected checks
write separate products. A plot-only follow-up removed overlapping minor
temperature labels; the earlier PDF is preserved as `previous_layout.pdf`.
The initial controller stopped on the test failure. Both corrected
validation and the plot-only controller completed; none remains running.

## Remaining physical decisions

Nine selected boundaries now have a usable published mean-heat treatment.
At completion of this nine-reaction step, twelve additions remained.
The subsequent [further validation](kf96_total_recombination_validation.md)
implements original N4 and identifies Fe4's qualified Table 2 total, while
Fe generation remains blocked by recombination review. Ten higher-charge
reactions still require a physical prescription.
Higher-charge capture estimates do not establish product branching or
transfer probabilities. KF99's higher-charge `2.86q` wording and Cloudy's
`2.86(q-1)` evaluation are distinguished: the latter follows the product
Coulomb charges. No external source was edited and no new high-charge rule
was introduced.

Fe V/VI low-temperature DR, N IV energy inconsistencies, cooling convergence,
Cl VII low-temperature collisions, collisional terminal losses, and Auger
ionization still require review. No full regression suite, matched nebular
comparison, production activation, or profile calculation was performed;
this scope tests the newly generated rates and their charge/heat consumers.

## Subsequent high-charge adoption and physical correction

The user adopted the ten high-charge Langevin upper limits after this
lower-charge report. Fe q=3--6 recombination is now verified and implemented.
The separate high-charge mean heat follows the 9.5-bohr exit-channel Coulomb
estimate 2.86(q-1) eV. The literal final KF99 prose 2.86q conflicts with its
initial-charge definition and exit potential. See
[Langevin adoption and new native checks](langevin_adoption_validation.md).
Earlier unresolved-choice statements in this report describe the historical
audit, not a current activation blocker. Cooling and CI/Auger work remain.
