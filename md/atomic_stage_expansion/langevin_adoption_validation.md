# Adopted high-charge hydrogen transfer approximation

Date: October 10, 2026. The user explicitly adopted the Langevin capture
upper-limit approximation. This resolves the previously pending rate choice
for the ten selected q >= 5 boundaries; it does not establish measured
species-specific electron-transfer probabilities.

## Physical definition and scope

For X^(q+) + H(1s) -> X^((q-1)+) + H+, use
k = 2 pi sqrt(2 C4 / mu), where U(R) = -C4/R^4 and mu is the actual
ion--hydrogen reduced mass. The H(1s) polarizability is 4.5 a0^3.
The model assumes classical trajectories and unit electron-transfer probability
for every trajectory crossing the centrifugal barrier. It is temperature
independent; neither the final Ferland et al. (1997) article nor the numerical
sampling below establishes a temperature interval of physical accuracy.
Existing measured/published KF96 fits are retained at lower charges.

The separate thermal prescription uses the KF99 9.5-bohr avoided-crossing
model with the physical exit-channel charge: mean kinetic release =
2.86(q-1) eV. The initial charge is q, while the products have charges +1
and +(q-1); their Coulomb energy is (q-1)/9.5 hartree. The final KF99 prose
prints 2.86q, which does not agree with that potential or the repeated
lower-charge Table 1 estimates. Cloudy c25 uses the consistent q-1 index.
The initial staging attempt copied the literal KF99 prose; this confirmed
physical inconsistency was corrected before production activation. It would
otherwise overstate each selected reaction's mean heat by 2.86 eV (25% at
q=5, 20% at q=6 relative to the corrected estimate). Capture rates are unchanged.
The approximation assumes ground-state reactants and low-density radiative
relaxation. It excludes cascade-photon reabsorption and product collisional
redistribution and does not establish species-specific product populations.
No state-specific reverse rate is inferred: CXI remains zero at these boundaries.

| Element | Initial charge q | Computed k (cm^3 s^-1) | Mean kinetic release (eV) |
|---|---:|---:|---:|
| N | 5 | 9.861879492e-09 | 11.44 |
| S | 5 | 9.673682519e-09 | 11.44 |
| S | 6 | 1.160841902e-08 | 14.30 |
| AR | 5 | 9.644553455e-09 | 11.44 |
| AR | 6 | 1.157346415e-08 | 14.30 |
| FE | 5 | 9.610714650e-09 | 11.44 |
| FE | 6 | 1.153285758e-08 | 14.30 |
| CL | 5 | 9.659586950e-09 | 11.44 |
| CL | 6 | 1.159150434e-08 | 14.30 |
| CA | 5 | 9.644168538e-09 | 11.44 |

## Implementation and measured validation

`charge_exchange_capture.py` computes the reduced-mass bound and selected
avoided-crossing heat. `make_element_data.py` emits constant form-1 CXR coefficients and
CXHEAT records with their approximation limits. `validate_langevin_expansion.py`
prepares all 11 elements, 64 stages and 53 adjacent transitions in staging.
The generator still rejects boundaries outside this reviewed target.

The actual linked debug Fortran reader, reaction rates, direct/cached fractions,
electron charge, reaction sums, and mean kinetic heat consumers passed.
Across ten reactions and 101 temperatures (10--10^7 K), 1,010 samples have
maximum relative rate difference 3.743798548949684e-7 from the direct bound.
Thirty heat checks at 5,000, 10,000, and 30,000 K have maximum difference
1.7763568394002505e-15 eV. This sampling checks numerical implementation and
temperature independence, not physical accuracy. Thirteen focused Python tests
for capture, generation safety, and published lower-charge totals pass.

Evidence: [native result](validation_2026-10-10/langevin64_exit_heat/products/native_result.json)
and [controller state](validation_2026-10-10/langevin64_exit_heat/controller/state.json).
Production still declares 43 metal stages and 32 transitions. No expanded
production activation or nebular comparison has occurred. Missing cooling,
terminal CI and nonadjacent Auger accounting remain required before activation.

Final source provenance and the original lower-charge audit remain in
[charge-transfer sources](remaining_charge_transfer_decision.md).

The previous `langevin64` products and their controller remain preserved as a
superseded diagnostic using the inconsistent literal 2.86q heat. Do not use
those heat records as accepted outputs; `langevin64_exit_heat` replaces them.
