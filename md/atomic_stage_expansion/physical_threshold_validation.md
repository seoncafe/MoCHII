# Physical Ionization Energies and Photon-band Validation

Date: October 10, 2026 (KST).

## Physical judgment and implemented policy

The energy removed from the radiation field and used to ionize an atom must
use its physical ionization potential. An energy parameter fitted to a
collision rate is not interchangeable with that potential. MoCHII now uses
the physical energies in Cloudy's identified NIST ASD snapshot for `ETH`,
the PHOTO cutoff, and the outermost occupied PHSHELL cutoff. Inner PH1 edges,
all cross-section shape parameters, and Voronov CI fit parameters retain
their independent source values. The generator does not shift the photon
energy supplied to the analytic cross-section formula.

This follows the inspected implementation in Cloudy c25.00:
`source/atmdat_adfa.cpp::getEthresh` returns `getIonPot(...)*EVRYD` for
`set_vshell_index(nelec,Z)` and the original PH1 edge for other shells;
`phfit` invokes that threshold check before evaluating the source fits.
`source/opacity_createall.cpp` also prevents the discretized photon energy
from falling below the physical edge. These are implementation checks, not
new verification of the original fit papers.

The source is `~/CLOUDY/c25.00/data/ionization_potentials.dat`, magic number
`20221115`. Its SHA-256 is recorded in the generated
`src/metal_ionization_potential_mod.f90` header and source inventory.
`tools/fitting/ionization_potentials.py` converts wavenumbers using
`hc/e = 1.2398419843320026e-4 eV cm`, checks ion identities, duplicate rows,
finite positive energies, and the snapshot identity, and stops before the
source bibliography. This is an identified local snapshot; its values have
not been compared with a newly queried NIST release.

## Production changes and preserved products

All 32 existing metal transition energies were updated in the eleven
`element_*.txt` files. No additional ion stage was activated: the production
roster remains 43 stages and 32 transitions. Cross-section thresholds and
charge-transfer energy budgets were regenerated consistently. Numeric record
comparison allowed only the physical energies and dependent energy budgets;
CI, RR, DR, charge-transfer rate coefficients, and cross-section shape
parameters were unchanged during this step. The later KF96 temperature-range
correction is recorded separately below.

Representative measured energy changes are:

| Absorbing ion | Previous ETH (eV) | Physical ETH (eV) |
|---|---:|---:|
| S III | 34.830000 | 34.855677706 |
| Cl III | 39.610000 | 39.795208171 |
| Cl IV | 53.470000 | 53.242534333 |

For example, Ar IV's Voronov fit energy remains 52.3 eV while its physical
ionization potential is 59.5818464 eV. That fit energy is preserved when the
corresponding added transition is eventually generated; Ar IV is currently
the terminal production stage. Existing element products were copied with
their original timestamps and SHA-256 hashes to
`validation_2026-10-10/physical_thresholds/products/previous_atomic/` before
replacement. Older nebular results that depend on these energies are stale.
No nebular calculation has been rerun.

## Conservative terminal-stage photon check

`main.f90` invokes `species_terminal_photoionization_check(par%eion_max)`
after `species_setup` and before `ion_setup`. It rejects a transported band
extending above the physical ionization potential of any active element's
terminal stage. Zero-abundance elements do not constrain the band. Fully
stripped carbon has no further photoionization edge.

With the current roster, a 30 eV upper edge passes this particular check,
whereas 100 eV is rejected at C IV, whose next ionization would create the
unrepresented C V stage. Passing at 30 eV is not a claim that the full
chemistry, line cooling, or collisional-ionization truncation is complete.
The check is deliberately conservative: it tests the transported interval,
without deciding whether a particular source has a negligible spectral tail.
Existing 100 eV production inputs will therefore stop until physically
complete stage data are supplied. Lowering a photon-band limit changes the
physical input and must not be used to disguise missing ions.

The approved target has the following terminal photoionization edges:

| Element | Target terminal stage | Next ionization potential (eV) |
|---|---:|---:|
| C | V | 392.090518 |
| N | VI | 552.067329 |
| O | V | 113.898952 |
| Ne | V | 126.246910 |
| S | VII | 280.954393 |
| Ar | VII | 124.411944 |
| Mg | IV | 109.265414 |
| Fe | VII | 124.967145 |
| Si | V | 166.767426 |
| Cl | VII | 114.201349 |
| Ca | VI | 108.783736 |

All exceed 100 eV. The physical-energy criterion still selects 64 stages and
53 transitions. Collisional ionization and nonadjacent Auger products require
additional physical checks; this photon test does not implement them.

## Existing charge-transfer temperature correction

KF96 Table 1's thirteen existing higher-boundary rates had been evaluated
with one common 5000--50000 K clamp. Their published intervals differ.
They now use the existing form-8 log-log table calculation with individual
intervals and original published coefficients. Seven element files changed
in this separate correction; other numeric atomic records were checked for
equality. The private recombination-table capacity is 1024; the maximum
generated table has 653 nodes. Heating tables retain their separate capacity
of 40. Nonfinite, nonpositive, and nonincreasing rate tables are rejected.

N + H+ and N+ + H use their published 100--50000 K base-fit interval;
the endothermic activation factor uses the actual temperature. At 300 K,
the corrected N+ recombination fit is 0.652411 times the former clamped
evaluation. The Mg forward interval remains 5000--50000 K. Mg's inferred
reverse and the Fe first-boundary rates still require independent source
validity checks. Holding an endpoint rate outside a published fit interval
is an explicit continuation approximation, with no asserted accuracy there.

See [the charge-transfer audit](charge_exchange_expansion.md) for the
published rows, the Cloudy Si2 coefficient discrepancy, and unresolved
thermal-energy provenance. KF96 Table 1 has no heat or product-level column.
Existing comments attributing Cloudy's stored heat to that table were
corrected; this attribution correction does not establish those product
states. The nine verified added-boundary rate candidates do not authorize
their heating records or activation.

## Measured validation and limits

- The complete native debug executable builds with bounds and floating-point
  checks. The generated physical-energy module is compiled and linked into
  the species calculation.
- All 77 compiled potential entries were checked: 76 physical values and
  fully stripped carbon. Against a separately parsed source converted using
  exact SI constants, the largest relative energy difference is 2.22045e-16.
- All 32 loaded transitions have zero subthreshold absorption, positive
  absorption immediately above threshold, consistent photoelectron energy,
  and agreement between direct and ordered packet cross sections.
- A 100 eV current-roster run is rejected with the terminal-stage diagnostic;
  all eleven target terminal potentials exceed 100 eV.
- Compiled charge-transfer rates were checked at 26,049 samples. The thirteen
  table fits differ from the original published expressions by at most
  5.62851e-6 relative; nitrogen's first-boundary calculation differs by at
  most 1.43219e-14. Three malformed rate tables are rejected.
- The 25 focused Python tests pass. Adjacent equilibrium, restart mapping,
  charge conservation, cell convergence, detailed balance, H/He exchange,
  and charge-transfer heating tests pass after production changes. Measured
  worst ionizing/recombining reaction-conservation differences are
  3.565e-16 and 3.169e-16, respectively.

Results are in `validation_2026-10-10/physical_thresholds/`,
`charge_transfer_ranges/`, and `final_photon_energy/`. Finite controllers
record commands, source hashes, child PIDs, return codes, and JSON acceptance
checks. Two failed validation-script attempts are preserved with recovery
notes; both stopped before their associated production updates.

No expanded production model, matched comparison with another code, full
regression suite, or temperature/density/composition profile comparison was
run. These checks address the changed energy and rate calculations and their
consumers. Missing charge-transfer product energies, Fe V/VI low-temperature
DR, N IV energy consistency, cooling convergence, and Auger physics still
prevent acceptance of the complete expanded physical model.
