# Excitation convergence across the cooling-table temperature grid

## Physical assessment

The earlier reduced Ar V, Mg IV, Fe VII, and Ne V atoms satisfy their measured
1,000--100,000 K convergence checks. They do not satisfy the same 0.3% target
throughout the numerical cooling grid, 100--1,000,000 K. An accepted finite
interval cannot establish convergence outside it. Full collision-connected
bound source atoms have therefore passed compiled source validation; the
additional states change internal atomic storage, not the approved 64-ion
roster or the transport photon band.

The values below were read directly from the completed extended diagnostic
JSON files. Each maximum covers 401 logarithmic temperatures and four electron
densities. The density-ratio difference compares each atom's cooling at the
specified density divided by that atom's cooling at 1e-6 cm^-3. It is distinct
from the absolute cooling difference.

| Ion | Reduced / full levels | Maximum cooling difference (%) | Maximum density-ratio difference (%) |
| --- | ---: | ---: | ---: |
| ArV | 110 / 197 | 4.373428 | 1.419933 |
| MgIV | 110 / 157 | 17.926576 | 2.530718 |
| FeVII | 170 / 189 | 5.327865 | 4.799698 |
| NeV | 25 / 304 | 51.292409 | 1.752487 |

## Matched input conditions

| Condition | Full source | Reduced source |
| --- | --- | --- |
| Initial ion and raw source | Same identified CHIANTI ion, or coherent Stout Fe VII source | Same raw source and physical energy selection |
| Temperature | 401 logarithmic samples, 100--1,000,000 K | Identical samples |
| Electron density | 1e-6, 1e2, 1e6, 1e10 cm^-3 | Identical samples |
| Electron distribution | Maxwellian; detailed balance at the actual temperature | Same |
| Collision interpolation | Source cubic/PCHIP or log-temperature Stout PCHIP | Same original retained rows |
| Low-temperature Stout treatment | Adopted constant first-node Upsilon | Same |
| Radiative processes | Spontaneous emission, optically thin | Retained spontaneous transitions, same A values |
| Excluded processes | Photoexcitation, photon trapping, autoionization, excited-state ionization, other colliders | Same |
| Spatial domain, gravity, velocity, spectrum, boundaries | No spatial or transport calculation | Same |
| Intentional mismatch | All collision-connected bound states | First selected states and their transitions only |

The main session inspected all complete cooling profiles and density-ratio
profiles for all four densities. Ne V's high-temperature cooling deficit is
visible in the complete curves; high-density ratio errors also become visible
outside the earlier interval for all four ions. These are atomic statistical-
equilibrium diagnostics and do not establish a nebular temperature,
composition, density, or velocity profile.

- ArV: [cooling profiles](validation_2026-10-10/excitation_temperature_range/ArV/cooling_profiles.pdf), [density-ratio differences](validation_2026-10-10/excitation_temperature_range/ArV/density_ratio_differences.pdf), [measured values and source hashes](validation_2026-10-10/excitation_temperature_range/ArV/result.json).
- MgIV: [cooling profiles](validation_2026-10-10/excitation_temperature_range/MgIV/cooling_profiles.pdf), [density-ratio differences](validation_2026-10-10/excitation_temperature_range/MgIV/density_ratio_differences.pdf), [measured values and source hashes](validation_2026-10-10/excitation_temperature_range/MgIV/result.json).
- FeVII: [cooling profiles](validation_2026-10-10/excitation_temperature_range/FeVII/cooling_profiles.pdf), [density-ratio differences](validation_2026-10-10/excitation_temperature_range/FeVII/density_ratio_differences.pdf), [measured values and source hashes](validation_2026-10-10/excitation_temperature_range/FeVII/result.json).
- NeV: [cooling profiles](validation_2026-10-10/excitation_temperature_range/NeV/cooling_profiles.pdf), [density-ratio differences](validation_2026-10-10/excitation_temperature_range/NeV/density_ratio_differences.pdf), [measured values and source hashes](validation_2026-10-10/excitation_temperature_range/NeV/result.json).

## Native validation and limits

The `atomic-full-bound-gnu` controller completed all eight declared source
and registered-table checks using isolated GNU Fortran O2 objects. All four
Intel O2 source/profile checks also passed, along with the Ar V, Mg IV, and
Fe VII tables. The remaining duplicate Intel Ne V table was deliberately
interrupted after its GNU counterpart passed; its controller records an
interruption rather than a successful step. Partial products are preserved
and no automatic restart is scheduled. Neither full-bound controller is running.
See [deliberate duplicate interruption](validation_2026-10-10/full_bound_native/deliberate_duplicate_stop.json).

The matched 401-temperature/four-density Intel/GNU cooling differences are
3.33e-16 (Ar V and Mg IV), 3.30e-13 (Fe VII), and 2.22e-16 (Ne V).
The complete profiles and population differences are recorded in
[compiler comparison](validation_2026-10-10/full_bound_native/compiler_profile_comparison.json).
Maximum GNU table application errors were 0.233029%, 0.145363%, 0.115808%,
and 0.135671%, respectively. These are numerical checks against the identified
source model, not independent certification of its physical accuracy.
No dependencies were installed and no production compiler setting changed.
A matched benchmark of 100 solves for the full 197-state Ar V atom at
10,000 K and 1e4 cm^-3 gave exactly the same printed cooling coefficient,
7.340666730132349e-20 erg cm^3 s^-1. CPU times were 1.315972 s (Intel) and
0.250334 s (GNU). This is one fixed-state timing, not the complete registry
or a production transport runtime. Build flags and hashes are preserved in
`validation_2026-10-10/full_bound_native/`.

Source completeness means completeness within the identified bound excitation
model. It does not certify low-temperature scattering resonance positions,
missing forbidden radiative data, neglected autoionization, an ionization
closure at 1,000,000 K, or extrapolation accuracy of the separate cooling fits.
The fit interval remains 1,000--100,000 K. Production activation and the
requested matched nebular comparisons remain separate unfinished steps.
