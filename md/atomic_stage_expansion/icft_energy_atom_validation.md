# O V and Cl VI excitation energy and atom validation

Date: October 10, 2026. Production remains at 43 metal stages; no expanded
production data or transported nebular calculation was activated.

## Physical energy correction

Original CHIANTI mixed observed/theoretical energies reverse five O V and
one Cl VI radiative transitions. O V also contains rounded zero collision
gaps. Indices identify levels and do not define physical energy ordering.
External files were preserved; no A value, collision strength or transition
identity was reversed, discarded or replaced to hide an inconsistent gap.

`atomic_excitation_energies.py` now owns the common physical energy selection
for cooling, detailed balance and n-level generation. The existing checked
N IV selection remains in its source-specific module. The obsolete general
selection function in that module was removed, and callers use the shared
physical energy module.

For O V, levels 1--40 retain observed energies, exactly matching the Stout
NIST subset. All remaining 126 retained CHIANTI levels are mapped to the
original ICFT target energies in the high-precision Stout table. Each map is
unique within the CHIANTI rounding interval, and checks statistical weight,
orbital occupation and term identity; parent-term labels are normalized
without changing a level identity. The map is one-to-one. At CHIANTI level
87 the correct Stout index is 89, not 87. A direct index copy is wrong.
The source hash, complete level and transition counts and every directed
physical energy gap are checked before a subset is selected.

For Cl VI, the first 44 observed CHIANTI levels are retained; above them,
the original CHIANTI ICFT target energies are used consistently. The Stout
header identifies this contiguous low subset as NIST, plus separate high
observations that can reverse mixed source ordering. The selected splice
preserves the low-level observations and avoids inconsistent high-level
splices. It is an explicit approximate physical model, not a new measured
energy dataset. All 171 levels, 3,760 radiative rows and 14,535 collision
pairs pass physical energy-direction checks.

Both models use the final 2014 Be-like/Mg-like calculations documented in
[reference provenance](reference_sources.md). Their scattering resonance
positions are not recalculated using the observed energy splice. Low-T
resonance accuracy and convergence beyond the supplied target expansion
remain unquantified. The 1,000--100,000 K fit interval is numerical and does
not certify atomic accuracy. No new state-resolved ionization, photoexcitation,
radiation trapping or finite-density suppression of recombination is supplied.

## Measured results

| Ion | Maximum cooling-fit error (%) | Selected levels | Radiative rows | Collision rows | Maximum native/full-atom difference (%) |
|---|---:|---:|---:|---:|---:|
| OV | 0.2815946 | 29 | 133 | 406 | 0.1512383 |
| ClVI | 0.7186875 | 29 | 130 | 406 | 0.1191395 |

Each rounded cooling fit passes the unchanged 1% criterion on 801 temperatures.
The actual debug Fortran readers/evaluators pass. Both 29-level models fit
within native MAXLEV=40, MAXTR=600 and 16 Chebyshev coefficients. Every
collision fit meets the existing 0.2% target at its acceptance nodes.
No native capacity or public data layout changed.

An independently assembled full statistical-equilibrium matrix uses all
166 O V or 171 Cl VI levels and source spline collision strengths. Complete
801-temperature profiles were compared with the native 29-state atom at
n_e=1e-6, 100 and 1e6 cm^-3 (2,403 samples for each ion). The measured maximum
differences remain below 0.5%. Native populations are finite, nonnegative
and normalized. A separate scalar full-matrix assembly also checks nine
states for each atom; representative samples alone were not used as the
full convergence acceptance.

| Input or physical condition | Full-source calculation | Native selected atom | Mismatch |
|---|---|---|---|
| Ion and energy vector | Same selected O V / Cl VI vector | Same, rounded energy serialization | Numerical rounding only |
| Temperature | 801 logarithmic samples, 1e3--1e5 K | Same grid | None |
| Electron density | 1e-6, 100, 1e6 cm^-3 | Same | None |
| Electron distribution | Maxwellian | Maxwellian | None |
| Collision strengths | Source splines | Chebyshev fits of those splines | Controlled numerical fit |
| Bound-level count | 166 / 171 | 29 | Intentional truncation under test |
| Spontaneous rates | Source A values | Same selected A values | Omitted high levels |
| Radiative transfer | Optically thin escape; no photoexcitation | Same | None |
| Geometry, gravity, velocity | No spatial dynamics | Same | Not applicable |
| Ionization, charge transfer | Excluded from level-population comparison | Same | Not a nebular ionization comparison |

Complete source, fit and three-density curves were visually inspected. No
clipped labels, overlapping plot elements or unexplained discrepancies were
found. The earlier twelve-ion cooling figures were also inspected and their
[workflow record](remaining_cooling_workflow.md) was updated.

Evidence and vector figures:

- [Cooling controller](validation_2026-10-10/icft_energy_cooling/controller/state.json)
- [Atom controller](validation_2026-10-10/icft_energy_cooling/atom_controller/state.json)
- [Density profile controller](validation_2026-10-10/icft_energy_cooling/density_controller/state.json)
- [Complete density comparisons](validation_2026-10-10/icft_energy_cooling/density_profiles_result.json)
- [O V cooling](validation_2026-10-10/icft_energy_cooling/products/OV/OV.pdf) and [density curves](validation_2026-10-10/icft_energy_cooling/products/OV/OV_density_profiles.pdf)
- [Cl VI cooling](validation_2026-10-10/icft_energy_cooling/products/ClVI/ClVI.pdf) and [density curves](validation_2026-10-10/icft_energy_cooling/products/ClVI/ClVI_density_profiles.pdf)

Three focused ICFT source/independent-balance tests and six existing stationary
cooling tests pass. Two Cl VII source-domain tests pass separately. No full
regression or nebular transport comparison was run: these changes affect
excitation energy selection and staged coefficient generation, and production
has neither added O V nor added Cl VI.
