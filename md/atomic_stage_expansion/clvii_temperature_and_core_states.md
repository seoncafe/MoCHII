# Cl VII temperature domain and core-excited states

Final outer-shell and inner-shell articles were downloaded to
`~/MoCHII/references/`; their methods and final journal identities were read.
See [provenance](reference_sources.md). Their supplied temperature interval
is 200(q+1)^2--2e6(q+1)^2 K with q=Z-11. For chlorine q=6, this gives
9,800--98,000,000 K, matching the actual Stout temperature grid. No
1,000--9,800 K collision-strength accuracy is established by those papers.

The inspected Stout atom has 161 levels and 9,397 collision/radiative rows.
Its 32 outer-shell levels contribute 496 collision and 496 radiative rows.
The first core-excited level is at 204.55780946898412 eV; the physical
Cl VII ionization potential is 114.20134924002703 eV. All 129 core-excited
levels lie above that threshold. They cannot be treated as an ordinary
bound-level cooling atom without autoionization and competing radiative
branching. The source configuration set matches a 32+134-5 combination;
its header cites only the outer-shell paper. Numerical provenance for every
core-excited row remains unverified. External source files are preserved.

The user adopted a constant-Upsilon continuation below 9,800 K. The point
of use in `stout_excitation_cooling.py` states that this is an unverified
low-temperature resonance approximation; the numerical use interval is
1,000--100,000 K. For T below 9,800 K only Upsilon is fixed at its lowest
source value. The Maxwellian exponential, the T^-1/2 prefactor, and detailed
balance retain the actual temperature. Cooling and excitation rates are not
held constant. Above 9,800 K, positive strengths use PCHIP in log T/log
Upsilon; temperatures above the final source node are rejected.

The shared excitation reader selects the 32 bound Stout levels explicitly;
it does not create substitute CHIANTI files. All 496 retained radiative and
496 collision rows preserve their original level identities. Type 2 supplies
only an invertible coordinate for native Chebyshev storage; the source
interpolation remains log-PCHIP. Native file layout and capacity are unchanged.
The approximation applies to all retained excitation/de-excitation pairs.

Four focused tests verify the physical bound-state scope, invalid-temperature
rejection, endpoint continuation with actual-temperature Maxwellian cooling,
and preservation of all published collision-strength nodes. Six existing
stationary low-density balance tests also pass. Published-temperature samples
are preserved at [table](validation_2026-10-10/clvii_published_temperature_cooling.txt).
New numerical products and checks are staged separately under
`validation_2026-10-10/clvii_constant_strength/`; production activation still
requires the remaining physical and convergence checks.

The finding concerns final-ion bound excitation, not a measured Auger yield.
Core-state autoionization and photoionization Auger branching belong to the
separate ionization/energy audit; one must not infer those rates from this
bound cooling table.

## Measured validation of the adopted continuation

The staged low-density fit and native evaluation passed over 801 logarithmic
temperatures from 1,000 to 100,000 K. The 32-level native statistical balance
was checked at every temperature for all three densities, 2,403 states total.
All populations are finite, nonnegative, and normalized within 2e-13. Native
cooling is positive throughout the tested interval. These results establish
numerical representation of the adopted source model, not cold resonance accuracy.

| Quantity | Directly measured maximum |
| --- | --- |
| Rounded low-density cooling fit / adopted source | 0.5728601711129544% |
| Native cooling-fit evaluation / Python coefficient evaluation | 4.440892098500626e-16 relative |
| Serialized collision strength / source, all 496 pairs and 1,601 temperatures | 0.1999403980527781% |
| Serialized strength variation at 1,000, 2,000, 5,000, and 9,799 K | 0 relative |
| Native finite-density cooling / independent full 32-level source | 0.1348220130716005% |
| Native population / full source population | 1.2323475573339238e-14 absolute |

The cold strength check evaluates the serialized data with the same clipping
and type-2 descaling equations inspected in `src/nlevel_mod.f90:ups_eval`.
The independent native matrix check exercises the linked Fortran routines
`nlevel_load`, `nlevel_cooling_prepT`, `nlevel_cooling_solve`, and
`nlevel_populations`. A constant strength does not remove the physical
Boltzmann suppression; neither rate floors nor cooling floors were introduced.

| Input condition | Independent source | Native calculation | Mismatch |
| --- | --- | --- | --- |
| Levels | All 32 outer bound states | Same 32 states | None |
| Energies and A values | Original Stout | Serialized Stout | Numerical decimal rounding |
| Upsilon | Log-PCHIP source table, cold endpoint held | Type-2 Chebyshev, lower coordinate held | Bounded fit and serialization error |
| Temperatures | Same 801 values, 1e3--1e5 K | Same grid | None |
| Electron densities | 1e-6, 100, 1e6 cm^-3 | Same | None |
| Electron distribution | Maxwellian | Maxwellian | None |
| Processes | Electron excitation/de-excitation and spontaneous emission | Same | None |
| Radiation trapping and photoexcitation | Excluded | Excluded | None |
| Core-state autoionization and excited-state ionization | Excluded | Excluded | None |
| Geometry, gravity, and boundaries | Fixed-state calculation; inapplicable | Same | None |

The entire cooling and excited-level population profiles were visually
inspected from these vector PDFs; the source and native curves agree within
the measured bounds:

- [Low-density fit and residuals](validation_2026-10-10/clvii_constant_strength/products/ClVII/ClVII.pdf).
- [Three-density cooling and residuals](validation_2026-10-10/clvii_constant_strength/products/ClVII/ClVII_density_profiles.pdf).
- [All excited-level population profiles](validation_2026-10-10/clvii_constant_strength/products/ClVII/ClVII_population_profiles.pdf).

Evidence: [cooling audit](validation_2026-10-10/clvii_constant_strength/products/ClVII/ClVII.json),
[native cooling evaluation](validation_2026-10-10/clvii_constant_strength/native_cooling/native_result.json),
[atom balance](validation_2026-10-10/clvii_constant_strength/atom_result.json), and
[serialized strengths](validation_2026-10-10/clvii_constant_strength/serialized_strength_result.json).
The final `atom_native_controller/state.json` records completion. No Cl VII
controller remains running. The earlier failed controller records and the
initial global-fit atom are preserved and are superseded for acceptance.

## Corrections encountered during validation

The initial single-polynomial fit included the cold join and reached 1.33%
strength error; 188 pairs exceeded 0.2%. The accepted representation instead
uses the original 9,800 K table endpoint as its lower coordinate and the
existing native coordinate clipping. All pairs now pass without increasing
coefficient capacity or changing the native format.

Two defects in the new validation utility were corrected before acceptance:
a NumPy Boolean was converted to a Python Boolean before JSON serialization,
and the native executable output directory is created before compilation.
Their failed execution records were preserved; the corrected linked native
calculation and complete numerical checks passed afterward.

Twenty focused tests passed across the cold Stout source, stationary balance,
ICFT energy consistency, excitation underflow, selective generation, and
N IV collision fitting. No full regression or nebular calculation was run:
the affected path is staged excitation-data generation and native statistical
balance. Production remains at 43 metal stages; the Cl VII atom and cooling
products remain staged. Remaining atoms' finite-density convergence and the
terminal-ionization/Auger accounting must be completed before activation.
