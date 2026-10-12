# Positive subnormal excitation cooling

## Confirmed defect and physical judgment

A positive representable cooling coefficient must remain positive unless
arithmetic underflow makes it unrepresentable. The previous registered
cooling path divided by `max(reference, tiny())` during table storage and
again during application. `tiny()` is the smallest normal number, not the
smallest positive representable number. This changed physically positive
subnormal coefficients and did not merely regularize exact zeros.

An actual compiled two-level diagnostic at 100 K and electron density
1e-6 cm^-3 returned direct cooling 1.54961155261344e-310 erg cm^3 s^-1.
The previous registered value was 1.07919831731e-312, and application returned
7.515877166e-315. These are measured results in
[original diagnostic](validation_2026-10-10/cooling_subnormal_diagnostic/result.json).
The normal floating-point boundary is 2.2250738585072014e-308.

## Installed correction

The installed module stores absolute Lambda at each registered temperature
and electron density. Positive endpoint interpolation uses logarithms;
exact zero endpoints retain linear interpolation. A source-converged atom
uses `fit * (Lambda / reference)` when the reference is positive, and the
absolute computed Lambda when the reference is exactly zero. For a partial
atom, the existing approximation becomes `Lambda + (fit - reference)` when
`fit >= reference`, and `(fit/reference)*Lambda` otherwise. No positive
reference is replaced with a normal-number floor. Public routine signatures
and the original 401 by 161 grid were unchanged by this correction.
The later source-limited grid extension preserves those rows for unlimited
BT/legacy atoms while respecting shorter physical-temperature tables; see
[its validation](validation_2026-10-10/source_temperature_registry/README.md). The registry's private storage changes
from ratios and a reference column to absolute cooling.

The correction now replaces `src/nlevel_cooling_mod.f90`. All earlier
controllers completed or were deliberately stopped before shared objects were
rebuilt. Their historical results remain attributable to archived source and
linked-object hashes. The actual `make DEBUG=1 -j4` build passed. The rebuilt
production objects passed all four native fixtures (48 rows each), linked
explicitly with `-no-ftz`; default production compiler settings are unchanged.
See [installation](validation_2026-10-10/cooling_subnormal_installed/installation.json).

## Scoped acceptance

The isolated candidate preserved 24 positive subnormal samples in each of
the complete and partial atom fixtures. Four fixtures also cover exact zeros,
temperature/density interpolation, and both partial-atom coverage branches.
The maximum sampled table/direct difference is 3.313589e-5. At the original
100 K node the registered value retains the direct 1.54961155e-310 result.
Actual S I, Fe IV, and Ca IV registration and public application checks each
built all 64,561 cells. The maximum checked application differences are
0.168139%, 0.152081%, and 0.163713%, respectively. Source hashes, driver
hashes, fixture measurements, and acceptance results are preserved in
[candidate result](validation_2026-10-10/cooling_subnormal_candidate/result.json).

These checks exercise modified storage, interpolation, and the public
cooling combination. The excitation rate equations and population solver
are unchanged. No spatial nebular model, ionization-stage activation, or
complete cross-code comparison was run for this numerical correction.

## Separate existing partial-atom limitation

Correcting the normal-number floor does not resolve the existing partial-atom
residual approximation. A corrected nine-ion public cooling check passed
table interpolation and the low-density limit, but failed the combined-cooling
criterion for seven ions. The greatest measured relative discrepancy is
27.7687% for C II. After collisional de-excitation, small errors in the
interpolated low-density reference can dominate the retained `fit - reference`
residual. These are relative differences in the combined ion coefficient,
not a measured change in a nebular cooling budget. Source-converged density
dependence is being investigated; the failure is not an acceptance result.
See [actual failed diagnostics](validation_2026-10-10/cooling_ne_corrected_diagnostics/result.json).
