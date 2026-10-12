# Density-dependent cooling from complete or converged excitation atoms

## Physical correction

A complete atom has no omitted excitation cooling. The difference between a
fitted low-density cooling coefficient and an interpolated reference cannot
be interpreted as a physical contribution from missing levels. The previous
combination used `fcov = min(1, lamref / lam_fit)` even for a complete atom,
leaving `lam_fit * (1 - fcov)` unsuppressed at arbitrarily high density.
This turns a fit or interpolation error into spurious emission.

The linked S I calculation established this defect directly. At 101,158 K
and 8.91e9 cm^-3, the reference-column interpolation error was approximately
2.98e-5 relatively, but the applied cooling differed by 1.975%. The raw
registered cooling table itself remained within its 0.5% interpolation target.
The correction is required because the residual has no physical emitting
levels, independently of its numerical size.

For a complete or independently converged atom, the corrected expression is

```text
Lambda(T, ne) = Lambda_fit(T) * Lambda_atom(T, ne) / Lambda_atom(T, ne_floor)
```

The explicit optional version-2 record `DENSITY_RATIO 1` selects this model.
The generator writes it for the full bound collision-connected source. A
reduced atom requires full-source convergence before receiving that record;
format version alone does not establish convergence. The atom, source hashes,
full profiles, and acceptance record remain traceable. Existing files without
the record retain their partial-atom approximation; no user parameter changed.

S III's earlier six-level atom is partial: the actual bound source contains
53 levels. It was therefore replaced with the full source, rather than being
marked complete merely because a new representation was available. Its
53 levels, 355 radiative rows, and 1,378 collision rows passed native source,
population, cooling, and emissivity checks at 801 temperatures and four
densities. Its maximum relative native cooling difference is 1.60e-13.

## Inputs and measured acceptance

| Condition | Native atom and independent source |
|---|---|
| Temperature | 801 logarithmic samples, 1,000--100,000 K |
| Electron density | 1e-6, 100, 1e6, 1e10 cm^-3 |
| Radiation and spatial boundaries | None; local statistical equilibrium only |
| Processes | Electron excitation/de-excitation and spontaneous emission |
| Energies | Identical checked bound-level energies and Maxwellian gaps |
| Collision interpolation | Original source nodes and adopted cubic/PCHIP rule |
| Excluded processes | Photoexcitation, line trapping, excited-state ionization |

The measured maximum density-ratio differences of selected versus full source
are 0.114098% for Ar V, 0.0731532% for Mg IV, 0.188543% for Fe VII, and
0.0842984% for Ne V. These establish the marker only over the stated
temperature/density comparison. Full-source curves and all four density
profiles were inspected; they do not certify a larger physical domain.
The independent 100--1,000,000 K source-range controller tests that extension
separately. Its completion status and model-acceptance status are distinct.

The unchanged numerical table contains 401 temperatures and 161 densities.
S I's actual table passed with maximum interior raw-table error 0.163364%
and maximum corrected-application error 0.168139%. S III's full source passed
with 0.0951391% and 0.0949361%, respectively. Fe IV's 37-level full source
passed with 0.155935% and 0.152081%. The calculation and resource records are
under `validation_2026-10-10/registered_source_nodes/density_ratio/`.

Fe IV's native calculation completed successfully. Its first Python validation
failed to parse a valid Fortran real printed without `E` at a three-digit
exponent. The validator now reads that notation and checks the complete saved
output and native exit status. The 64,561-cell table was not recalculated
merely to change printed notation; the failed controller is preserved.

## Installation and remaining work

The accepted S I cooling fit and five-level atom, and S III's full 53-level
atom, are installed. Original S III data are preserved. The aggregate cooling
table now has 29 ions and 152 terms, retaining primary coefficient precision.
The installed S I file was opened and evaluated by the actual `cooling_load`
and `cooling_eval` routines at 801 temperatures; its maximum difference from
the rounded coefficients is 2.818e-12. See
`validation_2026-10-10/installed_source_density_atoms/installation.json`.

Ne V and the larger Fe/Ar/Mg registered tables remain under validation in
`atomic-registered-remaining`. Missing cooling sources, terminal electron
ionization, nonadjacent Auger rates and energy accounting, expanded-roster
activation, and matched nebular comparisons remain separate incomplete work.
No nebular model was rerun. Previous density-dependent S III model products
are stale for the new full source and must be recalculated before reuse.
