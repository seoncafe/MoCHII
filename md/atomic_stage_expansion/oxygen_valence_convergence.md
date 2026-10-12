# Oxygen bound-level radiative consistency

The currently selected O II and O III hybrid energies cannot be paired with
the unmodified full CHIANTI radiative arrays. In the full bound atoms, 31 O II
and 48 O III positive radiative rates have the opposite energy order. The
installed O III six-level atom has already received the absolute Storey and
Zeippen (2000) ground optical rates; this study does not change that atom.

The machine-readable diagnostics are in
`validation_2026-10-10/oxygen_valence_convergence/`. The input source is local
CHIANTI v11.0.2 `.elvlc`, `.wgfa`, and `.scups`; the selected energies are
calculated by `tools/fitting/atomic_excitation_energies.py`. The script
`assess_clean_prefixes.py` and its `ground_excitation_prefixes.npz` evaluate
all 801 temperatures from 1,000 to 100,000 K. The PDF displays the entire
temperature profile and was rendered and visually checked. These values are
the ground-excitation coronal power, not a finite-density solution. They do
not depend on uncertain high-level radiative branches when all excited levels
can eventually radiate to the ground. The latter condition holds in both
complete source graphs after excluding the inverted directions.

| Atom | Lowest inverted positive A row | Small clean prefix | Maximum omitted ground-excitation power against all bound levels |
| --- | --- | ---: | ---: |
| O II | 2 to 3, source 1.32e-7 s^-1, selected gap -19.8 cm^-1 | No prefix above two levels | 82.169% for 5 levels; 7.884% for 64 levels |
| O III | 80 to 90, source 2.16e4 s^-1, selected gap -2634.7 cm^-1 | 89 levels | 1.001% for 89 levels; 2.038% for 64 levels |

Thus a simple contiguous prefix does not meet the approved 0.3% coronal
omission criterion. Its low-temperature agreement alone is insufficient.

## Transition identity and a coherent energy option

[CHIANTI Technical Report 10](https://www.chiantidatabase.org/tech_reports/10_wgfa/chianti_report_10.pdf)
defines the `.wgfa` wavelength, weighted oscillator strength, and spontaneous
probability. Its negative wavelength indicates a theoretical transition
energy when an observed energy is unavailable. The Einstein identity for an
electric-dipole oscillator strength is
`A_ki = 6.6703e15 gf/(g_k lambda_A**2)`; see the final printed
[Wiese and Fuhr (2007)](https://doi.org/10.1063/1.2740642) compilation
archived in `~/MoCHII/references/Wiese_Fuhr_2007_JPCRD_36_1287.pdf`.

All 48 inverted O III rows and 30 of 31 inverted O II rows have nonzero `gf`
and connect opposite-parity configurations. The frequency inferred from each
source A/gf agrees with its theoretical energy gap to within 0.315% for
O III and 0.171% for O II. This supports preserving their electric-dipole
line strengths when observed and theoretical energy ordering differs:

`A_physical = A_source * (abs(deltaE_physical)/deltaE_source)**3
              * g_source_upper/g_physical_upper`.

Here `deltaE_source` must be the actual frequency of the source A and gf,
which can be reconstructed from their Einstein identity. Across all positive
`gf` rows, each inferred frequency matches either the theoretical gap or
the recorded wavelength to 2%; no row is unidentified by that test. However,
nonzero `gf` by itself does **not** establish electric-dipole character:
the O III ground 2-to-3 forbidden line has `gf=1.947e-8` and
`A=9.632e-5 s^-1`. Opposite parity and the actual level configurations must
also be checked. The isolated `audit_e1_reciprocity.py` records every source
row, frequency match, reversal, and conditional electric-dipole correction in
`e1_reciprocity.json`. It does not write a candidate atom.

The 3,978 O II and 2,949 O III positive-A opposite-parity rows all have
nonzero gf and satisfy the E1 angular momentum selection condition
`abs(delta J) <= 1`, excluding `J=0` to `J=0`. Their inferred A/gf
frequencies agree with either source theoretical gaps or recorded wavelengths
to within 0.779% for O II and 0.387% for O III. All inverted opposite-parity
rows agree with theoretical gaps within 0.171% and 0.315%, respectively.
These are source-identity checks, not independent accuracy estimates for
the atomic wavefunctions.

The remaining inverted O II row is the `2D` fine-structure transition 2 to 3,
which has `gf=0`. The Cloudy c25 Stout O II source lists the lower level as
its level 2 (`2D5/2`, 26810.550 cm^-1) and the emitting upper level as its
level 3 (`2D3/2`, 26830.570 cm^-1). Stout gives
`1.30e-7 s^-1 M1` and `1.49e-20 s^-1 E2`. Mapping those identities to the
CHIANTI selected levels would reverse the CHIANTI index order and require
the distinct third- and fifth-power frequency factors for the observed
19.8 cm^-1 gap if the Stout probabilities are adopted. The CHIANTI source
itself records a 1.32e-7 s^-1 rate with a positive wavelength corresponding
to that observed 19.8 cm^-1 gap; its direction must be established from the
observed energies, as the independent Stout level identities confirm. This
is the only inverted row without a supported dipole-strength correction.

A simpler source-coherent high-level option retains the observed first six
O III levels and theoretical energies above level 6. Its 3,164 radiative
rows and all collision pairs have positive source order. O II observed first
five plus theoretical higher levels has only the identified 2-to-3
fine-structure inversion; two low-level collision pairs also reverse and
must use the physical energy order. In this option, all high-to-high
transitions retain their original theoretical energy gaps. Low-to-high
transitions still mix observed low energies with theoretical high energies:
among same-parity positive A rows, 59 of 132 O II rows and 5 of 93 O III rows
have a gap difference above 1%, with maxima 2.19% and 2.65%, respectively.
The source does not identify every M1 and E2 component, so these rates cannot
yet be corrected by a universal exponent. This is an explicit approximation,
not an exact radiative reconstruction.

An implementable observed-energy model would retain the current hybrid
energies, reconstruct every opposite-parity electric-dipole A from its source
line strength and the selected physical photon frequency, and reorient the
48 O III and 30 O II inverted dipole lines. It would orient the O II
fine-structure M1 branch using the identified observed 2D levels and the
original CHIANTI observed-gap A=1.32e-7 s^-1, with the Stout multipole
components as an independent direction and scale check. Every remaining
same-parity A would remain its published atomic estimate. This last choice
is a declared approximation to an unspecified M1/E2 mixture, not a claim that
its original wavelength always equals the source calculation frequency.
It avoids discarding observed high-level energies, while requiring a direct
finite-density sensitivity calculation before acceptance.

`audit_multipole_frequency.py` records a conditional frequency-sensitivity
envelope in `multipole_frequency.json`: for each positive same-parity rate,
take `r = abs(selected gap) / (1e8/abs(wgfa wavelength))`. If the rate were
pure M1 calculated at that wavelength, preserving its matrix element would
multiply A by `r**3`; pure E2 would use `r**5`. The actual source A may
already use a different calculation frequency, and configuration mixing or
other multipoles can change the matrix element. Therefore these factors
bound only the stated M1/E2 frequency model, **not** the uncertainty of the
published atomic rate.

| Selected energies | O II same-parity rows; M1/E2 conditional maximum change | O III same-parity rows; M1/E2 conditional maximum change |
| --- | ---: | ---: |
| Current hybrid | 188 rows; 17.52% / 30.87% | 215 rows; 28.40% / 51.69% |
| Observed low, theoretical high | 188 rows; 20.00% / 35.50% | 215 rows; 18.75% / 33.16% |

The hybrid choice has 9 O II and 7 O III same-parity rows above 10% under
the conditional M1 factor. The source-order choice has 60 and 21. Under the
conditional E2 factor the respective counts are 9 and 11 for the hybrid,
and 128 and 33 for source order. These row counts show why replacing observed
high energies with theoretical values does not automatically yield a better
radiative model: some CHIANTI positive wavelengths already use observed
high-level gaps. They do not determine the effect on total finite-density
cooling without solving the level populations and branching.

`compare_energy_choices.py` evaluates the full bound ground-excitation power
for this source-order energy option against the current hybrid selection on
the same 801-temperature grid, without changing collision strengths or
radiative data. The maximum absolute cooling changes are 2.019% for O II and
3.520% for O III; at 10,000 K they are -0.000187% and -0.000234%, and at
100,000 K they are -1.836% and -2.782%, respectively. The entire curves are
in `energy_choice_ground_excitation.pdf` and the numerical arrays are in
`energy_choice_ground_excitation.npz`. The PDF was rendered and visually
checked. These changes measure the different selected excitation energies;
they do not establish which high-level energies are physically more accurate.

| Input condition | Current hybrid | Source-order option |
| --- | --- | --- |
| Bound levels | O II 196; O III 177 | Identical |
| Level energies | Observed wherever available; otherwise theoretical | Observed O II 1–5 or O III 1–6; theoretical higher levels |
| Effective collision strengths | Original CHIANTI `.scups` curves | Identical |
| Temperature and density limit | 801 points, 1,000–100,000 K; ground-population coronal limit | Identical |
| Excitation loss | Ground excitation rate times selected physical level gap | Same equation, different selected gap |
| Radiation | Not evaluated in this comparison | Not evaluated in this comparison |

The full oxygen source atoms were installed after the private full-bound
candidates under `validation_2026-10-10/oxygen_valence_convergence/candidate_atoms/`
passed the independent balance and source checks.
`build_full_candidates.py` creates one baseline per ion and conditional M1/E2
frequency alternatives. It chooses the source theoretical gap or `.wgfa`
wavelength nearest the frequency inferred from source A/gf for each
opposite-parity E1 row, then preserves its line strength at the selected
physical gap and the physical emitting statistical weight. The inferred
rounded frequency is evidence for the choice, not silently treated as an
exact source value. For the inverted E1 rows, the maximum difference between
the chosen and inferred frequencies is 0.171% (O II) and 0.315% (O III),
giving at most 0.513% and 0.951% in the corresponding cubic A factor from
that identification alone. The underlying source A and gf errors remain
unquantified. The O II identified M1 rate remains 1.32e-7 s^-1; the O III
ground optical rates are exactly 6.791e-3 and 2.046e-2 s^-1 in all candidates.

The native source-node file retains all original collision nodes and scaling
parameters. Its level energies set the collisional Boltzmann factors, photon
energies, and detailed balance. The baseline keeps every other same-parity
A as an approximate positive effective rate. The conditional variants multiply
only non-ground same-parity rates by `r**3` or `r**5`, with `r` based on the
source `.wgfa` wavelength. This is a sensitivity test, not a frequency-based
accuracy bound: a `.wgfa` wavelength need not be the frequency used in the
source atomic A calculation, especially when gf is zero. Where same-parity
gf is nonzero, the A/gf identity can supply an additional source-frequency
check, but it does not identify the M1/E2 mixture.
Eight same-parity rows in each ion have nonzero gf. All are within the low
five O II or low six O III levels. Their A/gf frequencies match the recorded
wavelength frequencies to within 0.031% for O II and 0.223% for O III;
several differ from the theoretical gaps by more than 10%. The private
M1/E2 sensitivity variants preserve all of these low-level rates exactly.

| Candidate inputs | O II | O III |
| --- | ---: | ---: |
| Bound levels | 196 | 177 |
| Radiative pairs | 4,166 | 3,164 |
| Collision rows | 19,110 | 15,576 |
| Collision source nodes | 286,690 | 233,850 |
| E1 rows transformed | 3,978 | 2,949 |
| Inverted E1 rows reoriented | 30 | 48 |

The private `candidate_manifest.json` records each candidate SHA-256 and
source-frequency choice. The unchanged CHIANTI input SHA-256 values are:

| Ion | `.elvlc` | `.wgfa` | `.scups` |
| --- | --- | --- | --- |
| O II | `ae169f011f2564828640380b85367ad330c2cff77b06c65249b659a02306685d` | `570ea65bfc4c646fde8f1a7cb008eb9f206088c2eb765fb43c99b3d1d6e23620` | `0614d5b3ebe5e655e82b89a38fc5cfdfb6d294faddc02b592f6bb37be6a78c01` |
| O III | `2bd4da43c855325c8d4c5c6640f539945106219e51231095596797801ec7e297` | `19b9363a630a092eb24963f72ee9b201b8354ed3b0535ec0c4da0c3be922d2e3` | `c45303a6d6a4dd931aece7c1023ac80d19692c35d78de534b4c145ddfe7e9333` |

Finite-density population, emissivity, cooling, stationary balance, and
collision/emission energy closure checks are required before accepting these
private candidates. The relevant thermal domain is 1,000–100,000 K at
electron densities 1e-6, 1e2, 1e6, and 1e10 cm^-3; these tests do not
establish the accuracy of the original atomic collision calculation.

The bounded native calculation completed all 401 temperatures and four
densities for the baseline and both conditional M1/E2 variants of each ion.
`native_result.json`, `native_full_profiles.npz`, and
`same_parity_native_sensitivity.pdf` preserve the full results; the PDF was
rendered and visually checked. The largest population-normalization error
was 1.78e-15 for O II and 4.44e-16 for O III. Relative to the baseline, the
largest total-cooling changes for conditional M1/E2 scaling were
1.1025e-8%/1.8752e-8% for O II and 1.0946e-8%/1.8746e-8% for O III.
These tiny total-cooling sensitivities do not establish the accuracy of an
individual forbidden line probability. Independent source-node, population,
and thermal-balance checks completed in a separate persistent successor
workflow. All original collision x/y nodes were exactly preserved: 286,690
O II and 233,850 O III nodes. Five-temperature native Upsilon comparisons
over every 19,110 O II and 15,576 O III collision row had maximum relative
errors 3.55e-15 and 8.22e-15. At 401 temperatures and four densities, the
independent Python and native cooling differed by at most 2.96e-14 and
2.04e-14 relative, and populations by at most 2.44e-15 and 1.92e-15
absolute. The maximum stationary residual relative to level outflow was
3.70e-15 and 5.33e-15. Collision-to-radiation energy residual relative
to gross collisional power was 5.05e-14 and 5.23e-14; relative to the
smaller net radiation it reached 2.79e-9 and 7.05e-9. These are numerical
balance checks, not uncertainties in the atomic input data.

`independent_balance_result.json`, `independent_balance_profiles.npz`, and
`independent_balance_profiles.pdf` preserve the independent result and full
curves. The PDF was rendered and visually checked. Both persistent
controllers completed with all frozen source and output hashes matching.
`line_sensitivity_result.json` also checks individual emitted lines whose
share of total cooling exceeds 1e-6 at a sampled condition. The largest
conditional M1/E2 relative changes were 5.73e-6/9.64e-6 for O II and
4.54e-5/7.37e-5 for O III. These checks do not cover weaker lines.

The shared reader now applies the identified E1 transformation and the
published O III ground optical rates. A scoped test verifies the original
and corrected low-density coronal cooling are exactly equal on 401
temperatures because both radiation graphs have a unique ground sink.
Thus the installed `cooling_o_2.txt` and `cooling_o_3.txt` coefficients were
preserved. The new full-atom numeric sections are byte-identical to the
privately validated atoms: SHA-256 `b465bc862f5dda7f17fe44ccd793affc57818ab1fc9e8176eab58b07a597bf2f`
for O II and `9e59299ce699b550e6bc31187da2e93b743ff6bdb6229ecfa68f99d335e549b5`
for O III. `installation_manifest.json` gives the archived prior files and
final installed SHA-256 values. An isolated Intel `-O2` direct solver gave
exactly identical cooling and population text values to the previously
validated Intel debug driver at all 1,604 sampled conditions for each ion.
The full installed registered-table consumer was tested separately because
it builds 401 by 161 temperature-density tables, or 64,561 table cells per
ion. The driver checked 40 specified grid and midpoint conditions for each
ion. Both installed atoms passed that consumer with an isolated Intel `-O2`
build from the current unchanged
Fortran source. Registration required 106.45 s for O II and 76.09 s for
O III. Across the driver's exact grid and midpoint samples, the largest
grid-point relative errors were 2.30e-14 and 2.78e-15; midpoint table errors
were 0.1961% and 0.1375%, and final applied-cooling errors were 0.1978%
and 0.1436%. All are below the specified 0.5% midpoint tolerance.
The test also samples numerical table temperatures outside the 1,000–100,000 K
validated domain. Those samples test interpolation behavior only and do not
establish atomic-data accuracy outside that domain.

`registered_O2.txt`, `registered_O3.txt`, and their reviewed result JSON files
contain the measured registration outputs. The original controller source,
state, and JSON results are retained byte-for-byte in
`registry_execution_archive/`, with their hashes verified against the
controller state. A post-run metadata correction renamed the tested domain
field to `validated_numerical_temperature_range_K`; no numerical result was
recalculated. The isolated optimized direct solver and the current Intel
debug-object solver had exactly equal text values for all 1,604 cooling and
population conditions in each ion before the optimized registry build.
`optimized_build_manifest.json` records the compiler version, build commands,
source and binary hashes, installed atomic files, and the distinction between
the optimized registered-table results and the current Intel debug-object
direct profiles. The build used the unchanged Fortran source and left the
production objects untouched.
