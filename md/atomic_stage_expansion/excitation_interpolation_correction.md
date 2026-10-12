# Bound excitation and nonnegative collision interpolation corrections

## Physical judgment

A nonnegative effective collision strength cannot be represented by a negative interpolated value. Clipping a cubic overshoot to zero removes a real collision channel over a finite temperature interval. The corrected evaluator preserves the original tabulated nodes: it retains cubic interpolation if its exact source-domain extrema are nonnegative and uses PCHIP when nonnegative nodes otherwise produce a negative interior extremum. Type-6 logarithmic ordinates remain cubic. Negative source nodes are preserved for diagnosis; requested temperatures giving negative strengths are rejected rather than assigned a rate floor. This corrects interpolation, not unresolved scattering resonances.

Ordinary bound-state excitation now excludes states above the physical ionization potential and bound states disconnected from the populated collision/radiation network. Autoionization, photoexcitation, excited-state ionization, and line trapping are outside this model. Scattering scaling energies remain distinct from actual excitation and emitted-photon energy gaps.

## Verified production cooling corrections

Both corrections affect source physics, so numerical regeneration was required. Eleven low-density coefficients were generated in staging, compared with direct stationary source curves at 801 temperatures from 1,000 to 100,000 K, checked through linked native evaluation, and installed after inspecting the complete source/fit/previous-product profiles. Original production products are preserved. Installed file hashes match the accepted staged products. The source-fit target remains 1%; the maximum measured error is 0.948463%. The maximum native difference from rounded coefficient evaluation is 9.05210e-12. These are representation errors, not source atomic-data accuracy bounds.

| Ion | Maximum rounded source-fit error (%) | Maximum native relative difference | Full source/fit profile |
|---|---:|---:|---|
| CII | 0.295932 | 4.440892e-16 | [PDF](validation_2026-10-10/bound_cooling_correction/products/CII/CII.pdf) |
| MgII | 0.220633 | 4.440892e-16 | [PDF](validation_2026-10-10/bound_cooling_correction/products/MgII/MgII.pdf) |
| NIII | 0.694403 | 4.440892e-16 | [PDF](validation_2026-10-10/bound_cooling_correction/products/NIII/NIII.pdf) |
| NeIII | 0.335046 | 4.440892e-16 | [PDF](validation_2026-10-10/bound_cooling_correction/products/NeIII/NeIII.pdf) |
| SiIII | 0.438875 | 5.906386e-14 | [PDF](validation_2026-10-10/bound_cooling_correction/products/SiIII/SiIII.pdf) |
| SiIV | 0.895859 | 4.440892e-16 | [PDF](validation_2026-10-10/bound_cooling_correction/products/SiIV/SiIV.pdf) |
| CI | 0.948462 | 4.440892e-16 | [PDF](validation_2026-10-10/positivity_cooling_correction/products/CI/CI.pdf) |
| NI | 0.589366 | 4.440892e-16 | [PDF](validation_2026-10-10/positivity_cooling_correction/products/NI/NI.pdf) |
| SII | 0.379285 | 9.052092e-12 | [PDF](validation_2026-10-10/positivity_cooling_correction/products/SII/SII.pdf) |
| SIII | 0.229389 | 4.440892e-16 | [PDF](validation_2026-10-10/positivity_cooling_correction/products/SIII/SIII.pdf) |
| SIV | 0.608741 | 4.440892e-16 | [PDF](validation_2026-10-10/positivity_cooling_correction/products/SIV/SIV.pdf) |

The comparison holds element, level-energy selection, temperature grid, optically thin radiation, stationary populations, and the zero-density limit fixed. Earlier production coefficients were shown for context; their approximation error is distinct from the change in the direct source. This is a local cooling diagnostic, with no spatial domain, velocity, gravitational potential, or incident spectrum, and is not a nebular code comparison.

## Numerical checks and incomplete native work

The aggregate coefficient writer previously printed eight decimal places for
amplitudes and six for excitation temperatures. Those reductions are unsafe
for signed sums with cancellation. It now retains 16 decimal places, matching
the primary files. The aggregate was rewritten without regenerating primary
cooling products; all 131 coefficients for 28 ions match their primary parsed
values exactly. The previous aggregate and the
[comparison result](validation_2026-10-10/aggregate_coefficient_precision/result.json)
are preserved. This listing is not read by the native cooling loader.

Three physical positivity tests passed: the Ne V 6--15 channel remains positive below the previous artificial zero crossing, original source nodes are retained, and negative type-6 ordinates remain positive physical strengths while a genuinely invalid negative-source region raises an error. Seven stationary cooling tests passed, covering multiple stable states, energy orientation, parser fields, and the independent finite-density limit. No unrelated regression suite or nebular model was run.

The [independent scope audit](validation_2026-10-10/positivity_correction_scope/scope_assessment.md) measures actual changes to direct source curves and lists affected native transitions. C I's earlier corrected five-level native atom passed 3,204 linked checks at 801 temperatures and four electron densities and was installed at that stage. It has since been replaced by the full 42-level bound source; see [compact-atom correction](existing_compact_atom_correction.md). The historical maximum cooling difference is 0.149659%, the maximum absolute population difference is 8.22690e-5, and the maximum rounded transition-fit error is 0.193060%. The source's scaling coordinates remain unchanged; optimized type-2 coordinates apply only to its native numerical representation. Two coordinate-storage tests passed. The complete [native difference profiles](validation_2026-10-10/native_positive_atoms/products/CI/native_differences.pdf) were visually inspected. [Installation hashes and numerical evidence](validation_2026-10-10/native_positive_atoms/products/CI/installation.json) are retained. This validates the existing level set, not full-source truncation or a nebular calculation.

The earlier S III single-series native coefficients required replacement. Full S III source-node storage is now installed and compiled checks passed. Source-node interpolation can exceed the accuracy target of a single 16-coefficient fit, including Ne V and S III. This limitation must not be hidden by accepting the old clipped source or relaxing the transition target. The approved [storage implementation](excitation_storage_design.md) uses direct source nodes and dynamic atoms. Historical single-series failures remain preserved.

## Historical S I single-series limitation and completed replacement

The general cooling fitter now attempts a constrained signed-exponential
minimax expansion only after its existing shorter expansions fail. The
excitation base is derived from the physical source ladder; numerical basis
offsets are not additional atomic levels. Original successful fits retain
their existing calculation path. Rounded coefficients must pass 801- and
1,601-temperature source checks, sampled positivity at 10,001 temperatures
from 1 to 1e9 K, and a sampled cancellation condition below 1e5. Positivity
outside the source interval does not establish extrapolation accuracy.

S I generation through the normal fitter produced 21 nonzero terms, with
maximum source-fit error 0.496964521%, maximum native/source error
0.496964521%, and maximum native/rounded-coefficient difference 2.818e-12.
The largest sampled cancellation condition is 5.419e4. Its five bound states,
nine radiative rows, and ten collision rows are retained in the source.
The [completed calculation controller](validation_2026-10-10/si_cooling_fit_diagnostic/normal_generator_workflow/progress.md)
ran generation followed by validation. The
[main-session review](validation_2026-10-10/si_cooling_fit_diagnostic/normal_generator/main_review.json)
records original source hashes and the accepted staged product hash.
The complete [21-term source and residual profiles](validation_2026-10-10/si_cooling_fit_diagnostic/normal_generator/si_cooling_profile.pdf)
were visually inspected by the main session. The
[detailed assessment](validation_2026-10-10/si_cooling_fit_diagnostic/si_cooling_fit_assessment.md)
records the controller's execution limitation and the subsequent comment-only
hash change; no numerical regeneration followed that wording edit.

The attempted five-level native atom failed before its population solve.
Transitions 3--4, 3--5, and 4--5 contain consecutive zero source nodes
followed by positive values; the adopted PCHIP source has exact zero
intervals. A single finite polynomial cannot reproduce these intervals.
The [source-zero diagnostic](validation_2026-10-10/si_cooling_fit_diagnostic/si_native_collision_zeros.json)
identifies them explicitly. No S I product was installed at that earlier
failed single-series step. The approved source-node representation subsequently
preserved these zero intervals, passed native population and registered
cooling checks, and was installed together with the accepted cooling fit.
See [storage implementation](excitation_storage_design.md).

The approved 64-stage roster is still staged. These corrections do not implement nonadjacent Auger balance, terminal-stage ionization losses, or the expanded production roster.
