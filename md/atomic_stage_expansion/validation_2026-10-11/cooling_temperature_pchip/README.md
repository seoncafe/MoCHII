# Cooling interpolation at fixed electron density

## Physical and numerical correction

The original cooling table interpolated ln Lambda linearly in inverse
temperature. For N III's expanded source atom, the resulting density ratio
exceeded the existing 0.5% interpolation criterion. Doubling the temperature
grid reduced that discrepancy, but doubles the number of temperature solves.
The production correction instead retains all 401 temperature nodes and
uses shape-preserving cubic interpolation of ln Lambda in -1/T. A linear
Boltzmann exponent is reproduced, positive values remain positive, and
turning profiles remain bounded by adjacent coefficients. At a zero endpoint,
the nonnegative linear rule in log T is retained without taking log zero.
Two-node positive intervals use their secant slope. Source temperature ceilings
and the electron-density interpolation rule are unchanged.

`src/nlevel_cooling_mod.f90` computes the local slopes from existing table
values; no additional table arrays or temperature solves are introduced.
The atomic format, public routines, user inputs, and 401-by-161 grid are
unchanged. This corrects interpolation; it does not add a physical cooling
process or establish atomic accuracy. `metal_freefree=.true.` remains.

## Direct verification and scope

The [actual private-interpolator check](interpolation/result.json) passed
797 samples across eight profiles: a Boltzmann exponent, constant cooling,
a thermal peak, very small positive coefficients, a two-temperature interval,
zero endpoints, an isolated positive pair, and identically zero cooling.
Comparison with SciPy's independent PCHIP and the explicit linear zero limit
gave maximum relative difference 1.7068473908370947e-13. All outputs were
finite, nonnegative, and bounded by the adjacent values. This is not a claim
about subnormal handling under every runtime's underflow configuration.

The full production build passed, as recorded in [build_result.json](build_result.json).
The [controller](controller/state.json) completed all four steps. Actual
thermal/emission consumers passed 1584 thermal samples, 1584 emitted-energy
samples, and 24 zero-density checks for six prepared shared-ion atoms.
The maximum sampled thermal error against the existing fit times an independent
density ratio was 0.023494%; the maximum emitted-energy relative error was
3.2218672174622043e-13. Root verified source and artifact hashes and visually
inspected all six pages of [the full profiles](native/thermal_emission_profiles.pdf).
[Root review](root_review.json) records completion. No production atomic
replacement occurred. These tests compare each calculated quantity with its
own reference; they do not establish equality of an old-source thermal fit
and a different new-source atom's emitted energy. The paired source proposal
addresses that separate issue.

The preceding read-only [54-model assessment](../cooling_interpolation_source_scope/source_interpolation.json)
compared 200 temperature midpoints from 1000 to 100000 K, shortened wherever
the source ceiling requires it, at four fixed electron densities. The maximum
density-ratio difference was 0.009766524503168839 for linear interpolation
(C II), and 0.00453541579242045 for PCHIP (N I). This measures the sampled
temperature interpolation at the density nodes; it is not a continuous error
bound, a test of intermediate electron densities, or a matched nebular run.
It supplements the actual compiled interpolator checks without changing
source rates or regenerating unrelated products.

## Earlier diagnostic source change

The separate 801-node diagnostic generated a numerical result, but its
controller rejected the step when the production interpolation source
changed. Its workflow remains failed. The exact prechange source was
preserved as `../shared_atom_expansion_refined_complete/frozen_nlevel_cooling_mod.f90`
and verified against that controller's startup hash. No failed workflow is
relabeled as successful. The current controller validates the corrected
401-node production method from its own frozen sources and new outputs.

## Subsequent documentation corrections

The module comment was updated to describe PCHIP correctly. Noncomment
Fortran text matches the tested source exactly, so no numerical rerun was
needed. Exact source is preserved in the later paired validation record.
The plotting title was also clarified for future figures: dotted ion-fraction
curves are native low-density results; only the coefficient panels use an
independent reference. Preserved figures were not regenerated.
