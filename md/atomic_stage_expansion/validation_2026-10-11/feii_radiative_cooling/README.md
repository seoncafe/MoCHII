# Fe II source-consistent cooling representation

The candidate combines the 142 supplied MOCASSIN v2 Fe II levels with 1144
positive allowed-transition probabilities and the 40 existing MoCHII forbidden
probabilities. It is a prepared physical model, not an installed production
replacement or a complete inventory of all Fe II forbidden transitions.

## Matching low-density coefficient

The 29-term exponential representation uses the same energies, collision
strengths, and radiation probabilities as the candidate atom. Its maximum
relative error against the independent strict low-density statistical-equilibrium
solution is 0.0299785842% on 1601 temperature samples from 1000 to 100000 K.
The actual Fortran cooling reader and evaluator passed 801 samples, with a
maximum relative source error of 0.0299785841%.

The measured maximum cancellation factor is 99999.99999911059, close to the
existing acceptance limit of 100000. This sampled check is not a continuous
bound or an atomic-accuracy measurement. The fit is qualified only over the
stated temperature interval and is not evidence for extrapolation beyond it.

## Evidence and limits

[Preparation](prepare/result.json), [native evaluation](native/result.json),
and [root review](root_review.json) preserve the measured results and hashes.
The exact tool used for these completed calculations is retained in
[frozen_preparation_tool.py](frozen_preparation_tool.py). Subsequent source-hash
checks and finite-density diagnostics did not regenerate these products.
The [full low-density profiles](prepare/low_density_profiles.pdf) were visually
inspected: the coefficient and its relative fitting error are shown over the
full interval.

The 1604 values in `derived_thermal_profiles.txt` combine the qualified fit
with previously computed finite-density ratios. Their maximum derived
thermal/emitted-energy difference is 0.0299785625%. These are derived values,
not new native thermal-balance invocations. Actual matching-fit application is
tracked separately in [the finite-density continuation](../feii_density_cooling_checked/README.md).

The physical default choice remains unaccepted. No production atom, cooling
coefficient, public input, or data layout changed. This work does not establish
full spatial agreement with MOCASSIN v2 or complete the 77-stage objective.
