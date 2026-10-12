# Sulfur charge transfer and the selected product-ion excitation model

## Physical judgment

The selected S II excitation atom must determine the product-term populations
used in the density-dependent reverse charge-transfer rate and its mean
kinetic-energy exchange. A selector used only for cooling and emitted lines
would leave an inconsistent charge-transfer calculation. Larger source models
cannot be assumed to have negligible effects without measurement.

Actual native calculations compared the supplied v2 5-level S II atom with
the installed 70-level atom, using exactly the same sulfur reaction file.
Each model passed 404 samples: 101 temperatures from 1000 to 100000 K at
four electron densities. An independent sum of the forward term rates,
detailed-balance reverse rates, and direct statistical-equilibrium term
populations reproduced the actual `species_cx_coefficients` and
`species_cx_kinetic_energy` results.

| Quantity | Measured result |
| --- | ---: |
| Native samples, both models | 808 |
| Maximum reverse-rate check error, 5 levels | 1.58651e-13 |
| Maximum reverse-rate check error, 70 levels | 1.09246e-13 |
| Maximum mean reverse-energy check error | 5.36793e-14 eV |
| Installed/reference reverse-rate ratio range | 0.7563384623-1.4464932194 |
| Maximum relative reverse-rate source difference | 44.64932194% |
| Largest difference location | 19952.62315 K, 100 cm^-3 |
| Maximum mean reverse-energy source difference | 0.2184102926 eV |

The forward rates are identical in both calculations. Reverse-rate differences
include both collision-data changes and additional excitation levels; this
experiment does not isolate the effect of adding levels alone. It does not
measure a total thermal cooling difference or a spatial nebular response.
No atomic-accuracy claim or physical-default selection follows from this
numerical comparison.

## Product terms and input conditions

The existing `CXSTATEK` ranges identify 4S (level 1), 2D (levels 2-3), and
2P (levels 4-5) in both source atoms. Statistical weights are exactly the
required 4, 10, and 6. Maximum weighted term-energy differences from the
reaction file are 4.48108e-6 eV for the supplied 5-level model and
2.65616e-7 eV for the installed 70-level model. Both satisfy the actual
2e-3 eV reader tolerance; no remapping of these three ranges is needed for
these particular candidate files.

| Input condition | Both comparisons |
| --- | --- |
| Reaction file | Identical installed `element_s.txt` |
| Resolved channels | Identical forward rates into 4S, 2D, and 2P |
| Temperature samples | 101 values, log10(T/K)=3+0.02*i, i=0-100 |
| Electron density | 1e-6, 100, 1e6, 1e10 cm^-3 |
| Excitation calculation | Electron collisions and radiative de-excitation |
| Changed input | Product-ion atom: supplied 5 levels versus installed 70 levels |
| Source accuracy | Unverified; source-data difference is unavoidable and explicit |
| Radiation spectrum, boundaries, gravity, and spatial domain | No photoionization or spatial model |
| Charge-transfer heating in thermal balance | Not invoked; diagnostic mean energies only |

The 70-level atom's three represented terms have a minimum combined population
of 0.9904478973 over the samples. The atom remains normalized over all levels;
renormalizing just these three terms would overstate the represented reverse
rate. Rates from the other terms are absent from this reaction table, so this
is not evidence that those physical channels are zero or complete. The existing
implementation correctly retains this distinction. Its inaccurate comment
claiming that the represented term populations sum to one was corrected.
Noncomment production statements were verified unchanged; no additional
numerical run was needed for that comment correction.

## Evidence and current state

The controller `mochii_sulfur_cx_models_20261011` completed with finite time
limits. Root verified its source hashes and the native artifacts before the
comment correction. The exact tested physics source is retained in
[native/species_mod.f90](native/species_mod.f90). See
[native results](native/result.json), [controller state](controller/state.json),
and [root review](root_review.json).

All four pages of [full temperature profiles](charge_transfer_profiles.pdf)
were visually inspected. Each shows the complete sampled reverse-rate,
mean-energy, and three-term population profiles at the stated density.
The smallest-density mean-energy panel uses an explicit offset because both
models approach the same ground-term energy.

Production atomic files, input syntax, and rate equations were not changed.
The public excitation-model selector remains unimplemented and needs the
previously requested interface decision. `metal_freefree=.true.` is retained.
This qualification does not complete the 77-stage production objective.
