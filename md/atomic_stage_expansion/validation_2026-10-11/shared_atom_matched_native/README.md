# Actual combined atom and cooling-fit invocation

The finite controller completed both the actual thermal/emission calculation
and the figure step. The twelve paired candidate files were read in isolated
fixtures. The 33 temperatures span 1000-100000 K at ne=1e-6, 100, 1e6, and
1e10 cm^-3, with local-density and low-density thermal choices. All 1584 thermal
samples, 1584 emitted-energy samples, and 24 zero-density checks passed.
The maximum local-density thermal/direct-emission discrepancy is
0.0045777284511308025 (0.457773%), for Si II, within the 0.5% representation
criterion. This independently confirms the preceding derived predictions by
actually invoking the combined candidates. Production was not modified.

Root verified recorded sources and artifacts and visually inspected all six
pages of [the full profiles](native/thermal_emission_profiles.pdf). There was
no clipping or overlap. The title's independent dotted curves refer to the
thermal and emitted panels; lower-right dotted curves are native low-density
ion fractions, not independently validated ionization. That wording has been
corrected in the plotting tool for future figures without rewriting the
preserved figure. Exact source copies precede the subsequent comment/caption
correction. See [numerical result](native/result.json) and
[root review](root_review.json). Identical source atoms and temperature/density
conditions are used by the independent matrix. Current imposed photo rates
supply stage fractions, which are divided out; no spatial, velocity, spectral,
boundary-condition, or atomic accuracy equivalence is claimed.

The companion [O I hydrogen-excitation check](../shared_oi_hydrogen/README.md)
also passed. All prerequisites are concrete and reviewable; production paired
replacement remains subject to the recorded approval requirement.
