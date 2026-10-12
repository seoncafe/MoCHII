# Candidate MOCASSIN finite-level atoms

These 54 files contain supplied MOCASSIN v2 levels, positive radiative
probabilities, and effective collision-strength nodes in native version-2
format. The packed Fe II emitting direction is corrected. Original source
files remain unchanged. Mg I and Ca VII are excluded for identified source
limitations; `manifest.json` records the reason.

The actual native reader and downstream level, line, and cooling consumers
passed 5,408 temperature/density cases. See
[the numerical validation](../../validation_2026-10-11/mocassin_source_atoms/README.md).
Positive log-log PCHIP differs from the original natural cubic interpolation;
`interpolation_comparison.json` quantifies individual strength differences.
Below the first supplied node the user-adopted constant-strength approximation
is used with actual Maxwellian temperature. No extrapolation above the final
node is provided. Supplied table bounds are not published accuracy limits.

This directory preserves the 54 original candidate files. Twenty-three
previously absent atoms have now been installed byte for byte; existing shared
models were not replaced. See [the installation report](../../validation_2026-10-11/missing_atom_installation/README.md).
These files do not establish completed exponential cooling fits. Native numerical validation does not certify atomic accuracy. All 54 source-profile pages have been visually inspected. Physical
review, default/shared-ion model selection, and production integration remain
necessary. When no cooling fit is supplied, the production consumer can now
use absolute atom cooling; see [the linked validation](../../validation_2026-10-11/atomic_cooling_without_fit/README.md).
