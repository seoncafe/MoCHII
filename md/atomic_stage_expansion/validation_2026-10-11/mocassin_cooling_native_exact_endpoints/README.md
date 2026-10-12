# Compiled cooling reader validation

52 candidate fits passed the actual linked production `cooling_load` and
`cooling_eval` routines at 801 temperatures each: 41,652 scalar samples.
The measured maximum relative difference from the strict supplied-model
low-density source is 0.0049953117071.
The acceptance condition is less than 0.5 percent within the declared
numerical interval. Signed terms describe total cooling; their cancellation
condition is bounded by 1e5 on the dense validation grid. Absolute atomic
accuracy and extrapolation outside the declared interval are not certified.

Root verified all controller source hashes, the result artifact, and every
native output hash. See [root_review.json](root_review.json) and
[native/result.json](native/result.json). The source fits are listed in
[the collected manifest](../mocassin_cooling_native/manifest.json).
No candidate file has been installed in production. Consumers must enforce
the numerical temperature interval before these fits can be installed.

Four coefficient preparation series preserved all successful earlier
products and retried only unfinished ions. The series containing unfinished
fits remain marked failed; they are not reported as complete workflows.
Ne VI and Fe V still have no accepted exponential fit. Mg I and Ca VII
are separate atomic-source requirements, excluded from these 54 coherent
source models.

The first native-reader attempt failed at the Fe II table ceiling because
its computed final temperature exceeded the supplied limit by roundoff.
The diagnostic now writes the exact first/last temperatures. Failed products
remain in `../mocassin_cooling_native/`. The shared Python temperature-table
error previously named Cl VII and its interval for every caller; its message
now reports the actual supplied ceiling, without changing interpolation.

Validation uses fixed-ion low-density excitation with electron collisions
and spontaneous cascades. It is not a spatial nebular comparison or a
production 77-stage acceptance check. The 142-level Fe II candidate remains
an allowed-line model with the forbidden-transition omission documented in
[the source-atom report](../mocassin_source_atoms/README.md).
