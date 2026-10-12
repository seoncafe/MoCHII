# Diagnostic selection migration

Three native diagnostic drivers now read the same explicit excitation model
as the registered cooling table. The obsolete Fe selector was removed. The
single-atom registered-cooling validator creates an isolated model directory
without regenerating the source atom or cooling fit. Existing saved-output
validation retains the original source-file convention.

All three updated drivers compiled against the current production modules.
No numerical table, plot, or unrelated regression output was regenerated for
these path changes. Historical driver text before the public-input migration
is preserved with hashes. The first compile attempt encountered an unrelated
compiler-incompatible module in shared `/tmp`; an isolated module directory
resolved it. The prior full production build and relevant model-consumer
validation remain the numerical evidence for the current implementation.
