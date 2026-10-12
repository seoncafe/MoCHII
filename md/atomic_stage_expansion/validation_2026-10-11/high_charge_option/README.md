# High-charge transfer option validation

The production default for `high_charge_exchange` is `.false.`. Explicit
`high_charge_exchange=.true.` enables the adopted q >= 5 unit-transfer
Langevin upper-limit coefficients; these are not measured reaction rates.
The direct population solver, cached rates, hydrogen charge transfer, and
mean kinetic energy use the same selection. Low-charge reactions remain
unchanged by this option.

The persistent controller completed the linked debug build and the native
77-stage/66-transition staged-rate check. Both namelist choices passed
198 reaction/temperature samples. The enabled choice contains 66 high-charge
samples at three temperatures; the disabled choice has zero corresponding
rates and energy transfers. Mg4+ retains its published fit and 8.60 eV mean
heat. The maximum relative error against the independent capture expression
is 3.894832152351313e-7, including printed-output rounding.

Root verified 98 recorded source/product hash entries and both controller
artifact hashes. See [root_review.json](root_review.json) and
[native/checks_result.json](native/checks_result.json). No spatial nebular
comparison or atomic-accuracy validation was performed: this check concerns
the input option and its rate/charge/energy consumers. The staged rate files
have not replaced production data. The full 77-stage line-atom goal and
other additional-physics controls remain incomplete.
