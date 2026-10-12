# Terminal charge diagnostic at original recombination nodes

The production target remains 64 stages and 53 adjacent transitions. This
diagnostic changes no production data, charge roster, temperature limit, or
acceptance tolerance. It separates a recombination fit discrepancy from the
physical approximation of omitting higher charge states.

## Calculation and validity

For a stationary adjacent electron-only chain with exactly one additional
state, let `r = C_terminal / alpha_reverse`. If the retained chain has total
unnormalized population `S` and terminal weight `w <= S`, the additional state
has fraction `r*w/(S+r*w) <= r/(1+r)`. The cancellation of electron density
requires both rates to be proportional to that density. Neither photons nor
charge exchange nor nonadjacent transitions are included. This inequality
does **not** bound all further omitted stages or thermal errors.

The calculation uses Voronov terminal ionization and the original ground-state
Badnell RR/DR columns at their exact common printed temperature nodes. It does
not interpolate those original columns or extend their coverage. The terminal
CI HYBRID correction is uncalibrated. Fe VIII requires its separate experimental
DR source and was not evaluated in this original-node comparison.

The directly executed results are in [native_node_bounds.json](native_node_bounds.json).
The recorded maxima refer only to the sampled nodes, not to the full
1,000--100,000 K interval:

| Element | Common nodes | Sampled temperature range (K) | Maximum one-state fraction bound | Maximum total reverse fit discrepancy |
| --- | ---: | ---: | ---: | ---: |
| C | 1 | 50,000 only | 1.95022e-39 | 0.465314% |
| N | 1 | 72,000 only | 4.99525e-39 | 0.431276% |
| O | 6 | 1,250--50,000 | 5.64753e-11 | 2.54968% |
| Ne | 6 | 1,250--50,000 | 8.62740e-12 | 2.40344% |
| S | 6 | 2,450--98,000 | 1.15534e-13 | 1.76119% |
| Ar | 6 | 2,450--98,000 | 3.25734e-6 | 5.91499% |
| Mg | 6 | 1,600--80,000 | 6.41067e-5 | 0.864961% |
| Si | 6 | 1,250--50,000 | 3.21167e-15 | 0.668296% |
| Cl | 6 | 2,450--98,000 | 1.30820e-4 | 1.49294% |
| Ca | 6 | 1,800--72,000 | 4.83056e-7 | 2.56585% |

The discrepancy compares the evaluated analytic total reverse coefficient to
the rounded original RR-plus-DR columns. It is not an estimate of atomic
calculation accuracy. In particular, the O/Ne maxima at 50,000 K cannot replace
the prior analytic-fit diagnostics at 100,000 K. A fit's stated accuracy range
and the existence of original computed nodes are different kinds of evidence.

The script checks positive finite total reverse rates, nonnegative original
DR values, exact common temperature identities, finite bounds in `[0,1]`, and
source hashes. It refuses to overwrite existing results. Execution completed
in a short local run; no persistent calculation controller was started. No
whole-code or nebular comparison was performed because no production routine
was changed and this model excludes terms required for that comparison.
