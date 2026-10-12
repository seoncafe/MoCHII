# Positive interpolation of bound collision strengths

The current CHIANTI scaled-strength cubic spline can become negative between
nonnegative source nodes. Clipping the result to zero removes a physical
collision channel over finite temperature intervals. A positive source-node
interpolant is therefore needed before finalizing native fits for affected
transitions. A 16-coefficient fit of the clipped curve, documented in
`diagnosis.md`, repairs representation error only and is not a source-physics
correction.

## Exact inventory

`positivity_inventory.py` uses the approved 64-stage limits matching
`validate_terminal_ionization.py`, calls the current
`bound_excitation_atom()` selection for each available ion, and checks every
stationary point of each cubic polynomial inside its tabulated scaled-
temperature interval. It also checks extrema within the image of
1000–100000 K. Comparison with SciPy's independent derivative-root finder on
150 sampled transitions gave a maximum difference of `1.92e-15` in scaled
strength. All affected transition identities and extrema are in
`positivity_inventory.json`.

Of 64 stages, eight have no excitation source: Ar I/II, Fe I, Cl I/V, and
Ca I/III/IV. The other stages contain 331,963 bound collision transitions:
313,701 CHIANTI scaled-strength transitions and 18,262 Stout transitions.
The Stout transitions use an independent temperature-table evaluator and are
excluded from the cubic-spline diagnosis. Two CHIANTI type-6 transitions,
O I 2–3 and 4–5, store logarithmic scaled ordinates; negative ordinates are
expected there and cannot be classified as negative collision strength.

Among the 313,699 remaining CHIANTI transitions, 1,397 have nonnegative
source ordinates but an exact negative interior cubic extremum. These occur
in 24 ions. Fifty have negative extrema within the 1000–100000 K fitting
interval. Their identities, grouped by ion, are:

| Ion | Negative cubic extrema in fit interval | All negative interior extrema |
| --- | ---: | ---: |
| C II | 7 | 8 |
| N I | 9 | 64 |
| Ne V | 2 | 59 |
| S I | 3 | 9 |
| S II | 15 | 155 |
| S IV | 1 | 3 |
| Ar VI | 1 | 97 |
| Mg IV | 1 | 2 |
| Fe II | 4 | 14 |
| Fe III | 7 | 166 |

The other affected ions have negative extrema outside the fitting interval:
C I (70), C IV (2), C V (89), N III (4), N IV (287), N VI (36), O II (22),
O III (26), O IV (1), O V (1), S III (110), Ar IV (4), Ar V (150), and
Si I (18). The full 50-pair list and source hashes are in the JSON; the
Ne V pairs within the fitting interval are 6–15 and 6–20. Some excursions
are extreme: S II 41–45 reaches `-624.17` between positive table nodes.
These numbers are exact extrema of the implemented interpolant, not sampled
approximations.

Three further type-2 transitions already have a negative source ordinate at
their final scaled-temperature node: Mg I 2–12 (`-0.009492`), Mg I 13–15
(`-0.0007207`), and Fe II 1–9 (`-0.001309`). None has a negative spline value
in the 1000–100000 K interval. PCHIP through these signed nodes cannot make
the source nonnegative. Their high-temperature source limits require a
separate source decision. The inventory excludes type-6 logarithmic
ordinates from all three counts.

## Ne V 6–15 full profile and fitting

The 15 CHIANTI nodes for this transition are nonnegative. Its ordinary cubic
spline reaches `-8.3216173394e-5` near 1729.12 K and crosses zero at
3622.576354 K. The current evaluator clips this to zero at 448 of 1601
logarithmically spaced temperatures from 1000 to 100000 K. PCHIP through the
same original nodes is nonnegative throughout the tabulated domain; it gives
`7.7375955e-6` at 1000 K, `6.1139165e-5` at 3000 K, and
`8.5378893e-5` at the cubic zero crossing. It gives `9.9853883e-4` at
10000 K, compared with `9.9413190e-4` from the clipped cubic. The full
curves and their fractional difference appear in `nev_pchip_vs_cubic.pdf`,
which was rendered and visually inspected.

`nev_pchip_profile.py` minimizes relative error of a degree-15 Chebyshev
polynomial on the 1601-point PCHIP curve, testing both linear and logarithmic
storage and 42 positive type-2 coordinate constants from approximately
`1e-3` to `1e2`, including the original `C = 0.9492`. The best tested
coordinate is `C = 0.023713737`, with 1.05836% maximum relative error on
1601 points and 1.05868% on 160001 points after native field rounding. With
the original coordinate, the best logarithmic fit has 4.37261% error.
Thus the current one-piece 16-coefficient representation does not meet the
0.2% target for PCHIP in this coordinate scan. This is not a proof that no
other positive interpolation or coordinate can meet that target.

As an exploratory smooth alternative, a cubic spline with its first-node
derivative set to the first secant slope is nonnegative across the entire
tabulated domain for this transition; its only global minimum is the zero
first node. Its best tested 16-coefficient logarithmic fit has 0.13920%
error on 1601 points and 0.13930% after rounding on 160001 points.
The first-secant boundary condition is supported by adjacent source nodes,
but its applicability to the other 1,396 overshoot transitions has not been
established. It should not be applied indiscriminately as a global fix.

The calculations in this directory are diagnostics only. No shared evaluator,
source table, production atom, or native data layout was changed.
