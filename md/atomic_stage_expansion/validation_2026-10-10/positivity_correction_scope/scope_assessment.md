# Reach of selective positive collision-strength interpolation

The current source evaluator selects PCHIP for the 1,397 bound CHIANTI
transitions whose nonnegative scaled-strength nodes nevertheless produce a
negative interior cubic-spline extremum. An independent scan of every
available stage in the approved 64-stage roster selected exactly the same
1,397 transitions: no extras or omissions. Type-6 logarithmic ordinates
remain cubic. The three transitions with a negative final source node remain
cubic and now reject a requested negative effective strength rather than
silently clipping it. This report measures affected source and installed
native atoms; it does not change either.

`assess_scope.py` compares a separately implemented previous
`CubicSpline` plus zero-clipping evaluator with the current source evaluator.
It uses all 1,397 selected transitions, including those whose negative cubic
region is outside 1000–100000 K. The low-density calculation calls the
implemented `cooling_low_density_limit()` on 801 logarithmically spaced
temperatures in that interval. This checks the radiative-cascade and
collision-closed-state source path, not the fitted cooling tables or a full
thermal-equilibrium simulation. Exact values and the complete file/pair list
are in `scope.json`.

Thirty installed `data/atomic/nlevel_*.txt` files were inspected. Four
contain transitions selected by the policy:

| Installed atom | Included transition identities | Reason |
| --- | --- | --- |
| C I | 1–2, 1–3, 2–3 | Nonnegative source nodes; negative cubic extremum |
| S III | 2–6, 3–6, 4–6 | Nonnegative source nodes; negative cubic extremum |
| Fe II ordinary | 1–9 | Negative final source node |
| Fe II full | 1–9 | Negative final source node |

The two Fe II files require attention to the new invalid-temperature
rejection rule, although their 1–9 source strength stays nonnegative on the
evaluated 1000–100000 K interval. The source-level rule also affects atoms
not currently installed; a native file for Ne V, for example, would include
its selected 6–15 and 6–20 collisions if those levels are retained.

The production `fit_cooling.IONS` registry has 28 ions. Twelve contain at
least one selected source transition and were computed with both evaluators;
the other 16 have no selected source transition. Five of the twelve had a
nonzero difference in the 801-point steady low-density cooling curve:

| Ion | Largest signed relative change | Temperature [K] | Largest absolute change [erg cm^3 s^-1] |
| --- | ---: | ---: | ---: |
| C I | -6.90496% | 1251.70 | `1.02077e-23` |
| N I | +0.199408% | 100000 | `4.76066e-22` |
| S II | +4.53132% | 100000 | `5.38257e-20` |
| S III | -0.000169229% | 91727.6 | `1.51109e-24` |
| S IV | +0.584073% | 100000 | `7.45022e-21` |

Each absolute maximum may occur at a different temperature from the relative
maximum. The other seven measured registered ions had bitwise identical
cooling arrays on this grid. This can occur when affected pairs do not
connect a populated radiatively closed state or where their contribution is
outside the sampled temperature range.

Thirteen additional approved stages contain selected source pairs but are
not in the current cooling registry. Five have nonzero 801-point source
changes: Ar VI (+`5.25e-8` relative), C V (+1.3011%), N IV (+1.4053%),
N VI (-3.5491%), and S I (+30.1948%). The C V and N VI relative maxima
occur where baseline cooling is only `3.71e-65` and `2.51e-103`
`erg cm^3 s^-1`; their largest absolute changes are `1.25e-36` and
`3.45e-42` in the same units. S I and N IV have much larger absolute changes,
`2.01e-22` and `6.69e-23`. The remaining eight additional stages had
identical cooling arrays on this grid. Ne V's strict low-density cooling is
unchanged because the selected 6–15 and 6–20 collisions do not contribute
from a populated radiatively closed state in that limit; finite-density
populations still require a separate validation.

Only this diagnostic directory was written. Installed native and production
cooling files retain their previous numbers until independently regenerated
and validated after the source-correction design is settled.
