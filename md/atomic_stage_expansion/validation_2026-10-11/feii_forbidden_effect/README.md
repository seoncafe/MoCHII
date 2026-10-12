# Fe II forbidden-transition cooling diagnostic

Forbidden Fe II cooling cannot be assumed negligible. At 10,000 K and electron
density 100 cm^-3, adding the 40 radiative probabilities already present in the
installed MoCHII atom increases this controlled Fe II cooling coefficient by
5.332443822766384. This is not the ratio of total nebular cooling or a predicted
change in equilibrium temperature.

## Matched inputs

| Input | Allowed source model | Add installed forbidden transitions |
| --- | --- | --- |
| Levels | Same 142 source levels | Same 142 source levels |
| Energies and statistical weights | Unchanged | Unchanged |
| Collision strengths and interpolation | Same candidate tables | Identical tables |
| Electron density | 1e-6, 100, 1e6, 1e10 cm^-3 | Identical values |
| Temperature | 401 points, 1,000--100,000 K | Identical points |
| Spontaneous probabilities | Source allowed transitions | Add 40 absent pairs from installed 13-level atom |
| Spatial transport and ionization fractions | Not calculated | Not calculated |

The first 13 states match statistical weights and energies to within
0.004 cm^-1. Both calculations retain the 142-level energies; installed A
values are estimates and are not recalibrated from the tiny energy difference.
This controlled difference does not complete all missing Fe II probabilities.

| Electron density (cm^-3) | Cooling ratio at 10,000 K | Cooling ratio at 100,000 K |
| ---: | ---: | ---: |
| 1e-6 | 5.417326614441275 | 1.2442885521915117 |
| 100 | 5.332443822766384 | 1.2394794385077854 |
| 1e6 | 1.0764351424783747 | 1.00448848385113 |
| 1e10 | 1.0000080231682442 | 1.0000004483528693 |

The low-temperature ratio is very large because the allowed-only coefficient
is extremely small. Absolute coefficients, rather than that ratio alone, are
shown in [the full curves](products/forbidden_cooling_comparison.pdf).
All four panels were visually inspected; both 1,604-sample outputs have finite,
normalized populations and matched input columns. [Root review](root_review.json)
records hashes and acceptance. This verifies the numerical diagnostic, not
published absolute atomic accuracy or an entire nebular thermal calculation.

No production Fe II data were replaced. Default selection awaits the user's
choice between explicitly selecting forbidden emission and retaining it as a
physical exception to the requested default-off rule. Metal free-free remains
enabled under the separate explicit user instruction.
