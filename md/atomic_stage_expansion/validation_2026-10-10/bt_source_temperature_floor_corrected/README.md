# BT source-temperature continuation validation

The optional `BT_TMIN_K` record holds the source effective collision strength
constant below the identified lower source temperature. Both the Burgess-Tully
coordinate and its inverse scaling use `max(T_actual, T_min)`. Electron-impact
rates, Boltzmann factors, and statistical equilibrium use `T_actual`. The
original scaled nodes, scaling energy, and collision type remain unchanged.

| Ion | Minimum source T (K) | Selected levels | Radiative rows | Collision rows | Maximum native/source cooling relative error |
| --- | ---: | ---: | ---: | ---: | ---: |
| C IV | 3,200 | 24 | 74 | 66 | 1.78e-15 |
| Mg III | 1,800 | 27 | 117 | 351 | 2.00e-15 |
| Si V | 5,000 | 169 | 3,778 | 14,196 | 2.74e-14 |

Each actual-source comparison used 801 temperatures from 1,000 to 100,000 K
and electron densities of 1e-6, 1e2, 1e6, and 1e10 cm^-3. The complete
temperature profiles are in each ion's `cooling_profiles.pdf`; all four Si V
profiles were visually inspected and the compiled curves overlap the source
curves. These comparisons check numerical agreement with the selected source
atom, not the physical completeness of its radiative data. Mg III uses the
first 27 states because the remaining radiative rows have unresolved source
direction and statistical-weight consistency.

`result.json` checks two-level atoms for all six BT types at the minimum and
neighboring temperatures, including detailed balance with the actual
temperature. Its largest relative cooling error is 6.67e-16. Eleven malformed
metadata cases fail with explicit native reader errors. Old version 2 files
without the record and legacy Chebyshev files retain their original behavior.
The source hashes in `result.json` and `actual_source_workflow/state.json`
identify the build used for these positive-temperature comparisons. A later
Python guard rejects nonfinite and nonpositive requested temperatures; its
focused tests passed, and the positive-profile products were preserved.

The source-temperature continuation is limited to the three identified ions
and the source ranges recorded in their generated atom headers. This directory
contains validation products only; it does not install production atoms.
