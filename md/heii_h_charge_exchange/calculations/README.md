# Rate Formula Comparison

`rate_comparison.py` evaluates the radiative expression and the KF96 term
inspected in `src/ion_balance_mod.f90` on 2026-09-25. It uses Python's
standard library and writes CSV to standard output:

```bash
python md/heii_h_charge_exchange/calculations/rate_comparison.py
```

[rate_comparison.csv](rate_comparison.csv) is the recorded output. The
script does not modify simulation inputs or products.

## Columns

| Column | Meaning |
|---|---|
| `T_K` | Gas temperature in K |
| `kf96_usage` | Whether the code extrapolates below the fitted range, evaluates within it, or clamps to the upper limit |
| `radiative_cm3_s` | `1.25e-15 * (T/300)^0.25` |
| `nonradiative_as_implemented_cm3_s` | KF96 term with `t4 = min(T, 1e5)/1e4` |
| `total_as_implemented_cm3_s` | Sum of the preceding two rates |
| `total_with_quarter_nonradiative_cm3_s` | Radiative rate plus 0.25 times the KF96 term |
| `quarter_total_ratio` | Ratio of the sensitivity scenario to the current total |

The fit interval is the KF96 interval, not a validated interval for the
separate radiative fit. The quarter-scaled rate is hypothetical. Its values
do not establish that a correction should be applied.

This script is a frozen transcription of the formulas, not a Fortran test,
a simulation result, or an independent cross-section calculation. It does
not automatically follow future source edits. Inspect the live source
before reusing it for a later implementation.

The rough Figure 8 integration discussed in the review is excluded from
these data. It used visually estimated and assumed inputs and did not
establish a physical rate suitable for adoption. A new calculation requires
verified cross sections and controlled interpolation and integration.
