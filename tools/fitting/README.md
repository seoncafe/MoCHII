# tools/fitting

Python pipeline: CHIANTI v11.0.2 (`~/RT_Codes/CHIANTI/dbase`) -> rate/cooling
fits. Seeded from EXHALE `cooling_data` (`chianti_cooling.py`). The auditable
source of every coefficient in `data/atomic/` — Fortran never parses CHIANTI
at runtime.

## Scripts

| Script | Product | Role |
|---|---|---|
| `chianti_cooling.py` | (library) | direct CHIANTI ASCII reader (.elvlc/.wgfa/.scups), Burgess-Tully descaling, cooling builders; copied from EXHALE |
| `fit_cooling.py` | `data/atomic/cooling_<ion>.txt` | Tier-1 cooling fits Lambda(T) = T^-1/2 sum A_i exp(-T_i/T), 1e3-1e5 K, ground-term Boltzmann population, low-density limit; for the G2 thermal loop |
| `fit_nlevel_tier2.py` | `data/atomic/nlevel_<ion>.txt` | Tier-2 n-level data: level energies/weights, A-values, Upsilon(T) as Chebyshev fits in Burgess-Tully scaled space; for output-time diagnostics |
| `verify_nlevel_pyneb.py` | (report) | end-to-end check: reads the fitted files, solves statistical equilibrium as the Fortran will, compares diagnostic ratios against PyNeb |

## Current ion set (2026-07-11)

O II, O III, N II, S II, S III, Ne II, Ne III (docs/PLAN.md section 8 order).
Adding an ion = add one entry to the `IONS` dict of each fit script and rerun.

## Fit quality (recorded per file in the provenance headers)

- Tier-1 cooling: max fit error <= 2.3% over 1e3-1e5 K (most < 1%).
- Tier-2 Upsilon: worst transition <= 1.7%; the residual outliers are weak
  intercombination transitions (O III 1D/1S - 5S2) and the coarse 10-point
  S II .scups tables.
- End-to-end vs raw CHIANTI (same solver, fitted file vs direct spline):
  [O II] 3726/3729 agrees to 0.05%, [S II] 6717/6731 to 0.5% — the fitted
  files are faithful to CHIANTI v11.
- Cross-check vs PyNeb 1.1.31 (independent atomic data, not CHIANTI):
  [O III] 5007/4363 within 1.2%, [N II] 6584/5755 within 0.3%,
  [O II] 3726/3729 within 6%, [S II] 6717/6731 within 10% (worst at
  n_e = 1e4 cm^-3, near the critical density, where the ratios become
  A-value-limited; PyNeb uses Kal09/TZ10 collision data and Z82-WFD96/RGJ19
  A-values there, so this is database spread, not fit error).

## Requirements

numpy, scipy; PyNeb only for `verify_nlevel_pyneb.py`.
