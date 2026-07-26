# tools/fitting

Python pipeline: CHIANTI v11.0.2 (`~/RT_Codes/CHIANTI/dbase`) -> rate/cooling
fits. Seeded from EXHALE `cooling_data` (`chianti_cooling.py`). The auditable
source of every coefficient in `data/atomic/` — Fortran never parses CHIANTI
at runtime.

## Scripts

| Script | Product | Role |
|---|---|---|
| `chianti_cooling.py` | (library) | direct CHIANTI ASCII reader (.elvlc/.wgfa/.scups), Burgess-Tully descaling, cooling builders; copied from EXHALE |
| `fit_cooling.py` | `data/atomic/cooling_<ion>.txt` | Tier-1 cooling fits Lambda(T) = T^-1/2 sum A_i exp(-T_i/T), 1e3-1e5 K, in the n_e -> 0 limit of statistical equilibrium; for the G2 thermal loop |
| `fit_nlevel.py` | `data/atomic/nlevel_<ion>.txt` | Tier-2 n-level data: level energies/weights, A-values, Upsilon(T) as Chebyshev fits in Burgess-Tully scaled space; for output-time diagnostics |
| `make_cooling_table.py` | `data/atomic/cooling_fit_parameters.txt` | collects the Tier-1 coefficient files into one listing |
| `verify_nlevel_pyneb.py` | (report) | end-to-end check: reads the fitted files, solves statistical equilibrium as the Fortran will, compares diagnostic ratios against PyNeb |

## Current ion set (2026-07-26)

28 Tier-1 cooling ions and 28 Tier-2 n-level models across the 11-element
registry (H, C, N, O, Ne, S, Ar, Mg, Fe, Si, Cl, Ca) — every cooling ion has
a model, so the runtime density suppression covers all of them.
Adding an ion = add one entry to the `IONS` dict of each fit script and rerun.

## Tier-1 population model

The fitted quantity is the n_e -> 0 limit of the n-level statistical
equilibrium: at n_e << n_crit every ion sits in its lowest level, collisional
excitation happens out of that level only, and each excitation radiates its
full excitation energy through the cascade, so Lambda_ion is the excitation
power out of the ground level.  The T_i of the fit are therefore excitation
temperatures, and `ground_excitation_ladder` derives the starting ladder from
the level energies rather than from hand-tuned values.  Transition energies
are the observed `.elvlc` level differences in both the Boltzmann factor and
the emitted photon; the `.scups` dE is the Burgess-Tully scaling energy and is
used only to build the scaled temperature.

The limit carries no density suppression, so the low-n_crit lines
([C II] 158 um, [O III] 52/88 um, the [O II] and [S II] optical doublets)
are overestimated in dense gas; the density-dependent emissivities used for
the diagnostics come from the n-level solve.

## Fit quality (recorded per file in the provenance headers)

- Tier-1 cooling: max fit error <= 0.7% over 1e3-1e5 K for the H/C/N/O/Ne/S
  set, except C III (5.5% full range, 0.50% where the ion cools) and S II
  (2.2%); the heavily depleted Si, Cl, Ca ions reach 4-28%, concentrated at
  the ends of the fitted range.
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
  A-values there, so this is database spread, not fit error),
  [S IV] 10.51um/1406 within 9% (PyNeb carries DHKD82 collision strengths,
  CHIANTI v11 carries Tayal 2000; the [S IV] 10.51um emissivity itself is
  1.34x PyNeb's, and the Tier-1 [S IV] cooling matches the Cloudy c25.00
  channel of the HII40 benchmark to 3%).

## Requirements

numpy, scipy; PyNeb only for `verify_nlevel_pyneb.py`.
