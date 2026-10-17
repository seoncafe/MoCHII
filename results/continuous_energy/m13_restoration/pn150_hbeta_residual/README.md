# PN150 L(Hbeta) residual against Cloudy c25.00 (2026-10-14, extended 2026-10-17)

MoCHII `../pn150_p27_nodrift` gives L(Hbeta) 0.893 of Cloudy's. The
findings and the plan are in
`md/METAL_RECOMBINATION_DIFFUSE_PLAN_2026-10-14.md` (Sections 1.2-1.6).

Result of 2026-10-17: Cloudy's heavy-ion recombination photons (one at 0.75
and one at 0.25 of the IP for every metal recombination) have the same
effect in both codes -- R_B(H)/Q +0.074 in Cloudy (`faint` against
`nofrac`), +0.073 in MoCHII (`l6c_metals` against `l6vf_metals`); L(Hbeta)
+7.1% and +7.3%. Without them Cloudy gives R_B(H)/Q 1.0444 and MoCHII
1.0068, so the remaining 3.6% is a separate metal-related difference (with
the metals off the codes agree to 0.8%).

| script | what |
|---|---|
| `budget.py` | MoCHII vs Cloudy: ionized volume, EM, R_B(H), Hbeta per recombination, He recombinations |
| `budget_runs.py <base>...` | R_B per Q(H) for MoCHII runs |
| `metal_absorption.py`, `absorption_runs.py <rates.h5>` | photoionizations by absorber (MoCHII) |
| `cloudy_absorption.py` | the same for Cloudy, from `gam.gam` (save gammas) |
| `euv_lines.py` | MoCHII metal lines shortward of 911.75 A |
| `cloudy_variants.py <base>...` | R_B per Q, L(Hbeta), printed EUV metal lines for Cloudy decks |
| `final_solve_residuals.py <base>...` | leaves the final solve moves, their R_B(H) share, leaves above 30000 K |
| `compare_l6.py --base pn150_hbeta_residual/<run>` | the `../compare_pn150.py` line and front tables for a run here (`compare_<run>.txt`, `<run>_profiles.pdf`) |
| `valence_photon_budget.py [rates.h5...]` | Cloudy's heavy-ion recombination photons per Q, evaluated on Cloudy's structure and on MoCHII states |
| `ledger_table.py <log>...` | the absorption ledgers (`par%absorption_ledger`) per Q |
| `vf_check/channel8_check.py <rates.h5>` | independent channel-8 luminosity |
| `make_cloudy_state.py <template> <out>` | Cloudy's converged structure as a MoCHII restart file (`cloudy_state_rates.h5`) |
| `supply_demand.py <rates.h5>...` | ionizations supplied by the field against recombinations demanded by the state, total and by radius |

Cloudy decks (`*.in`, run with `~/CLOUDY/c25.00/source/cloudy.exe -r <name>`):
`base`/`faint` (the reference deck), `iter` (iterate to convergence),
`nometal`/`nometal_f` (C N O Ne Mg Si S off), `nometal_ems` (Hbeta
emissivity profile), `nocx`, `noind`, `nolev2`, `noots`, `nooutl`,
`nooutc`, `nometal_noots`, `nometal_nooutc`, `gam`. `nofrac` is `faint` run
with a diagnostic build of Cloudy c25.00 in which `ion_recomb.cpp`
`FRAC_LINE` is 1e-20 instead of 1 (not 0, which fails the `xLyaHeavy > 0`
assertion of the same file): copy `~/CLOUDY/c25.00/source` without objects,
change that one line, `make`; the build was in the session scratch area and
the original tree is unchanged.

MoCHII L6 runs (grid `pn150_L6.fits`):
- `l6_metals.in`, `l6_nometals.in` (`run_l6.sh`, binary `MoCHII_l6.x`, the
  tree of 2026-10-13): 60 iterations from a fresh start, not converged.
- `l6c_metals.in`, `l6c_nometals.in` (`run_l6c.sh`): their continuations,
  converged after 8 and 20 iterations.
- `l6hc_*.in` (`run_l6hc.sh`): the collisional H I Balmer component on the
  saved states (one pass, state held).
- `l6vf_metals.in` (`run_l6vf.sh`, binary `MoCHII_vf.x`): Cloudy's
  photons as channel 8 (`par%metal_recombination_photons =
  'valence_fractions'`), from `l6c_metals`, converged after 45 iterations.
- `ledger_*.in` (`run_ledger.sh`, binary `MoCHII_ledger.x`): one pass with
  `par%absorption_ledger` on the three converged states.
- `fixed_cloudy.in`: one pass with the state held on `cloudy_state_rates.h5`
  (the fixed-structure test; channel 8 and the ledger on).
- `vf_check/`: the 2-iteration checks of the two binaries against
  `MoCHII_l6.x` (bit identity with the new options off) and of channel 8.
