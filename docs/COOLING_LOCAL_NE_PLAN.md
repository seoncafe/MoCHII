
# Plan: line cooling from level populations at the local electron density

Status: **implemented and the default since `4f9ddd3` (2026-07-26)**
(`par%cooling_model='local_ne'`; `'low_density'` keeps the `n_e -> 0` fits).
Outcome and the grading of section 8's predictions are in section 12 below.
Owner document for the last structural difference between MoCHII's thermal
balance and the other photoionization codes (MOCASSIN, TORUS, CMacIonize,
Wood's `ionize`), identified in `docs/MoCHII_cooling_analysis.pdf`.

## 1. Goal

The Tier-1 collisional line cooling is currently the `n_e -> 0` limit of
statistical equilibrium,

```
Lambda_ion(T) = T^(-1/2) sum_i A_i exp(-T_i / T)        (data/atomic/cooling_<ion>.txt)
```

with **no density suppression**: every excitation radiates, so lines whose
critical density is at or below the nebular electron density are overcooled.
Measured at fixed state against Cloudy c25.00 (`save cooling each`), the
total is 1.116x Cloudy with the low-density limit and 1.011x when the same
transitions are evaluated from the n-level solve at the local `n_e`
(`docs/MoCHII_cooling_analysis.pdf`, section on the remaining limitation).
The excess sits in the low-`n_crit` lines: [C II] 158 um, the [O III] 52/88 um
pair, and the [O II] / [S II] optical doublets.

After this work the cooling of every transition covered by an n-level model
comes from the level populations at the leaf's own `(T_e, n_e)` — the same
populations that produce the line diagnostics — so the thermal balance and
the emergent line list are one physics, not two.

Secondary consequences that make this a prerequisite rather than a polish
item:

- planetary nebulae and any `n_e >~ 10^3` model (the `HARD_SPECTRUM_PLAN.md`
  PN track) are dominated by suppressed lines; the low-density limit is not
  usable there;
- both Lexington benchmarks currently sit at the **cool** edge of the
  published spread (HII40 7894 K in 7720-8199; HII20 6376 K just under the
  6402 floor), consistent with a cooling excess of roughly ten percent.

## 2. Where the cooling is evaluated today (facts)

| site                 | file:line                       | what it does                                                                                                                              |
| -------------------- | ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| metal line cooling   | `src/species_mod.f90:588-604` | `sum over ions ne * n_el * frac(i) * cooling_eval(fit_i, T)`                                                                            |
| H I line cooling     | `src/cooling_mod.f90:257`     | `ne * nHI * cooling_eval(cool_HI_fit, T)`                                                                                               |
| fit evaluation       | `src/cooling_mod.f90:110-119` | `sum A_i exp(-T_i/T) / sqrt(T)`                                                                                                         |
| thermal bisection    | `src/thermal_mod.f90:161-184` | up to 60`net_rate` calls for each leaf (~23 in practice: 2 bracket + ~20 bisection + 1 final), ionization re-solved at each trial `T` |
| n-level solve        | `src/nlevel_mod.f90:161-225`  | `nlevel_emissivities(atom, T, ne, emis)`, `emis(k)` [erg/s per ion]; loaded and called only at output time by `lines_mod`           |
| H-impact FIR cooling | `src/species_mod.f90:470-495` | two-level low-density [C II] 158 / [O I] 63 excited by neutral H (PDR package)                                                            |

Data coverage: 28 Tier-1 cooling fits, 26 distinct n-level models (plus the
`fe_2_full`/`fe_3_full` variants selected by `par%fe_levels_full`).  **C I and
N I have a cooling fit but no n-level model** — the only two gaps.

## 3. Design decisions

### D1. A precomputed `(T, n_e)` table, not a live solve

Cost of a live solve inside the bisection: 23 `net_rate` calls x 15 active
ions x 262k leaves x 40 gas iterations = 3.6e9 evaluations of `Lambda`.  A
Gaussian elimination for a 5-40 level atom is ~2e4-5e4 flops plus the
Chebyshev `Upsilon` evaluations, i.e. ~50 us: 5e5 core-seconds.  Impossible.
A bilinear interpolation is ~10 flops, i.e. the same order as the present
`cooling_eval`.

The physical fact that makes tabulation exact rather than approximate: for a
given ion, the level populations normalized to the ion density depend **only**
on `(T_e, n_e)`.  Nothing else in the leaf enters (electron collisions +
spontaneous decay, optically thin, no photo-excitation).  So a table in
`(log T, log n_e)` carries the full physics and the only error is
interpolation error, which the gate measures.

The other codes solve live because they do it once per iteration for a
handful of ions, outside a temperature bisection.  MoCHII's bisection makes
tabulation the only route — and it is the better one, since it is exact and
covers all 28 ions.

### D2. Suppression ratio, with the transitions outside the n-level model kept in their low-density limit

Store two precomputed quantities for each ion:

```
Ssup(T, ne)  = Lambda_SE(T, ne) / Lambda_SE(T, ne_floor)      (2-D, in (0,1])
f_cov(T)     = Lambda_SE(T, ne_floor) / Lambda_fit(T)          (1-D, clamped to [0,1])
```

where `Lambda_SE(T, ne) = sum_k emis_k(T, ne) / ne` from
`nlevel_emissivities`, and `Lambda_fit` is the existing Tier-1 fit.  The
runtime cooling coefficient is

```
Lambda(T, ne) = Lambda_fit(T) * [ f_cov(T) * Ssup(T, ne) + (1 - f_cov(T)) ]
```

Why this form and not a plain `Lambda(T, ne)` table:

1. **Exact reduction to the current default.**  `Ssup -> 1` as `ne -> 0`, by
   construction (the reference column is the table's own floor), so
   `Lambda -> Lambda_fit` exactly.  The low-density limit is unchanged, not
   re-derived, so the recorded fit quality and provenance still hold.
2. **Truncated n-level models cannot bias the total.**  The Tier-1 fits come
   from the full CHIANTI `.scups` set
   (`tools/fitting/chianti_cooling.py:232`, `cooling_low_density_limit`),
   while the n-level files are truncated models (Fe III: 14 of 34 levels).
   `f_cov` is exactly the fraction of the low-density cooling the n-level
   model accounts for; the remainder `1 - f_cov` is left unsuppressed.  That
   is the correct limit for it: the missing transitions are the
   high-excitation resonance lines, whose critical densities are orders of
   magnitude above nebular values.  A plain `Lambda_SE` table would silently
   drop them.
3. **Ions with no n-level model need no special case.**  `f_cov = 0` gives
   `Lambda = Lambda_fit`, i.e. C I and N I keep working unchanged (with a
   rank-0 note listing them).
4. **`Ssup` interpolates far better than `Lambda`.**  `Lambda` spans decades
   through `exp(-T_i/T)`; `Ssup` is a smooth O(1) sigmoid in `log ne` and
   varies slowly in `T`.  A coarse grid suffices, and positivity is automatic.

`f_cov > 1` is possible where the fit has error or the `.scups` sets differ;
clamp to 1 and report the worst offender for each ion at setup (it is a
useful data-quality number, printed once).

### D3. Built at startup from the n-level files, not shipped as a new data file

The project rule is that atomic **data** are fitted offline and only evaluated
online.  A `(T, ne)` table is not new data: it is a numerical evaluation of
`data/atomic/nlevel_<ion>.txt`, which is already the offline-fitted product.
Building it in Fortran at setup buys the property that matters physically —
the cooling and the line diagnostics come from the *same* solver on the *same*
data, so they cannot drift apart, and `par%fe_levels_full` changes both
consistently.  Python remains the independent reference for the gate
(the precedent is `tools/fitting/verify_nlevel_pyneb.py`).

### D4. Grid

| axis         | range    | spacing  | points |
| ------------ | -------- | -------- | ------ |
| `log10 T`  | 2 to 6   | 0.05 dex | 81     |
| `log10 ne` | -6 to 10 | 0.2 dex  | 81     |

`ne_floor = 10^-6` is far below every critical density (the lowest is
[C II] 158 um, `n_crit,e ~ 50`), so the reference column is the low-density
limit to ~1e-8.  `ne_max = 10^10` puts the clamp beyond any nebular
application, so edge clamping is harmless; below the floor `Ssup = 1`.

Cost: 81 x 81 = 6561 solves for each ion, ~1e5 solves for a full 28-ion
registry.  With the `Upsilon` values cached for each `T` (see C1) that is a
few seconds on one core; measured in gate G-h.  Memory: 81 x 81 x 8 B = 52 kB
for each ion, 1.5 MB for the whole registry — plain per-rank arrays, no
shared window needed.

The `T` range exceeds the Tier-1 fit range (1e3-1e5) at both ends and the
`Upsilon` fits are extrapolated (clipped) outside their own range by
`ups_eval`; that is already true of the current code and is noted, not fixed
here.

### D5. Switch

`par%cooling_model = 'low_density'` (initial default) | `'local_ne'`.

Character switch to match `par%recomb_model`, `par%ci_model`,
`par%line_case`.  The default flips to `'local_ne'` in C3 once the gates
pass, because that is the physically correct setting; `'low_density'` stays
available so the recorded gates and the published paper numbers remain
reproducible.

### D6. Alternatives considered and rejected

- **Per-line `1/(1 + ne/n_crit)`.**  Cheap and needs no solve, but wrong for
  multi-level ions: it ignores cascades and level-to-level coupling, and
  `n_crit` is itself temperature dependent.  Rejected — we have the exact
  solve.
- **Offline Python `Lambda(T, ne)` files.**  Would use the full CHIANTI model
  (slightly better than the truncated one) but guarantees nothing about
  consistency with the Tier-2 diagnostics and adds 28 files to maintain.  The
  `f_cov` construction of D2 recovers most of the accuracy difference.
- **Freezing populations from the previous gas iteration.**  Cheap, but the
  bisection sweeps `T` over decades; frozen populations would make the
  bisection inconsistent with its own root.

## 4. Stages

### C0 — scaffolding, no physics change

1. New module `src/nlevel_cooling_mod.f90`: table type, registry (`add_ion`
   returning an index, `suppression(idx, T, ne)`), setup and reporting.  Empty
   of physics at this stage; `cooling_model = 'low_density'` returns
   `Lambda_fit` unchanged.
2. `par%cooling_model` in `src/define.f90` with the comment block; validation
   in `src/setup.f90` (abort on an unknown value, joining the existing
   `distance_unit` / `grid_type` / `car_walk` / `ci_model` / `line_case`
   checks).
3. **Makefile order** (a known trap in this tree): `nlevel_mod.o` is currently
   at line 104, *after* `cooling_mod.o` (91) and `species_mod.o` (93).  It uses
   only `define`, so move it before them, then insert
   `nlevel_cooling_mod.o` between `nlevel_mod.o` and `cooling_mod.o`.  No
   circular use: the new module needs nothing from `cooling_mod`.

Acceptance: build clean; a Strömgren run bit-identical to the pre-change
binary (all rates datasets, `-fp-model precise`, `np=1`).

### C1 — build the table

1. `nlevel_mod`: add a public `nlevel_cooling(atom, T, ne)` returning
   `sum_k emis_k` so the sum is not duplicated, and an **optional**
   `Upsilon` cache argument to `nlevel_emissivities` (filled on the first call
   for a given `T`, reused for the whole `ne` column).  The cache must not
   change the arithmetic: `ups_eval` is a pure function of `(atom, k, T)`, so
   cached values are bit-identical.  Requirement: `lines_mod` output
   bit-identical after the refactor.
2. `nlevel_cooling_mod`: for each registry ion that has both a cooling fit
   and an n-level file, load the atom, fill `Ssup(T, ne)` and `f_cov(T)`,
   free the atom.  Rank-0 report: ions tabulated, ions left in the
   low-density limit (expect C I, N I), `min f_cov` and `max f_cov` for each
   ion, and the build time.
3. Bilinear interpolation in `(log T, log ne)` with edge clamping; linear in
   `log T` for `f_cov`.
4. Optional (C1b, a data operation): generate `nlevel_c_1.txt` and
   `nlevel_n_1.txt` with the existing pipeline so the two gaps close.  Their
   FIR lines ([C I] 370/609 um, `n_crit ~ 10^2-10^3`) do saturate in dense
   PDR gas, so this matters for the PDR track, not for the H II benchmarks.

Gates:

- **G-a (interpolation)**: `Lambda(T, ne)` from the table vs a direct
  `nlevel_emissivities` solve at 200 random `(T, ne)` points for each ion,
  inside 1e3-1e5 K and 1e-2-1e7 cm^-3.  Criterion: worst 0.5%, median 0.1%.
  Refine the `T` axis if it fails.
- **G-b (coverage)**: table of `f_cov(T)` at 5e3 / 1e4 / 2e4 K for every ion.
  Expect `f_cov -> 1` for the dominant coolants (O III, O II, N II, S II,
  S III, Ne III).  This is a *reported diagnostic*, not a pass/fail — a low
  `f_cov` says the n-level model omits cooling transitions, and the number
  goes into the physics document.
- **G-c (low-density identity)**: with `cooling_model = 'local_ne'` and `ne`
  forced to `ne_floor`, `Lambda` equals `cooling_eval(fit, T)` to
  interpolation accuracy (this validates the D2 algebra, not the physics).
- **G-h (cost)**: setup wall time added, and per-iteration cost of a
  Strömgren/HII40 run vs `'low_density'`.  Criterion: setup < 10 s at 28
  ions, transport-inclusive iteration cost within noise.

### C2 — wire it into the thermal balance

1. `species_mod::metal_cooling` (`species_mod.f90:601`): replace
   `cooling_eval(elems(ie)%cool(i), T)` with the D2 expression, keyed by the
   stored table index.  The `ne` argument is already in the signature.
2. `cooling_mod::cooling_total` (`cooling_mod.f90:257`): same for H I.
   `nlevel_h_1.txt` exists (25 levels).  Expect `Ssup ~ 1` (Ly-alpha
   `A = 6e8`), so this is a consistency change with a measurable-but-tiny
   effect; the caveat that the H I atom is the raw optically-thin (case-A)
   model already stands in `lines_mod` and applies here too.
3. `species_mod::species_setup`: acquire the table index for each ion with
   `has_cool`.

Gates:

- **G-d (off-path bit identity)**: `cooling_model = 'low_density'` reproduces
  the pre-C2 binary exactly (all datasets, `np=1`, `-fp-model precise`).
- **G-e (cooling = diagnostics)**: on one converged model, for each ion,
  `sum over leaves ne * n_ion * Lambda_table * V` vs the sum of that ion's
  line luminosities from `lines_mod` (`<base>_lines.txt`).  These are the same
  populations, so they must agree to interpolation accuracy plus the
  `1 - f_cov` remainder.  **This is the gate that has no counterpart today**:
  the current cooling and the current line list disagree by construction.
- **G-f (fixed-state vs Cloudy)**: repeat the fixed-state cooling-budget
  comparison of `docs/MoCHII_cooling_analysis.pdf` (Cloudy c25.00
  `save cooling each`, feedback free).  Prediction: total 1.116x -> ~1.01x,
  with the change concentrated in [C II] 158, [O III] 52/88, [O II] 3727,
  [S II] 6716/6731.
- **G-g (benchmarks)**: re-run `tests/g2_hii/hii40_bench.in` and
  `hii20_bench.in`.  Predictions in section 5.

### C3 — flip the default and re-baseline

Only after C2's gates pass.  `par%cooling_model = 'local_ne'` becomes the
default; then re-baseline everything the cooling touches:

- `tests/pdr/` references (already flagged in `CLAUDE.md` as needing a refresh
  after the 2026-07-25 cooling fix — do both re-baselines in one pass);
- `tests/g4_tng/` mean `T_e`, `tests/d_dusty/`, `examples/` figures that show
  temperatures (`isrf_cloud`, `pdr`, `hei_metastable` uses `te_fixed` and is
  unaffected);
- `paper/tables/make_wood_tables.py` output (both benchmark tables) and the
  radial-structure figures in `paper/python/`;
- the recorded numbers in `CLAUDE.md`.

### C4 — H-collider saturation (separate, PDR only)

`metal_cooling_H` (`species_mod.f90:470`) is a two-level **low-density**
formula, so it has the same defect for its own lines: [C II] 158 has
`n_crit(H) ~ 3e3 cm^-3` and the `isrf_cloud` example runs at `nH = 5e3`.  The
electron term will be density-resolved after C2 while the H term is not.

Fix without a third table axis, exact for a genuine two-level system
([C II] 158 is one; [O I] 63 is the strongest transition of a three-level
`3P` system):

```
Lambda_H = Lambda_2lev(ne q^e + nHI q^H) - Lambda_2lev(ne q^e)
Lambda_2lev(C) = A dE (gu/gl) C_lu / (A + C_ul + C_lu (gu/gl))     [schematic]
```

i.e. add the *increment* the H collider makes to the saturating two-level
population, so the electron part is not double counted against the SE table.
Reduces to the present expression when both colliders are far below `n_crit`.

Gate: the PDR shell of `tests/pdr/pdr_L5.in` (unchanged, `n_e ~ 2e-2`, both
colliders unsaturated) plus a dense variant at `nH = 5e3` where the H term
must saturate.

Out of scope, noted: proton-impact excitation (CHIANTI `.psplups`) is used by
neither the present code nor this plan.

### C5 — documentation

`docs/MoCHII_physics.tex` (cooling section: the D2 formula, the grid, `f_cov`
per ion from G-b, the new self-consistency gate G-e), `docs/MoCHII_fitting.tex`
(the fits are now the low-density *carrier*, not the cooling itself),
`docs/MoCHII_cooling_analysis.tex` (the "remaining limitation" section becomes
the resolution, with the G-f numbers, and the five-code comparison loses its
one MoCHII deficiency), `data/atomic/README.md` (the density caveat paragraph
is replaced by the runtime treatment), `docs/MoCHII_UserGuide.tex`
(`cooling_model` row), `README.md` features bullet, `CLAUDE.md`, and
`paper/mochii.tex` (the limitation paragraph, the benchmark numbers, the
cooling-verification section).

## 5. Predictions (so the gates can falsify them)

At the benchmark densities (`nH = 100`, `ne ~ 110`) the suppressed lines are
the FIR fine-structure pair and the optical doublets, which the fixed-state
comparison puts at ~10% of the total cooling.  With
`d ln Lambda / d ln T ~ T_i / T ~ 3.6` for the dominant [O III] optical lines
at 8 kK, a 10% cooling reduction raises `T_e` by roughly 2.8%:

| quantity                   | now      | predicted                             | published range              |
| -------------------------- | -------- | ------------------------------------- | ---------------------------- |
| HII40`<T[NpNe]>`         | 7894 K   | ~8100 K                               | 7720-8199 (Cloudy c25: 8210) |
| HII20`<T[NpNe]>`         | 6376 K   | ~6550 K                               | 6402-6749 (Cloudy c25: 6910) |
| fixed-state total / Cloudy | 1.116    | ~1.01                                 | -                            |
| `L(Hbeta)` HII40         | 1.946e37 | slightly lower (alpha_B falls with T) | 2.01-2.10e37                 |

If `T_e` moves the wrong way, or by much more than ~4%, the D2 wiring is
wrong (most likely `f_cov` mis-scaled) — that is the first thing to check.

The PN case is where the change is large, not subtle: at `ne ~ 1e4` the
optical doublets are suppressed by factors of a few.  No prediction is
recorded here because there is no PN gate yet
(`docs/HARD_SPECTRUM_PLAN.md`).

## 6. Risks and traps

1. **Truncated n-level models.**  Handled by `f_cov` (D2), but the reported
   values must be inspected: an ion with `f_cov ~ 0.3` at 1e4 K is telling us
   its n-level file is missing cooling transitions, and the honest response is
   to extend that file, not to accept the suppression of only a third of its
   cooling.
2. **`Upsilon` cache refactor changing `lines_mod` output.**  Requirement is
   bit-identity; verify before anything else in C1 is trusted.
3. **Makefile order** (section C0.3) — this tree has been bitten twice by
   dependency order (`species_mod` before `gas_opacity_mod`, `ion_band_mod`
   before `jtally_mod`).
4. **Interpolation at the `T` edges.**  `Ssup` is evaluated over the full
   bisection bracket `[te_min, te_max]`, including the `te_min = 100 K` PDR
   setting where the `Upsilon` fits are extrapolated.  The table must be
   built over at least that range and clamped, never extrapolated in `ne`
   below the floor (return 1) or above the ceiling (return the last column).
5. **`n_e` self-consistency.**  The `ne` passed to `metal_cooling` comes from
   the converged `n_e` fixed point at the trial `T`, including metal
   electrons when `par%metal_ne`.  The suppression must use that same `ne`,
   not `nH`-derived shortcuts.
6. **Paper ripple.**  Both benchmark tables and several figures change.  The
   published PDF is committed; regenerate the tables and figures in the same
   change so the manuscript never mixes pre- and post-fix numbers.

## 7. Delegation

Advisor keeps: D1-D6 (settled above), gate criteria, reading the diffs, and
running G-a/G-e/G-f/G-g personally before any default flip.

Workers (Opus, brief must start with "read `~/.claude/CLAUDE.md` and the
project `CLAUDE.md` and follow every rule", and must carry the forbidden-word
list and the physics-based naming rule):

- **W1 = C0 + C1** (module, table, `Upsilon` cache, G-a/G-b/G-c/G-h).  Include
  in the brief: the exact D2 formulas, the grid of D4, the Makefile reorder,
  and the bit-identity requirement on `lines_mod`.
- **W2 = C2** (the two call sites + `species_setup` index, G-d/G-e).  Depends
  on W1.
- **W3 = C4** independently of W1/W2 (it touches only `metal_cooling_H`), but
  gate it after C3 so the PDR re-baseline happens once.

C3 (the default flip and the re-baseline) and C5 (documents) stay with the
Advisor: they are judgment, not implementation.

## 12. Outcome (2026-07-30)

### 12.1 What the change delivered

`4f9ddd3` made `local_ne` the default. On the fixed Cloudy state the line total
went from 1.116x Cloudy's to 1.011x, and the converged benchmarks moved as
follows.

| quantity | before | predicted | delivered | published range |
|---|---|---|---|---|
| HII40 `<T[NpNe]>` | 7894 K | ~8100 K (+2.6%) | **8210 K** (+4.0%) | 7720-8199 (Cloudy c25: 8210) |
| HII20 `<T[NpNe]>` | 6376 K | ~6550 K (+2.7%) | **7024 K** (+10.2%) | 6402-6749 (Cloudy c25: 6910) |
| fixed-state total / Cloudy | 1.116 | ~1.01 | **1.011** | --- |
| `L(Hbeta)` HII40 | 1.946e37 | slightly lower | 1.9365e37 (-0.5%) | 2.01-2.10e37 |

What followed carried these to the current **8162 K** and **6917 K**, and none of
it was a cooling change --- every later step is on the radiation-field side:

| step | HII40 | HII20 |
|---|---|---|
| zero-density fits (this memo's starting point) | 7894 | 6376 |
| `local_ne` default (this plan) | 8210 | 7024 |
| C I/N I n-level models | 8211 | 7025 |
| continuous stellar photon energies | 8194 | 6916 |
| Milne ground continua | 8202 | 6922 |
| resolved He I diffuse cascade | **8162** | **6917** |

Note that for HII20 the single largest correction after this plan was the switch
to continuous stellar photon energies, which took 7025 K back down to 6916 K.
The `local_ne` gain and that loss are independent, and it is their sum, not this
plan alone, that sets where HII20 finally sits.

### 12.2 The predictions, graded

The fixed-state prediction was excellent: 1.011 against ~1.01. The converged
temperatures were not.

- **HII40: right.** +4.0% delivered against +2.6% predicted, and the result lands
  inside the published range and within 0.6% of Cloudy c25.00.
- **HII20: the magnitude was underestimated by a factor of four.** +10.2%
  delivered against +2.7% predicted. The sign was right and the endpoint is
  within 0.10% of Cloudy c25.00 --- but section 8 set an explicit acceptance
  guard, *"If T_e moves the wrong way, or by much more than ~4%, the D2 wiring is
  [suspect]"*, and 10.2% tripped it. **That flag was never followed up.** It is
  recorded here as an open item rather than retro-fitted into a pass.

  Why the two benchmarks diverge so much is the substance of the flag: the same
  fractional cooling reduction moved HII20 two and a half times further than
  HII40. HII20 is the cooler, softer-spectrum model, where the metal fine-structure
  channels carry a larger share of the cooling and the density suppression
  therefore bites harder. That is a plausible reading, not a verified one.

- **Consequence for the benchmark verdict.** HII20 now sits **above** the nine
  published Lexington codes (6402-6749 K). So does Cloudy c25.00 (6910 K), and
  MoCHII agrees with it to 0.10%, so the modern references agree with each other
  and both exceed the 1995-era range. Documentation that still reports HII20 as
  "0.4% below the published floor" is quoting the pre-`local_ne` value and
  inverting the verdict.
