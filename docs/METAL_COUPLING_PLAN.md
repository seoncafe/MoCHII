# Metal Coupling Plan: closing the gap to the reference codes

Kwang-il Seon — drafted 2026-07-30 (Claude session; reviewed plan, not yet
started).  Companion to `docs/PLAN.md`; same staged-gate discipline.

Owner document for the two metal feedbacks that `paper/mochii.tex` records as
genuinely absent from MoCHII and present in Cloudy — two-way charge exchange
and metal recombination continua as a diffuse source — and for the
measurements that size them.  The manuscript already points here by name.

The document exists because a source survey of five reference codes
(MOCASSIN, CMacIonize, TORUS, the Monte Carlo code of Wood, Mathis &
Ercolano 2004, and Cloudy c25.00) replaced argument with `file:line` fact
on what those codes actually do.  Every citation below was read in the
distributed source.  Trees:

| Tree | Path |
|---|---|
| MOCASSIN 2.02 | `~/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/source/` |
| CMacIonize | `~/RT_Codes/CMacIonize/CMacIonize-master/src/` |
| TORUS | `~/RT_Codes/TORUS/torus/src/` |
| Wood `ionize` | `~/RT_Codes/Wood_Codes/ionize/` |
| Cloudy c25.00 | `~/CLOUDY/c25.00/source/` |

## 1. What the survey established

**No code solves all elements as one coupled system.**  All five advance one
element at a time inside nested outer iterations, and four of the five use
the same nearest-neighbor ratio chain MoCHII uses.  Cloudy differs on
exactly two axes: an N x N matrix for each element (forced by Auger, which
jumps several stages at once, so a nearest-neighbor chain cannot represent
it), and charge transfer that is two-way in the bookkeeping.  The electron
temperature is an outer root find in all five.

Structure, with citations:

- **MOCASSIN** `ionBalance` (`update_mod.f90:1321`): ratio
  `ionRatio(elem,ion)` declared `:1356`, formed `:1495-1499`, cumulative
  products `:1506-1517`, one scalar `denominator(elem)` `:1348-1349`,
  `:1520-1527`, fractions `:1537-1544`.  Outer recursion `nIterateX`,
  `maxIterateX = 25` (`:84-85`, self-call `:1608`).  The only LU in the tree
  is `luSlv`/`lured` (`emission_mod.f90:2126-2164`), for bound-level
  populations, not ionization.
- **CMacIonize** `compute_ionization_states_metals`
  (`IonizationStateCalculator.cpp:323-501`), scheme documented `:280-310`;
  carbon `:392-403` is the hand-written chain
  (`C21 = jCp1/(ne*alphaC[0])`, `C31 = C32*C21`, `sumC_inv`).
- **TORUS** `solveIonizationBalance`: blocked by nuclear charge
  (`photoionAMR_mod.F90:6859-6867`), ratio `:6880-6881`, chain `:6891`,
  renormalization `:6902-6904`, under-relaxation 0.6 `:6853`, `:6914`.
- **Wood** chain in one pass with no iteration on the metal fractions:
  carbon `ionize2.f:68-77`, nitrogen/sulfur `:96-98`; each metal block reads
  only `ne`, `h0`, `he0` from the already-converged H/He solution.
- **Cloudy** `xmat(ion_range*ion_range)` (`ion_solver.cpp:70-72`), assembled
  in `fill_array()` `:538-748`, solved `:233` through
  `newton_step.cpp:346-400` with `nrefine = 3` refinement rounds.
  Deliberately not tridiagonal: `:727-736` loops
  `for(ion_to=ion_from+1; ion_to<=IonHigh; ...)` filling arbitrary
  super-diagonal entries because Auger ejects several electrons at once
  (`:661-668`).  Elements are still separate: `ConvBase` loops
  `for(nelem=ipHYDROGEN; nelem<LIMELM; ++nelem) ... ion_wrapper(nelem)`
  (`conv_base.cpp:557-566`).  The legacy chain solver `solveions()`
  (`:1317-1364`) is dead code in c25.

**MoCHII's position.**  Nearest-neighbor ratio chain for each element in
`species_fractions` (`src/species_mod.f90:819-849`), taking `T, ne, nHI,
nHII` as inputs.  So the chain is the standard method of the field, not a
MoCHII approximation, and the real differences lie in which feedbacks are
wired in.

## 2. The founding decision this plan engages with

`docs/PLAN.md:236-241` justifies the trace-species design:

> **Trace-species approximation enables incrementality.**  At HII abundances
> metals barely perturb `n_e` or the H/He opacity, so each metal element can
> be solved *after* the H/He+thermal iteration state of each cell.

and `docs/PLAN.md:315` lists as a non-goal:

> Metals feeding back on `n_e` / opacity (trace approximation) — revisit only
> if super-solar or dust-free regimes demand it.

This plan splits that decision in two and keeps one half.

**Kept: the solver structure.**  One element at a time, solved after the
H/He and thermal state, adjacent-stage ratios chained and normalized, adding
an ion a data operation.  The survey is direct confirmation, not a
concession: four of five reference codes do exactly this, and the fifth
departs only where Auger forces it.  Nothing in this plan proposes a joint
system over elements (Section 7).

**Retired: the quantitative premise and the feedback non-goal.**  "Metals
barely perturb `n_e` or the H/He opacity" is falsified by measurement.
Removing the metal photoionization opacity moves `r(He0 = 0.5)` inward by
**0.615 pc** at HII40 (`docs/HII40_CLOUDY_RESIDUAL_INVESTIGATION.md`,
"Metal opacity"), and `par%ion_metal_abs` has been `.true.` by default since
before the current benchmarks.  The non-goal at `docs/PLAN.md:315` is
therefore already overtaken by the code.  The module header of
`src/species_mod.f90` used to read

> Trace approximation: metals never feed back on `n_e` or on the ionizing
> opacity.

which contradicted the shipped default on both counts.  It was corrected when
this plan was drafted (`src/species_mod.f90:11-20`) and now records what is
actually one-way, namely charge exchange, and points here.  The parallel
comments at the `ion_metal_abs` declaration in `src/define.f90` and elsewhere
in `src/species_mod.f90` had already been replaced by the measured 0.615 pc
statement; that module header was the one site missed.  The remaining
documentation obligation is `docs/PLAN.md:315`, still listing the feedback as
a non-goal, which Stage M5 covers.

**Priority order.** M1 first, because it is the only item in this document
that is physically wrong today rather than merely approximate.  M2 next,
because it costs runs and no code changes.  M3 is a documentation obligation
with an implementation option.  M4 is forced but belongs to another plan.

## 3. Stage M1 — two-way charge exchange

### 3.1 Correctness verdict

**Physically wrong today.**  The reaction

```
X^{i+} + H^+  ->  X^{(i+1)+} + H^0        (metal ionized, hydrogen neutralized)
X^{(i+1)+} + H^0  ->  X^{i+} + H^+        (metal recombined, hydrogen ionized)
```

moves one unit of charge between hydrogen and a metal.  MoCHII records only
the metal's side.  `species_fractions` (`src/species_mod.f90:832-836`) adds

```
rion = rion + nHII*cx_rate(cxi_form(it), cxi(1:6,it), T)
rrec = rrec + nHI *cx_rate(cxr_form(it), cxr(1:6,it), T)
```

and the cached path does the same (`:672`, `:674`), while `solve_ion_cell`
(`src/ion_balance_mod.f90:42-96`) has no charge-exchange term at all — a
grep for `cx_` in `src/ion_balance_mod.f90` and `src/thermal_mod.f90`
returns nothing.  So hydrogen is neutralized or ionized by these reactions
without its balance knowing, and **charge is not conserved in the reacting
pair**.  That is a defect regardless of its size, and the fix is justified
regardless of what the measurement in 3.5 returns.

It is also the one structural axis on which Cloudy differs from all four
Monte Carlo codes.  Two-way, code by code:

| code | two-way? | evidence |
|---|---|---|
| Cloudy c25 | **yes** | `iso_charge_transfer_update()` `iso_ionize_recombine.cpp:18-50`, called at the top of every ion solve (`ion_solver.cpp:68`), accumulates hydrogen's CT rates over all heavier elements and all their stages `:38-45`; those totals enter hydrogen's own matrix rows at `ion_solver.cpp:687` (`ction`, off-diagonal pair `:696-697`) and `:709` (`ctrec`, `:723-724`), plus the H-like level matrix `iso_level.cpp:121-124`.  Both sides read the same `CharExcIonOf`/`CharExcRecTo` arrays, and `ChargTranEval()` runs once per outer loop precisely so that they do (`conv_base.cpp:498-502`). |
| MOCASSIN | no | `chex(nElements,nstages,4)` (`update_mod.f90:1358-1359`, populated `:1392-1451`) enters the metal denominator `:1498-1499`; `chex(1,...)` is never assigned, so hydrogen's equation reduces to the generic line at `elem=1, ion=1`.  Reverse rate `revRate` `:1477-1480` is nonzero only for N I and O I (`deltaE_k` zeroed `:1393`, two entries set `:1448-1449`, gated `:1475`). |
| CMacIonize | no | metal CT only in the metal expressions, e.g. oxygen `IonizationStateCalculator.cpp:438-443`; `grep charge_transfer` returns no hit above line 501, so the H/He solve `:649-753` has none. |
| TORUS | no | `getChargeExchangeRecomb` has N IV, O II, O III only (`photoion_utils_mod.F90:840-868`), `getChargeExchangeIon` N I, O I only (`:885-906`); the H balance calls the same routines with H II / H I and gets 0 (`photoionAMR_mod.F90:6872-6881`). |
| Wood | no | `ct1.f` (Kingdon & Ferland package) enters the metal denominators, e.g. O `ionize2.f:136-138`; `h0find.f:1` has no CT argument and its stated equation set `:6-13` has no CT term. |

### 3.2 The fix keeps hydrogen's closed form

Hydrogen's balance with the charge-exchange terms is

```
n_HI (Gamma_H + C_H n_e + R_cx_out)  =  (alpha_H n_e + R_cx_in) n_HII
```

with, summed over every registry element and every transition,

```
R_cx_in  = sum_el sum_it  k_cxi(el,it) n(X_el,it)       [s^-1]  (H+ -> H0)
R_cx_out = sum_el sum_it  k_cxr(el,it) n(X_el,it+1)     [s^-1]  (H0 -> H+)
```

Eliminating `n_HII = n_H (1 - x_HI)` gives

```
x_HI = (alpha_H n_e + R_cx_in)
       / (Gamma_H + (C_H + alpha_H) n_e + R_cx_in + R_cx_out)
```

**Both** charge-exchange terms appear in the denominator, for the same
reason `alpha_H` does: `x_HII = 1 - x_HI` was used to eliminate `n_HII`, so
every term multiplying `n_HII` moves into the denominator alongside its
appearance in the numerator.  Writing the denominator as
`Gamma_H + (C_H + alpha_H) n_e + R_cx_out` alone would be an error and would
break exactly the charge conservation the stage is fixing — the gate in 3.5
would catch it.

So the fix needs **no Newton solve and no matrix**.  It replaces
`src/ion_balance_mod.f90:76-78`

```
denH = gH + (cH + aH)*ne
xHI  = aH*ne / denH
```

with

```
denH = gH + (cH + aH)*ne + r_cx_in + r_cx_out
xHI  = (aH*ne + r_cx_in) / denH
```

Cloudy does the same physics with more work: it assembles hydrogen's totals
in `iso_ionize_recombine.cpp:38-45` and inserts them into hydrogen's own
matrix rows (`ion_solver.cpp:687`, `:709`).  MoCHII's hydrogen is one
equation, so the insertion is two added terms.

**The rate cache is already in the right place.**  `species_ne_prepare(il,
T)` (`src/species_mod.f90:636-655`) already evaluates `cch_cxi` and
`cch_cxr` at the leaf's temperature, and `species_ne_cached`
(`:662-688`) already walks the stage chain inside the same `n_e` fixed point
of `solve_ion_cell` (`src/ion_balance_mod.f90:74-95`) at the same
temperature.  The metal densities `R_cx_in` and `R_cx_out` need are exactly
the ones that walk produces.  Nothing new has to be iterated, and no new
module is needed — so no Makefile dependency reorder, the trap that has bitten
this tree three times (`docs/COOLING_LOCAL_NE_PLAN.md`, section 6.3).

`src/thermal_mod.f90:47` calls `solve_ion_cell` at every trial temperature of
the bisection, so the term enters the thermal balance with no separate
wiring.

### 3.3 Implementation

#### 3.3.1 One chain walk, not two

The stage densities `R_cx_in` and `R_cx_out` need are the ones
`species_ne_cached` (`src/species_mod.f90:662-688`) already forms: `fr` after
`fr = fr*r(i-1)` is the fraction of stage `i`, and stage 1 is `fr = 1/s`.  So
the two sums are three added lines inside the existing loop, not a second
walk.

**Do not write a separate routine that repeats the walk.**  An earlier draft
of this plan proposed exactly that, with a note to merge later if profiling
demanded it.  That is the wrong default in this tree: this session found seven
files carrying drifted copies of the same physics, which is why
`tests/rates/check_python_rates.py` exists at all, and the `n_e` chain and the
charge-exchange chain must agree by construction rather than by a gate.  The
change is therefore

```
species_ne_cached(nH, ne, nHI, nHII)                          [function]
  ->  species_metal_sums(nH, ne, nHI, nHII, nem, r_cx_in, r_cx_out)  [subroutine]
```

with `species_ne_cached` kept as a thin wrapper returning only `nem`, so the
`charge_neutrality_ne` call site (`src/ion_balance_mod.f90:125`) and the
electron-sum gates are untouched and stay bit-identical.  Gate the refactor on
bit-identity of `nem` **before** adding the new terms, as a separate commit.

Name rationale: `species_metal_sums` states the quantity; a name like
`species_ne_and_cx` would encode the caller's needs rather than the physics.

#### 3.3.2 Where the terms are evaluated, and the one-iterate lag

`solve_ion_cell` (`src/ion_balance_mod.f90:74-95`) is a damped fixed point on
`n_e` in which `x_HI` is recomputed at the top of each pass from the current
`n_e`, and `n_e` is then re-formed from the new `x_HI`.  `R_cx_in` and
`R_cx_out` depend on the metal stage fractions, which depend on `n_HI` and
`n_HII` through the charge-exchange terms of `species_fractions` itself, so
they cannot be evaluated from the `x_HI` of the pass that is being computed.

They are therefore evaluated **once per pass, from the previous iterate**, in
the same place and with the same arguments the electron sum already uses:

```
do it = 1, 200
   ne_old = ne
   call species_metal_sums(nH, ne_old, nH*xHI, nH*(1-xHI), nem, r_in, r_out)   ! <- xHI from the previous pass
   denH = gH + (cH + aH)*ne + r_in + r_out
   xHI  = (aH*ne + r_in) / denH
   ... He unchanged ...
   ne = nH*((1-xHI) + yHe*(xHeII + 2*xHeIII))
   if (with_metal_ne) ne = ne + nem
   ne = max(0.5*(ne + ne_old), 1e-12*nH)
   if (abs(ne - ne_old) <= 1e-10*ne) exit
end do
```

Two things to note.  The call moves **above** the `x_HI` update, where today's
`species_ne_cached` call sits below it; the electron sum consequently sees
`x_HI` from the previous pass rather than the current one, which is a change
of iteration order, not of the converged answer, but it means the refactor
commit of 3.3.1 is **not** bit-identical once the call moves — so move the
call in the same commit that adds the terms, and keep the bit-identity gate on
the refactor alone.  Second, `denH > 0` is guaranteed as before, since both
new terms are non-negative.

#### 3.3.3 The `metal_ne` trap

`species_ne_prepare` is called only under `with_metal_ne`
(`src/ion_balance_mod.f90:56-61`), which is
`par%metal_ne .and. par%use_metals .and. n_elements > 0 .and. present(il)`.
With `par%metal_ne = .false.` the cache — including `cch_cxi` and `cch_cxr` —
is never filled, so a naive implementation would silently read stale or zero
coefficients and produce no charge exchange while reporting success.

Charge exchange has nothing to do with whether metal electrons are counted.
The guard must therefore be split:

```
with_metal_ne = par%metal_ne .and. par%use_metals .and. n_elements > 0
with_metal_cx = two_way    .and. par%use_metals .and. n_elements > 0
if (present(il) .and. (with_metal_ne .or. with_metal_cx)) call species_ne_prepare(il, T)
```

with `nem` added to `n_e` only under `with_metal_ne` and the two
charge-exchange terms applied only under `with_metal_cx`.  Both call sites of
`solve_ion_cell` (`src/thermal_mod.f90:47`, `src/ion_balance_mod.f90:188`)
already pass `il`, so no interface change is needed; the `present(il)` test
stays as the guard it is.

#### 3.3.4 What must not change

- `species_fractions` (`src/species_mod.f90:819-849`) is correct as it stands
  and is not touched.  M1 adds hydrogen's side of a reaction the metal side
  already has; it does not re-derive the metal side.
- The metal coefficients themselves are not re-fitted.  `cch_cxi` and
  `cch_cxr` are read as they are, so both sides of every pair use one number
  by construction.
- `charge_neutrality_ne` (`src/ion_balance_mod.f90:123`) needs no new term: it
  solves `n_e` at **fixed** H/He fractions, and the fractions it is given
  already carry the charge exchange because they come from `solve_ion_cell`.
- Helium's balance is untouched, because MoCHII carries no helium charge
  exchange at all (3.7 item 1).

#### 3.3.5 Change list

| file | change |
|---|---|
| `src/species_mod.f90` | `species_ne_cached` -> `species_metal_sums` plus a thin wrapper; three added lines in the existing chain loop; export the new name |
| `src/ion_balance_mod.f90` | split the guard (3.3.3); move the call above the `x_HI` update and add the two terms (3.3.2); `hydrogen_balance_residual` for the diagnostic (3.4) |
| `src/define.f90` | `par%charge_exchange` declaration with the comment that `'metal_only'` is a reproduction value, not a physics option |
| `src/setup.f90` | validate the switch beside the existing `ci_model` / `line_case` checks; abort on anything else |
| `src/main.f90` | print the residual diagnostic beside `max |delta x_HII|` (`:158-176`) |
| `tests/charge_exchange/` | the new fixed-state gate of 3.6 |
| `docs/MoCHII_physics.tex`, `docs/MoCHII_UserGuide.tex` | the balance equation and the `par%` entry |

**Switch.**  `par%charge_exchange = 'two_way'` (new default) | `'metal_only'`
(the present bookkeeping, retained only to reproduce the recorded gates and
the published benchmark numbers).  Character switch to match
`par%recomb_model`, `par%cooling_model`, `par%ci_model`, `par%line_case`;
validated in `src/setup.f90` beside the existing `distance_unit` /
`grid_type` / `car_walk` / `ci_model` / `line_case` checks.  `'metal_only'`
is a reproduction value, not a physics option, and the declaration comment
must say so.

#### 3.3.6 Order of work

1. Refactor to `species_metal_sums` with the call left where it is; gate on
   bit-identity of `nem` against the current binary on `tests/g2_hii/smoke.in`.
2. Add the switch, the split guard, the moved call, and the two terms.
3. The fixed-state conservation gate (G-M1a) — before any benchmark run, since
   it is the one test that can fail for a sign or stage-index error.
4. The diagnostic print, then the benchmark and PDR runs (G-M1b to G-M1e).

`hydrogen_balance_residual(il, T, nH, xHI, ne)` in `ion_balance_mod` returns
the relative residual of the equation in 3.2 on the **written** state,
including both charge-exchange terms; it is the convergence diagnostic of 3.4.

**Switch.**  `par%charge_exchange = 'two_way'` (new default) | `'metal_only'`
(the present bookkeeping, retained only to reproduce the recorded gates and
the published benchmark numbers).  Character switch to match
`par%recomb_model`, `par%cooling_model`, `par%ci_model`, `par%line_case`;
validated in `src/setup.f90` beside the existing `distance_unit` /
`grid_type` / `car_walk` / `ci_model` / `line_case` checks.  `'metal_only'`
is a reproduction value, not a physics option, and the declaration comment
must say so.

### 3.4 The diagnostic Cloudy has and MoCHII would not

Cloudy checks the balance and acts on it: `conv_base.cpp:578-583` forms

```
netion[nlo][nhi] = CharExcRecTo[...]*dl1*ul2 - CharExcIonOf[...]*ul1*dl2
```

and declares the zone unconverged with the reason `"CX inconsistency"` at
`:611-620`.  It needs that check because within a single LU the partner's
density is a frozen coefficient, so exact charge conservation is reached by
the outer iteration plus the check, not inside one matrix solve.

MoCHII's situation after M1 is structurally better and the plan should say
so rather than copy the verdict.  Both sides are formed in the same leaf
solve from the same cached coefficients and the same stage densities, so the
bookkeeping mismatch is zero by construction *within* the solve.  What can
still be inconsistent is the **written** state, for two reasons already
documented in the tree: the H/He fractions are under-relaxed toward the
solved state (`par%ion_relax`, `src/ion_balance_mod.f90:193-195`) while the
metal fractions written alongside are evaluated at the written `n_e` and `T`
— which is why `charge_neutrality_ne` exists at all
(`src/ion_balance_mod.f90:109-114`).

So the diagnostic is the residual, not a rejection:

```
eps(il) = | n_HI (Gamma_H + C_H n_e + R_cx_out)
            - (alpha_H n_e + R_cx_in) n_HII |
          / [ n_HI (Gamma_H + C_H n_e + R_cx_out) ]
```

evaluated on the written state, reported each gas iteration as the maximum
over leaves and the volume-weighted mean, printed beside
`max |delta x_HII|` in `src/main.f90:158-176`.  At minimum it is printed; a
later decision may make it gate convergence, but not before the benchmarks
show what value it settles at.

### 3.5 Magnitude: measured on the frozen converged state

Measured this session from the converged HII40 production run
(`results/continuous_energy/production/hii40_continuous_rates.h5`), by forming
both sums leaf by leaf from the stored `Gamma_HI`, `n_e`, `T_e` and metal stage
fractions and recomputing `x_HI` with the closed form of 3.2.  This is a
**frozen-state** estimate: it does not re-converge, so it does not carry the
feedback `x_HI -> tau -> Gamma` or the shift in oxygen's own fractions.  It
does size the terms against hydrogen's own.

| r [pc] | `x_HI` | `R_cx_in / (alpha_A n_e)` | `R_cx_out / (Gamma_H + C_H n_e)` | `x_HI` ratio |
|---|---|---|---|---|
| 1.20 | 6.53e-05 | 0.0000 | 4.38e-05 | 0.999957 |
| 2.00 | 1.72e-04 | 0.0000 | 1.10e-04 | 0.999896 |
| 3.20 | 5.35e-04 | 0.0001 | 3.67e-04 | 0.999709 |
| 4.00 | 1.51e-03 | 0.0008 | 1.36e-03 | 0.999505 |
| 4.40 | 4.46e-03 | 0.0049 | 5.22e-03 | 0.999749 |
| 4.67 (`x_HI = 0.5`) | 0.5 | **1.1954** | **0.8205** | **1.1127** |

So the expectation is confirmed in the ionized zone and **understated at the
front**.  Through the ionized nebula the effect on `x_HI` is 0.005 to 0.05%.
At the ionization front the two charge-exchange terms are 1.20 times
hydrogen's own recombination rate and 0.82 times its own ionization rate, and
the frozen-state `x_HI` rises 11.3%.

At the front the term is **almost entirely O0/O+ against H+/H0**: 2.835e-11 and
2.838e-11 s^-1 for the two directions, over 99.99% of each.  That is the
near-resonance -- IP(O) = 13.618 eV against IP(H) = 13.598 eV, so
`dE = 0.020` eV and both directions are fast, which is what locks oxygen's
ionization to hydrogen's.  Note what the four-digit near-cancellation means
for the omission: what is missing from hydrogen's balance is not one small term
but **two comparable terms that would have largely cancelled**, so the net is a
small difference of large numbers.  That is an argument for computing it rather
than estimating it, and it is why G-M1a tests conservation to round-off rather
than to a tolerance.

Consequence for the HII40 helium-front residual (Section 8): **none**.  The
helium front sits near 3.8 pc, where the measured effect on `x_HI` is below
0.05%.  This closes by measurement what Section 8 previously argued.

Also note `abund_Mg`, `abund_Fe`, `abund_Si` are zero in the Lexington
benchmarks, so the three largest `CXI` coefficients at 8169 K -- Fe 4.00e-09,
Si 1.23e-09, Mg 1.20e-10 -- contribute nothing there.  A metal-rich or
depletion-free model would weight this term very differently, and the front
estimate above does not transfer to one.

The reaction census, taken from `data/atomic/element_*.txt` this session
(11 element files, `MAX_ST = 6` so at most five transitions each):

- **7 reactions take H+ to H0** (nonzero `CXI`): C, N, O, Mg, Si, S, Fe — all
  at the neutral stage, i.e. `X^0 + H^+ -> X^+ + H^0`.
- **22 reactions take H0 to H+** (nonzero `CXR`): C 3, N 2, O 2, Ne 1, Mg 2,
  Si 4, S 3, Ar 2, Fe 3.  Ca and Cl have none.
- Three of the 22 carry a nominal `1.0e-14 cm^3 s^-1` with no temperature
  dependence (Ne, S, Ar) — a placeholder, not a determination, and worth
  labeling as such in `data/atomic/README.md`.

So hydrogen's two sums run over 29 reactions.  Note the asymmetry: the
H0-consuming direction outnumbers the H0-producing one three to one and
reaches higher stages, which is why the term is expected to *ionize* hydrogen
on net in partially neutral gas.

### 3.6 Gates

- **G-M1a (exact pairwise charge conservation).**  **PASS, 2026-07-31**
  (`tests/charge_exchange/run_check.sh`, exit 0).  Fixed-state test on five
  converged cell states spanning `x_HI` = 0.025 to 1.0 and T = 300 to 10000 K.
  The census is confirmed from the shipped data: 7 CXI + 22 CXR = 29.  Every
  reaction agrees between hydrogen's side and the cascade's to **<= 6.50e-16**
  relative (worst CXR Fe transition 2; worst CXI 4.36e-16 on S), and the closed
  form of 3.2 is satisfied to **<= 4.96e-14**.
  **The denominator assertion fires as designed**: dropping `R_cx_in` from the
  denominator makes the residual jump to between **1.31e-3** and **1.75e+04**,
  and two states then return unphysical `x_HI` above 1.  The first draft of the
  brief for this stage carried exactly that truncated form, so this gate is the
  reason it did not ship.
- **G-M1b (reproduction path).**  **PASS, 2026-07-31.**
  `par%charge_exchange = 'metal_only'` reproduces the frozen pre-M1 binary on
  `tests/g2_hii/smoke.in` at `np = 1`: **31/31 datasets bit-identical**.
  Caveat that 3.3.2 predicts and this gate confirms: with `par%metal_ne = .true.`
  it is *not* bit-identical, because the call moved above the `x_HI` update.
  The deviations are at the fixed-point tolerance, <= 6e-11 relative on `x_HI`,
  `n_e` and `Gamma`.
- **G-M1c (benchmarks).**  **MEASURED 2026-07-31**, both benchmarks, same
  inputs and Sobol seed, only `par%charge_exchange` changed
  (`results/continuous_energy/m1_twoway/`):

  | | `<T[NpNe]>` K | `L(Hbeta)` | `r(He0=0.5)` | `r(C2=0.5)` | `r(H0=0.1)` | `r(H0=0.5)` |
  |---|---|---|---|---|---|---|
  | HII40 `metal_only` | 8168.9 | 1.9293e37 | 4.1278 | 4.0630 | 4.6111 | 4.6770 |
  | HII40 `two_way` | 8167.5 | 1.9301e37 | 4.1277 | 4.0631 | 4.6116 | 4.6775 |
  | delta | **-1.41 (-0.017%)** | +0.045% | -0.0002 | +0.0001 | +0.0005 | +0.0005 |
  | HII20 `metal_only` | 6924.5 | 4.6837e36 | 1.2408 | 1.1421 | 2.8352 | 2.8797 |
  | HII20 `two_way` | 6924.4 | 4.6840e36 | 1.2408 | 1.1421 | 2.8352 | 2.8798 |
  | delta | **-0.13 (-0.002%)** | +0.007% | -0.0000 | +0.0000 | +0.0001 | +0.0001 |

  **The recorded table values are NOT re-baselined on this, and the reason is
  quantitative rather than editorial.**  Both runs use `par%gas_tol_te = 1.0e-2`,
  so the solver's own convergence tolerance on `T_e` is 82 K at HII40 and 69 K
  at HII20.  The measured shifts are **58 times** and **533 times** smaller than
  that.  They are reproducible -- same seed, one switch -- but they are far
  inside the tolerance the runs were converged to, so carrying them into
  Tables 4 and 5 would quote a difference the calculation does not resolve.
  What the paper states instead is the bound.

  This is the expected outcome of 3.8 rather than a surprise: the two
  charge-exchange terms nearly cancel, the net is the departure from the
  oxygen-hydrogen lock, and the Lexington diagnostics are weighted by `n_p n_e`,
  which suppresses exactly the gas where the departure is largest.

  **Consequence for the helium-front residual: none, now measured.**  The He
  front moves by 0.0002 pc against the 0.17750 pc residual -- three orders of
  magnitude short.  Section 8's prediction is confirmed by direct measurement
  rather than by argument.

  The `hydrogen_balance_residual` diagnostic behaves as 3.8 predicts on these
  runs: maximum 4.6e-3, in the capped neutral cells, with the volume-weighted
  mean at 5e-4 and falling through the iteration.

- **G-M1d (the regime where it should matter most).**  `tests/pdr/pdr_L5.in`
  and `tests/pdr/pdr_dense.in`, where `x(H0)` approaches 1 and the metal
  electrons are the only electrons.  Report `x(H0)`, `n_e`, `T`, and the
  C II and O I fractions in the PDR shell.  Criterion: the run stays
  converged, the changes are reproduced by a hand evaluation of the two sums
  in one shell cell, and the PDR references are re-baselined in the same
  change.
- **G-M1e (written-state consistency).**  The `eps` of 3.4 at convergence,
  maximum and volume-weighted, for both benchmarks and both PDR cases, with
  `par%ion_relax = 1` and `< 1`.  Reported diagnostic, recorded in
  `docs/MoCHII_physics.tex`; expected to sit at the fixed-point tolerance
  with `ion_relax = 1`.

### 3.7 Two open items M1 does not close

1. **MoCHII carries no helium–hydrogen charge exchange at all.**  There is no
   `element_he.txt`, and `grep` for charge exchange in `src/` finds only the
   metal path.  MOCASSIN does carry it, one-way: `chex(2,1,:)` and
   `chex(2,2,:)` are set at `update_mod.f90:1395-1396`, though with
   `deltaE_k(2,*) = 0` it has no reverse rate.  Whether it matters is
   **UNMEASURED**, and there is no rate in the tree to evaluate it with, so
   this is a data item plus a measurement, not part of M1.  It is worth
   flagging for a separate reason: it is a **helium-specific** term, and it
   is not among the candidates the HII40 front investigation has closed
   (Section 8).  That does not make it the explanation — it makes it the one
   item the exclusion list never had a row for.
2. **Charge-exchange heating and cooling.**  The reactions are exothermic or
   endothermic and therefore exchange energy with the gas.  Whether the
   reference codes carry that term was **not read** in this survey, and
   MoCHII does not carry it.  Flagged, not planned.

### 3.8 The slow mode M1 introduces, and why it is structural

Implementing M1 exposed a convergence effect that the plan did not anticipate,
and it turns out to be a property of the physics rather than of the code.

**The measurement.**  With `two_way`, `solve_ion_cell` converges monotonically
but slowly in partially neutral gas.  At T = 1e4 K, n_H = 100 and the HII20
abundances, the cell iteration meets its exit test after 124 passes at
`x_HII = 0.18`, 239 at 0.10, and **1285 at 0.021**, against the hard 200-pass
cap in `solve_ion_cell`.  Capped, such a cell leaves the loop with `x_HII` high
by about 8% (2.252e-2 against 2.078e-2 converged) and a residual near 1e-3.
Repeated whole-cell solves converge it to <= 5e-12, and `metal_only` shows none
of it (2e-10 everywhere), so the slow mode comes from the lagged `R_cx` terms.

**Why.**  Where the near-resonant pair locks oxygen to hydrogen -- that is,
where `Gamma_O` and `n_e alpha_O` are negligible against the charge-exchange
terms --

```
f(O+)/f(O0) = (n_HII k_cxi)/(n_HI k_cxr)   =>   R_cx_out/R_cx_in = n_HII/n_HI
```

Substituting that into the closed form of 3.2 and solving gives, exactly,

```
x_HI = alpha_H n_e / (Gamma_H + (C_H + alpha_H) n_e)
```

**the charge-exchange terms cancel identically and return the no-CX solution**
(verified symbolically).  So M1's net effect is not the terms themselves but the
*departure* from that lock, which `Gamma_O` and `n_e alpha_O` set.  And the
slope of the substitution map at its fixed point is, in the same limit,

```
dF/dx = R_cx_in / (R_cx_in + alpha_H n_e)
```

which tends to 1 wherever `R_cx_in >> alpha_H n_e` -- exactly the neutral gas.
The iteration is therefore **marginally stable by construction**, and the pass
count follows `ln(tol)/ln(slope)`, which reproduces the measured 1285 to within
its order of magnitude.  This is the same near-cancellation that 3.5 flagged as
a numerical delicacy, now showing up in the solver rather than in an estimate.

**What it does and does not affect.**  Cells with `x_HII >~ 0.1` converge inside
the cap, so `r(H0 = 0.5)` and everything inside the front are unaffected.  The
Lexington diagnostics are weighted by `n_p n_e`, which is proportional to
`x_HII^2`, so a cell at `x_HII = 0.02` carries a relative weight near 4e-4:
`<T[NpNe]>` and `L(Hbeta)` are safe.  **What is at risk is the shadow beyond the
front and any PDR application**, which is precisely where M1 was expected to
matter, so this must be fixed before an M1 result there is quoted.

**The fix this analysis points at.**  Do not iterate on a quantity that mostly
cancels.  Replace the substitution step for `x_HI` with a Newton or secant step
on the scalar residual

```
R(x) = x (Gamma_H + (C_H + alpha_H) n_e + R_cx_in + R_cx_out) - (alpha_H n_e + R_cx_in)
```

which converges in a few passes whatever the slope is, since the slow mode is a
property of the substitution map and not of the root.  Raising the cap or adding
`|delta x_HI|` to the exit test are the cheap alternatives, but they treat the
symptom: the first costs 1285 passes for each of 23 temperature trials for each
leaf, the second only detects the failure.  **Not implemented; the cap and the
exit test were left as they are**, and the code site carries the measured
numbers.

**Gate G-M1f (new).**  A cell at `x_HII` between 0.01 and 0.05 must converge to
the same state as repeated whole-cell solves, to <= 1e-10 relative, inside the
pass cap.  Add it to `tests/charge_exchange/` when the fix lands.

## 4. Stage M2 — measure what the trace separation costs

### 4.1 Correction to the premise

`par%metal_ne` and `par%metal_heat` are both `.false.` by default
(`src/define.f90:580-581`), declared as "PDR-zone physics (each defaults off
so the recorded gates reproduce)" (`:568-570`).  **The published benchmark
runs do not use those defaults.**  `tests/continuous_energy/hii40_continuous.in`
and `hii20_continuous.in` set, at lines 49-52,

```
par%use_metals      = .true.
par%ion_metal_abs   = .true.
par%metal_ne        = .true.
par%metal_heat      = .true.
```

and `results/continuous_energy/production/` — the outputs
`paper/tables/make_wood_tables.py:23,33` reads — was produced from them.  So
is every run in the residual investigation
(`results/continuous_energy/metal_bound/*/in.in`,
`metal_sigma/*/in.in`).  The older grouped-mode gate inputs
`tests/g2_hii/hii{40,20}_bench.in` set neither, and
`tests/g2_hii/hii40_mfull.in` exists to set both.

**The consequence is the reverse of the obvious one:** the published physics
is the complete one, and the *default* configuration is the untested one.
Anyone building MoCHII and running it out of the box gets a different
treatment from the one the paper reports, and no number anywhere in the tree
says how different.

### 4.2 What to measure and why it is a consistency question

MOCASSIN always includes metal electrons, charge-weighted over every
switched-on element (`update_mod.f90:1564-1577`, written back `:1592-1594`),
and so does TORUS (`ion_mod.F90:1418-1421`, recomputed after every
ionization update, `photoionAMR_mod.F90:6916`, `:5941`).  CMacIonize does not
(`IonizationStateCalculator.cpp:134`, with the comment `:273-274` "we neglect
free electrons coming from ionization of coolants"), nor does Wood
(`ioneng.f:45`, stated as deliberate at `tbal1.f:39-41`).  MOCASSIN includes
metal photoheating (the loop runs `do elem = 1, nElements`,
`update_mod.f90:1128-1135`); CMacIonize does not
(`TemperatureCalculator.cpp:297`, plus the He Lyman-alpha, PAH and cosmic-ray
terms only), nor does Wood (`ioneng.f:53-54`).

So comparing MoCHII against MOCASSIN on the helium structure or on the
temperature while MoCHII has metal electrons switched off is a
code-comparison difference of known origin, and it should be quantified
rather than assumed small.  With `metal_ne` on, as the published runs have
it, that particular difference is absent — which is itself a statement the
paper can only make once the switch matrix is measured.

### 4.3 Gate

**G-M2: the switch matrix.**  Four combinations of
`(par%metal_ne, par%metal_heat)`, HII40 and HII20, everything else at the
production input, same Sobol seed and rank count as the recorded runs:

**MEASURED 2026-07-31**, on the pre-M1 binary as this section requires, from
`tests/continuous_energy/hii{40,20}_continuous.in` with only the two switches
changed (`results/continuous_energy/metal_switches/`).  Fronts from
`tests/continuous_energy/analyze_hii40_residuals.py` itself, so they are
comparable with the recorded residual measurements.

HII40:

| `metal_ne` | `metal_heat` | `<T[NpNe]>` K | `L(Hbeta)` | `r(He0=0.5)` | `r(C2=0.5)` | `r(H0=0.1)` |
|---|---|---|---|---|---|---|
| on | on | 8168.9 | 1.9293e37 | 4.1278 | 4.0630 | 4.6111 |
| on | off | 8147.9 | 1.9306e37 | 4.1248 | 4.0606 | 4.6090 |
| off | on | 8168.7 | 1.9291e37 | 4.1291 | 4.0639 | 4.6123 |
| off | off (shipped default) | 8147.6 | 1.9304e37 | 4.1259 | 4.0616 | 4.6100 |

HII20:

| `metal_ne` | `metal_heat` | `<T[NpNe]>` K | `L(Hbeta)` | `r(He0=0.5)` | `r(C2=0.5)` | `r(H0=0.1)` |
|---|---|---|---|---|---|---|
| on | on | 6924.5 | 4.6837e36 | 1.2408 | 1.1421 | 2.8352 |
| on | off | 6918.8 | 4.6839e36 | 1.2406 | 1.1421 | 2.8345 |
| off | on | 6924.4 | 4.6836e36 | 1.2408 | 1.1422 | 2.8358 |
| off | off (shipped default) | 6918.7 | 4.6838e36 | 1.2407 | 1.1420 | 2.8352 |

**The two switches separate cleanly, and they are not equally important.**

- `metal_heat` carries essentially the whole effect: **-21.1 K (0.26%)** at HII40
  and **-5.7 K (0.082%)** at HII20 when switched off.
- `metal_ne` is negligible for these diagnostics: **-0.2 K (0.002%)** and
  **-0.1 K (0.001%)**.  That is the expected size — the metal stages carry about
  1.1e-3 of `n_H` in charge against `n_e/n_H` near 1.1, so about 0.1% of the
  electrons, and `<T[NpNe]>` is weighted by `n_p n_e`, which suppresses the
  neutral gas where metal electrons actually dominate.
- The two are **additive to the digit**: at HII40, (-21.1) + (-0.2) = -21.3,
  which is the both-off value; at HII20, (-5.7) + (-0.1) = -5.8.  There is no
  interaction between them at this level.
- Fronts move by at most 0.003 pc, and `L(Hbeta)` by at most 0.07%, in the
  direction opposite to the temperature.

So the shipped default is **21.3 K (0.26%) cooler than the published HII40
number and 5.8 K (0.084%) cooler at HII20** — smaller than the published
cross-code spread (7720-8199 K) and than the 0.50% gap to Cloudy c25, but real,
one-signed, and until now unrecorded anywhere.

Note the caveat this measurement does NOT cover: `metal_ne` being negligible
here says nothing about neutral gas, where the metal stages are the only
electron source and `n_e` would collapse without them.  These two benchmarks
weight that region at essentially zero; `tests/pdr/` does not.

Acceptance is that the table exists and is recorded in
`docs/MoCHII_physics.tex` and in the manuscript, so both the paper and the
User Guide can state a number instead of an argument.  Add
`r(H0 = 0.1)` alongside for comparability with the recorded front
measurements.

**Decision the gate forces.**  Metal electrons and metal photoheating are
physically real, so switching them off is an approximation, and the house
rule is that a deliberate approximation is acceptable only when documented at
the code site **with its domain of validity**.  "Defaults off so the recorded
gates reproduce" is a reproducibility statement, not a domain of validity.
After G-M2, either flip both defaults to `.true.` (making the shipped
configuration the one the paper reports) or write the measured domain of
validity at `src/define.f90:568-583`.  Backward compatibility alone does not
settle it.  Note that flipping requires re-baselining the grouped-mode gate
inputs `tests/g2_hii/hii{40,20}_bench.in`, which would then change behavior
without changing their text.

## 5. Stage M3 — metal recombination continua as a diffuse source

### 5.1 Status: a justified approximation that is not documented

`grep metal src/diffuse_mod.f90` returns nothing.  The diffuse field has
exactly four channels (`dif_ch(4, nleaf)`, `src/diffuse_mod.f90:40`), all
H/He: the three ground continua (`:92`, `:94`, `:96`) and the He I
excited-level channel (`:101-106`).  Summing the metal ground-state
recombination continua whose photons lie above 24.587 eV — C IV, N III,
N IV, O III, O IV, Ne III, Ne IV, S III, S IV — gives **0.29%** of the He I
ground-continuum channel, about 0.001 pc on the front
(`docs/HII40_CLOUDY_RESIDUAL_INVESTIGATION.md`, "Metal opacity").

So the omission is a justified approximation.  It is also an **undocumented**
one: there is no comment at the emission site recording that metal
recombination continua are omitted, why, or over what range the omission
holds.  The house rule requires one or the other — implement it, or record it
at the code site with its domain of validity.  This stage is that choice.

### 5.2 What the reference codes do

- **Cloudy transports both, and they genuinely re-ionize.**  Metal
  recombination continua are emitted as diffuse radiation
  (`rt_diffuse.cpp:196-270`) from the stored Milne-inverted metal cross
  sections; metal recombination Lyman-alpha/Balmer photons go into the
  on-the-spot field (`rt_ots.cpp:362-406`).  Both enter
  `rfield.SummedCon[i] = rfield.flux[0][i] + rfield.SummedDif[i]`
  (`:588-592`), and `SummedCon` is what `GammaK`/`GammaBn` integrate against
  the cross sections (`cont_gammas.cpp:125`, `:143`, `:409`) to produce
  `ionbal.PhotoRate_Shell` (`ion_photo.cpp:100-110`) — including feedback
  into hydrogen.
- **MOCASSIN transports the continua and discards the metal line photons.**
  `gammaHeavies` is built in `fb_ff` over `do elem = 3, 19`
  (`emission_mod.f90:474`, with `:472-473` noting `Z > 19` is neglected)
  through the Milne relation with the true subshell cross section `:519-521`,
  folded into the continuum `:126-128`, into `sumDiffuseHI` `:983` and the
  cumulative diffuse PDF `recPDF` `:1175-1196`, sampled by `newPhotonPacket`
  (`photon_mod.f90:958`) and propagated `:461-466`.  The metal lines are
  summed (`emission_mod.f90:1143-1164`, combined `:1169`, normalized `:1302`)
  and the line-versus-continuum draw uses that probability
  (`photon_mod.f90:925`), but in production the line identity is thrown away
  (`:926-941`, `nuP = 0` outside `lgDebug`) and the packet is dumped into
  `linePackets` and terminated without transport (`:468-475`).  `lgDebug`
  defaults `.false.` (`set_input_mod.f90:44`).
- **CMacIonize drops the metals entirely.**  `PhotonType.hpp:43-55` defines
  only `PRIMARY`, `DIFFUSE_HI`, `DIFFUSE_HeI`, `ABSORBED`;
  `IonizationVariables.hpp:46-61` has the H probability and four He channels.
  Metal line emission is only an energy sink in the cooling function
  (`TemperatureCalculator.cpp:463`).
- TORUS and Wood were **not read** on this axis.

### 5.3 Implementation shape, if implemented

Adding metal ground continua as a fifth channel follows the existing pattern
exactly, and one property makes it cheaper than it looks: the channel draw
already uses one Sobol dimension for an inverse CDF over the channels
(`u(7)`, `src/diffuse_mod.f90:217`) and the first energy variate is `u(8)`
(`:218-219`), so widening the channel list consumes **no new Sobol
dimension** — the same argument that let the He I 584 A split in with none
(`CLAUDE.md`, diffuse field).  The energies come from the same Milne
free-bound spectrum `E^2 sigma_pi(E) exp[-(E - Eth)/kT]` already used for
the three H/He ground continua, evaluated on the transport cross sections
`species_sigma` (`src/species_mod.f90:800-813`) — the same construction, on
the same data, that `par%diffuse_energy_model = 'milne'` uses.

### 5.4 Gate

**G-M3.**  Whichever branch is taken, the numbers to record are the change in
`r(He0 = 0.5)` at HII40 and the change in the total diffuse band luminosity,
against the measured 0.29% of channel 2.

Quantitative acceptance, derived from numbers already in the tree: the front
response to the diffuse field is calibrated at **0.233-0.585 pc per unit
fractional change** in the total diffuse band luminosity
(`docs/HII40_CLOUDY_RESIDUAL_INVESTIGATION.md`, "What is excluded").  Since
channel 2 is at most the whole diffuse band, 0.29% bounds the front motion at
**0.0007-0.0017 pc**.  If an implementation moves the front by more than
about 0.002 pc, the implementation is wrong, not the estimate — most likely
double counting the recombinations already carried in the metal chain's
`alpha_rec`, which is the specific trap here: the diffuse channel must emit
only the **ground-state** capture fraction, and the packet count must come
from the same coefficient the cascade recombines with, never from one
determination minus another (`CLAUDE.md`, physics lesson 2).

If the branch taken is documentation rather than implementation, the code-site
comment at the `dif_ch` declaration must carry the 0.29%, the nine emitting
stages, and the domain of validity: HII-region stages and abundances, where
the metal recombination continua above 24.587 eV are three parts in a
thousand of the helium channel.  It must also say what would break it —
higher metal stages with harder recombination continua, which is exactly what
`docs/HARD_SPECTRUM_PLAN.md` HS2 introduces.

### 5.5 Sub-item, scoped out with a measurement named

Metal **line** photons: MoCHII creates none, treating metal line cooling as a
pure energy sink.  Whether any of the registry's lines fall above the
13.598 eV band floor is **UNMEASURED**; the measurement is a sum over the
line list in `<base>_lines.txt` against `par%eion_min`, and it is cheap.  It
becomes a real question in FUV runs, where `par%add_fuv` puts the band floor
at 6 eV and the metal stages are the only gas opacity there
(`src/define.f90:565-566`).  Named here, not planned.

## 6. Stage M4 — the N x N stage matrix, deferred to the hard-spectrum work

### 6.1 Why it is forced rather than optional

The nearest-neighbor ratio chain is **exact**, not approximate, as long as
only adjacent stages connect: single-electron photoionization, radiative and
dielectronic recombination of one electron, and charge exchange of one
electron all move a stage by exactly one.  Auger multi-electron ejection
breaks that: an inner-shell ionization of stage `i` lands at stage
`i + n_e(shell)`, not `i + 1`.

That is precisely why Cloudy's matrix is deliberately not tridiagonal.  The
ionization block at `ion_solver.cpp:727-736` loops
`for(ion_to=ion_from+1; ion_to<=IonHigh; ...)`, filling arbitrary
super-diagonal entries, with the Auger reason stated at `:661-668`; one row
is replaced by abundance conservation `:883-889`.  And it is why Cloudy's own
legacy nearest-neighbor chain solver `solveions()` is dead code in c25
(`:1317-1364`, with the tridiagonal matrix drawn in the comment
`:1290-1315`).  Cloudy did not choose a matrix for elegance; it outgrew the
chain for the same reason MoCHII will.

**So this item is not optional once `docs/HARD_SPECTRUM_PLAN.md` HS6
introduces inner-shell absorption and Auger yields, and it should be done
there rather than as separate work.**  HS6 item 3 already says so ("Auger
couples non-adjacent stages, so the cascade becomes a (small)
transition-matrix solve for each element — reuse the n-level partial-pivot
solver on the stage populations").  This section records the solver decision
and the cost accounting so HS6 inherits them; it is not a parallel plan.

Until then, the chain is exact, and the correct action in this plan is the
documentation one: `species_fractions` (`src/species_mod.f90:815-849`) should
carry at the code site the statement of *why* it is exact and what breaks it
— adjacent stages only, valid while no inner-shell channel is active.  That
is a domain of validity for a construction that is currently presented as
simply "the product chain".

### 6.2 Solver and cost

**Reuse `nlevel_mod`'s Gauss elimination with partial pivoting**
(`src/nlevel_mod.f90:328-350`), which already solves a rate matrix plus a
conservation closure row for level populations
(`M(1,1:n) = 1`, `M(1,n+1) = 1` at `:326-327`) — the same algebraic object as
the stage matrix with abundance conservation, which is the row Cloudy
replaces at `ion_solver.cpp:883-889`.  It is inline inside
`nlevel_cooling_solve` today, so HS6 must factor it out.  Suggested name for
the extracted routine: `statistical_equilibrium_populations`, since that is
what it computes in both uses; and `ionization_equilibrium_stages` in
`species_mod` for the routine that fills the stage rate matrix, including the
Auger super-diagonal entries, and calls it.

Cost accounting, so HS6 does not have to redo it:

- Stage counts are small.  `MAX_ST = 6` today (`src/species_mod.f90:52`), and
  `docs/HARD_SPECTRUM_PLAN.md` HS2 targets 6 stages for O and Ne, so the
  matrix is `n x n` with `n` of order 6-8 and `MAX_ST` must be raised.
- Gaussian elimination is then of order `n^3/3`, i.e. 70-170 multiply-adds,
  against of order `2n` = 12-16 for the chain — one to two orders of
  magnitude more arithmetic in the linear algebra.
- **But both paths share the expensive part**, and it is already hoisted out.
  `species_ne_prepare` (`src/species_mod.f90:636-655`) evaluates
  `voronov_ci`, `alpha_rec` and `cx_rate` — the exponentials and powers —
  once for each `(leaf, T)`, and `species_ne_cached` reduces the inner
  iteration to multiply-adds.  The matrix path must reuse that same cache, in
  which case its added cost lands entirely on the cheap side of that split.
- The multiplier that makes this worth measuring rather than asserting: the
  thermal bisection calls the cell solve about 23 times for each leaf and gas
  iteration (`docs/COOLING_LOCAL_NE_PLAN.md`, section 2), so wall time for
  each gas iteration is the number to record, not flop counts.

**Gate G-M4 (inherited by HS6).**  With no inner-shell channel active, the
matrix path must reproduce the chain path to round-off on every element and
stage — the chain is the exact adjacent-stage limit, so this is a real
identity and not a tolerance.  Plus the recorded wall time of both paths on
HII40.

## 7. Withdrawn: a global Newton over all species

An earlier iteration of this analysis proposed a global Newton solve for each
leaf — a Jacobian of roughly 70 x 70 coupling every stage of every element
plus `n_e` plus `T_e`.  **That is withdrawn as a target, and the reason is
the survey: no code does it.**

MOCASSIN, CMacIonize, TORUS and Wood all use chains, one element at a time.
Cloudy builds a matrix but one for each element
(`conv_base.cpp:557-566`), with the coupling between elements, `n_e`, the
opacity and on-the-spot field, the chemistry network and `T_e` all carried by
outer iterations.  Its hierarchy, outermost first: iteration/zone loop
(`cloudy.cpp:197-243`) -> pressure loop
(`conv_pres_temp_eden_ioniz.cpp:120`) -> Brent root find on `T_e` with
residual `CoolHeatError()` (`conv_temp_eden_ioniz.cpp:91-197`) -> Brent root
find on `n_e` (`conv_eden_ioniz.cpp:44`) -> `ConvIoniz()` calling `ConvBase`
up to ten times (`conv_ioniz.cpp:17-25`) -> `ConvBase` do-while over
`loop_ion` (`conv_base.cpp:486-901`) -> inner `ion_loop < 3` (`:539`) -> the
matrix for one element.  Even `n_e` is not a matrix unknown there: it is a
summed diagnostic (`eden_sum.cpp:51-62`) whose used value changes only in
`EdenError()` (`conv_eden_ioniz.cpp:271-276`), the residual of a root find.

So the outer-iteration structure MoCHII already has — a scalar `n_e` fixed
point inside a `T_e` bisection inside the Monte Carlo iteration — is the
structure the field uses, at every level.  Recorded here so the proposal is
not raised again.  A global Newton would also lose the property the founding
decision was built on: with each element an independent chain, adding an ion
changes no existing results path.

## 8. What this plan does not expect to buy

**The unexplained 0.17750 pc HII40 helium-front residual against Cloudy c25
is unlikely to be closed by anything in this document.**  It was 0.18561 pc
until the `alpha_1` correction moved the front outward by 0.00811 pc.  The
full exclusion accounting is
`docs/HII40_CLOUDY_RESIDUAL_INVESTIGATION.md`; the parts that bear on metals:

- **The metal opacity lever is already measured and it saturates upward.**
  Removing metal absorption costs 0.615 pc; adding it gains only 0.103 pc per
  unit.  Closing 0.17750 pc on the upward branch needs **+180% metal
  opacity (x2.8)**, and the saturation means more than that.  An independent
  Opacity Project comparison through sirocco's TOPbase tables puts OP/VFKY96
  at C II 0.946, N II 0.978, O II 0.863 — 2-14% **lower**, the wrong
  direction, worth -0.005 to -0.03 pc.  The earlier claim that a ~30% metal
  opacity error would be worth ~0.18 pc is withdrawn there: 30% is worth
  0.031-0.072 pc.
- **Metal recombination continua as a diffuse source are already bounded at
  0.29% of channel 2, about 0.001 pc** — so M3 cannot close it, by
  construction, and Section 5.4 turns that bound into M3's acceptance
  criterion rather than its motivation.
- **M1's value is in the front and the PDR region, not in the ionized core.**
  This is now measured rather than argued (Section 3.5): at 3.8 pc, where the
  helium front sits, the frozen-state effect on `x_HI` is below 0.05%, against
  11.3% at the hydrogen front near 4.67 pc.  The
  residual, by contrast, is a *helium* signature that compounds with depth —
  the `x(He0)` ratio to Cloudy runs from 1.007 at 1.6 pc to 1.196 at 3.8 pc.
  The two are different objects.

The one item in this document that touches the residual's own territory is
the missing **helium**–hydrogen charge exchange of Section 3.7, which is
helium-specific and has no row in the exclusion table.  It is named as an
open item with an unmeasured magnitude and no rate data in the tree; it is
not a claim, and it is not proposed as an explanation.

The surviving hypothesis for the residual remains what the investigation left
it as: several few-percent differences compounding together, which no
single-parameter run can produce or exclude.  Nothing in M1-M4 changes that.

## 9. Where MoCHII stands, axis by axis

"not read" means the survey did not cover that code on that axis; it is not a
statement that the feature is absent.

| axis | MOCASSIN | CMacIonize | TORUS | Wood | Cloudy c25 | MoCHII now | after this plan |
|---|---|---|---|---|---|---|---|
| stage populations | chain | chain | chain | chain, one pass, no iteration | N x N matrix + LU, one per element | chain | chain; matrix when an Auger channel is active (M4, in HS6) |
| elements coupled to each other | no | no | no | no | no | no | no (Section 7) |
| charge exchange two-way in the bookkeeping | no | no | no | no | **yes** | **no** | **yes (M1)** |
| metal electrons in `n_e` | yes | no | yes | no | yes (summed, root find) | switch, default off, **on** in the published runs | measured, default decided (M2) |
| metal opacity attenuates the ionizing field | yes | no | no | no | yes | yes, default on | unchanged |
| metal photoheating | yes | no | not read | no | yes | switch, default off, **on** in the published runs | measured, default decided (M2) |
| metal recombination continua transported | yes | no | not read | not read | yes | no (0.29% of channel 2) | implemented, or documented with its domain of validity (M3) |
| metal line photons transported | no, discarded in production | no, energy sink only | not read | not read | yes, on the spot | no | measurement named, not planned (M3.5) |
| `T_e` | outer root find | outer root find | outer root find | outer root find | outer Brent | outer bisection | unchanged |
| charge-exchange heating/cooling | not read | not read | not read | not read | not read | no | flagged (M1.7) |

## 10. Traps in the reference codes

Any validation against these codes must correct for the following.  Each is a
property of the distributed source, cited.

**MOCASSIN.**

- Missing dielectronic rates are substituted with **oxygen's**:
  `update_mod.f90:1648`,
  `if (diRec(elem,ion) == 0.) diRec(elem,ion) = diRec(8,ion)`.
- Convergence is tested on **H0, He0, He+ only** (`:1598-1604`,
  `limit = 0.01` at `:1337`) — no metal ion fraction is ever tested, so a
  "converged" MOCASSIN model carries no statement about its metal
  populations.
- `nstages = 7` (`set_input_mod.f90:82`) with `nElements = 30`
  (`constants_mod.f90:54`), so Fe truncates at Fe VII and no stage above 7
  can contribute anywhere, including to `n_e`.
- Metal **line photons are discarded** in production
  (`photon_mod.f90:926-941`, terminated `:468-475`); the line PDF and
  `grid%Jdif` exist only under `lgDebug` (`emission_mod.f90:1209`, `:1251-1289`;
  `photon_mod.f90:1491-1493`, `:1730-1732`), and `nPhotoDif` enters the ion
  balance only in the `lgDebug` branch (`update_mod.f90:1489` versus `:1495`).
  `lgDebug` defaults `.false.` (`set_input_mod.f90:44`).
- Ions with fraction below `1.e-10` are silently dropped from `n_e`
  (`update_mod.f90:1564-1577`).
- Charge transfer is **hard-zeroed outside 6000-50000 K** (`:1466`), and the
  reverse rate is nonzero only for N I and O I (`:1448-1449`, gated `:1475`).
- No excited-level H/He opacity (comment `ionization_mod.f90:398-401`), and
  stimulated emission is dead code because every `inOpacity` call passes
  `b = 0.` (`:467-472`).
- **`alpha_A(He I)` is Verner & Ferland (1996) with no special case for H/He,
  4.8394e-13 at 1e4 K — +11% against Cloudy's level-by-level Milne value**
  4.3583e-13, which MoCHII's Badnell 4.3719e-13 matches to 0.3%.  Subtract
  that 11% before attributing any helium front or line-intensity
  disagreement to transport or geometry (`CLAUDE.md`;
  `docs/HII40_CLOUDY_RESIDUAL_INVESTIGATION.md`, "The He ionization
  network").

**TORUS.**

- All four "Ar I-IV" ions are created with `n = 18`, i.e. **all neutral**
  (`ion_mod.F90:331-340`).  Argon results are not usable.
- The opacity omits **even He II**: `amr_mod.F90:16301-16309` calls
  `phfit2(1,1,1,...)` and `phfit2(2,2,1,...)` only, and
  `kappaAbsGas = kappaH + kappaHe` `:16311` sets the photon's optical depth
  (`photoionAMR_mod.F90:4762`, `:4741-4743`).  Metal cross sections exist
  (`ion_mod.F90:172`) but only build `photoIonCoeff` (`:6688-6692`), so metal
  photoionizations remove no photons from the beam.
- 22 ions hard-coded (`ion_mod.F90:201-294`), with higher stages commented
  out for lack of level data (`:245-247`, `:233`, `:262`, `:277`);
  `xraymetals` (Mg/Si/Ar/Fe) defaults `.false.`
  (`inputs_modV2.F90:3967`).
- Under-relaxation 0.6 on both the ionization and the thermal update
  (`photoionAMR_mod.F90:6853`, `:6914`, `:6246`), so a comparison must
  confirm the fixed point was reached and not merely damped.

**CMacIonize.**

- **No metal electrons** (`IonizationStateCalculator.cpp:134`, comment
  `:273-274`, also documented `:576-586`), **no metal opacity**
  (`DensityGrid.hpp:117-140`, task path identical at
  `DensitySubGrid.hpp:557-583`), **no metal photoheating**
  (`TemperatureCalculator.cpp:297`).  All three feedbacks are absent, so
  CMacIonize is the *least* coupled of the five on the metal axis, not the
  most.
- Mean-intensity integrals are nevertheless accumulated for every ion
  including the metals (`DensityGrid.hpp:162-181`), so metal
  photoionizations consume no ionizing budget — diagnostic-only sinks.
- Argon is enum-only dead weight (`ElementNames.hpp:142-151`, `HAS_ARGON`
  only under `ADDITIONAL_COOLANTS`, `:46-49`).
- The top stage of each element is implicit (`ElementNames.hpp:96-99`,
  recovered `TemperatureCalculator.cpp:363-367`) — stage counts are C 3,
  N 4, O 3, Ne 3, S 4.

**Wood `ionize`.**

- The thermal solve runs only `if(iter.ge.5)` (`itbal.f:27-41`), so the code
  is **effectively isothermal at `Tconst = 8000` K** (`input.params:6`,
  passed at `mcionize.f:115`) for its first four Monte Carlo iterations, and
  `ionize2.f:30` takes `T = temp(i,j,k)` as given.  A comparison must confirm
  the iteration count exceeded that.
- Metal fractions are one pass per cell per iteration with **no iteration at
  all**, each metal block reading only `ne`, `h0`, `he0` and never another
  metal's state.
- **S I = C I = 0 assumed** (`grid.txt:21-22`); no Fe, Si, Ar, Mg, Cl.
- No metal electrons (`ioneng.f:45`, `ionize2.f:54`; deliberate,
  `tbal1.f:39-41`), and the He term is single-ionization only, so He++
  contributes nothing.
- No metal opacity (`tauint.f:159-161`, sub-cell step-back `:169-171`) while
  the metal cross sections are computed every frequency (`opacset.f:20-55`)
  and accumulated into the metal rates in the same loop (`tauint.f:194-216`).
  No dust anywhere in the code.
- Heating is **H/He photoionization only** (`ioneng.f:53-54`, comment
  `tauint.f:180-182`); PAH heating off (`pahfac = 0`, `input.params:17`) and
  the Reynolds extra heating disabled (`G1 = 0.0e-25*ne*1.e20`,
  `ioneng.f:64`).
- Hard temperature floor: `if(t0.lt.4000.) t0=500.` with forced neutrality
  (`tbal1.f:109-115`).

**Cloudy c25.**

- Within a single LU the charge-transfer partner's density is a **frozen
  coefficient**, so exact charge conservation is reached by the outer
  iteration plus the `"CX inconsistency"` check
  (`conv_base.cpp:578-583`, `:611-620`), not inside one matrix solve.  Cite
  Cloudy as two-way in the bookkeeping, not as exactly conserving inside a
  solve.
- The legacy chain solver `solveions()` is dead code (`ion_solver.cpp:1317-1364`);
  reading it as Cloudy's method is an error.

## 11. Open questions and unverified points

Collected so a later session does not have to re-derive which claims rest on
measurement and which do not.

1. *(measured while drafting)* M1's magnitude was estimated on the frozen
   converged HII40 state (Section 3.5): 0.005-0.05% on `x_HI` through the
   ionized zone, 11.3% at the front, and over 99.99% of it O0/O+ against
   H+/H0.  What remains unmeasured is the **self-consistent** answer, since
   the estimate does not re-converge `x_HI -> tau -> Gamma` or shift oxygen's
   own fractions; that is what G-M1c measures.
2. **Helium–hydrogen charge exchange is entirely absent from MoCHII**, and no
   rate exists in the tree to evaluate it (Section 3.7).  Helium-specific,
   with no row in the residual exclusion table.  Magnitude unmeasured.
3. **Charge-exchange heating and cooling**: not read in any of the five
   codes; absent from MoCHII (Section 3.7).
4. **TORUS on metal photoheating, metal recombination continua and metal
   line photons**: not read.  **Wood on metal recombination continua and
   metal line photons**: not read (Section 9).
5. **Whether any registry line lies above the 13.598 eV band floor**:
   unmeasured; the measurement is named in Section 5.5.
6. **The switch-matrix numbers of G-M2 do not exist anywhere in the tree**
   (Section 4).  Until they do, no statement about what the trace separation
   costs is supported.
7. **Three of the 22 charge-exchange recombination entries are placeholders**
   (`1.0e-14 cm^3 s^-1`, no temperature dependence: Ne, S, Ar), which the
   provenance headers should label as such rather than as determinations
   (Section 3.5).
8. *(closed while drafting)* The `src/species_mod.f90` module header stated
   the opposite of the shipped configuration; it was rewritten and now records
   that charge exchange is the one-way term (Section 2).

## 12. Stage M5 — documentation

Not a physics stage; listed so the ripple is done once.

- ~~`src/species_mod.f90` header: replace the "metals never feed back"
  paragraph~~ — done while drafting (`:11-20`).
- `docs/PLAN.md:236-241` and `:315`: record that the trace *feedback*
  non-goal is retired by measurement while the trace *solver structure* is
  kept and confirmed by the survey (Section 2).  Do not delete the founding
  text; annotate it, the way `docs/COOLING_LOCAL_NE_PLAN.md` section 12
  annotates its own predictions.
- `src/species_mod.f90:815-849` (`species_fractions`): the domain of validity
  of the chain (Section 6.1).
- `src/diffuse_mod.f90:40` (`dif_ch`): the metal continuum omission, if M3
  takes the documentation branch (Section 5.4).
- `src/define.f90:568-583`: the outcome of the M2 default decision
  (Section 4.3).
- `data/atomic/README.md`: the placeholder charge-exchange entries.
- `docs/MoCHII_physics.tex`: the M1 formula, the G-M1a gate, the G-M2 switch
  matrix, and the `eps` diagnostic of 3.4.  `docs/MoCHII_UserGuide.tex`: the
  `par%charge_exchange` row.  `paper/mochii.tex`: the paragraph that already
  points here (Section 1) gains the M1 result and the M2 numbers.
- `CLAUDE.md`: status and next steps.

## 13. Order, effort, delegation

| order | stage | risk | blocking? |
|---|---|---|---|
| 1 | M1 two-way charge exchange | low (two added terms, one new routine, no new module) | nothing; it is a correctness fix |
| 2 | M2 switch matrix | none (runs only) | the manuscript's numbers; the default decision |
| 3 | M3 metal recombination continua | low (existing channel pattern, no new Sobol dimension) | closes a documentation obligation |
| 4 | M4 stage matrix | **high** | belongs to `docs/HARD_SPECTRUM_PLAN.md` HS6; do it there |
| 5 | M5 documentation | none | after M1-M3 land |

M2 can run in parallel with M1 provided the M2 runs are made on the
pre-M1 binary and labeled with it — otherwise the two effects mix and
neither is measurable.

**Advisor keeps:** the algebra of 3.2 (the denominator is the part that gets
this wrong), the M2 default decision of 4.3, the M3 branch choice, reading
the diffs, and running G-M1a, G-M1c and G-M2 personally.

**Workers (Opus).**  Each brief must open with "read `~/.claude/CLAUDE.md`
and the project `CLAUDE.md` and follow every rule", and must carry the
forbidden-word list and the physics-based naming rule verbatim.

- **W1 = M1**: the two routines of 3.3, the `solve_ion_cell` change, the
  `par%charge_exchange` switch and its `setup.f90` validation, the
  `tests/charge_exchange/` unit test, and gates G-M1a, G-M1b, G-M1e.  The
  brief must include the closed form of 3.2 **with both terms in the
  denominator**, the 29-reaction census, and the requirement that
  `'metal_only'` reproduce the pre-change binary bit for bit.
- **W2 = M2**: the eight runs and the table.  No source changes; if a source
  change looks necessary, stop and report.
- **W3 = M3**: depends on the branch decision, not on W1 or W2.
