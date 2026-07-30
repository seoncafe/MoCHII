# Cell-solve convergence with two-way charge exchange

- Date: 2026-07-31
- Owner item: `docs/METAL_COUPLING_PLAN.md` Stage M1, gate G-M1f
- Status: analysis and options; **nothing implemented**

Stage M1 put the metal charge-exchange terms into hydrogen's balance. That is
correct physics and it is gated (G-M1a, all 29 reactions to 6.5e-16), but it
made `solve_ion_cell` converge slowly in partially neutral gas. This document
records why, what the reference codes do, and what the options are, so the fix
is chosen on evidence rather than on the first idea.

## 1. The problem, stated exactly

`solve_ion_cell` (`src/ion_balance_mod.f90`) iterates hydrogen by substitution,

```
x <- F(x) = (A + R_in(x)) / (B + R_in(x) + R_out(x))
A = alpha_H n_e ,  B = Gamma_H + (C_H + alpha_H) n_e
```

with `R_in` and `R_out` evaluated from the previous iterate because the metal
stage fractions depend on `n_HI` and `n_HII` themselves.

**Measured.** At T = 1e4 K, n_H = 100, HII20 abundances, the iteration meets its
exit test after 124 passes at `x_HII = 0.18`, 239 at 0.10 and **1285 at 0.021**,
against the hard 200-pass cap. Capped, a cell leaves with `x_HII` high by about
8% (2.252e-2 against 2.078e-2) and a residual near 1e-3. `metal_only` shows
none of it (2e-10 everywhere).

**Why.** Where the near-resonant O0 + H+ pair locks oxygen to hydrogen -- that
is, where `Gamma_O` and `n_e alpha_O` are negligible beside the charge-exchange
terms -- detailed balance of the pair gives

```
f(O+)/f(O0) = (n_HII k_cxi)/(n_HI k_cxr)   =>   R_out/R_in = n_HII/n_HI = (1-x)/x
```

Substituting that into `F` and solving returns, **exactly**,

```
x = alpha_H n_e / (Gamma_H + (C_H + alpha_H) n_e)
```

the no-charge-exchange solution: the two terms cancel identically (verified
symbolically). So M1's net effect is only the *departure* from that lock, and
the slope of the substitution map at its fixed point is, in the same limit,

```
dF/dx = R_in / (R_in + alpha_H n_e)
```

which tends to 1 wherever `R_in >> alpha_H n_e`, i.e. exactly the neutral gas.
The iteration is **marginally stable by construction**, and the pass count
follows `ln(tol)/ln(slope)`.

**What it affects.** Cells with `x_HII >~ 0.1` converge inside the cap, so the
front and everything inside it are unaffected, and the Lexington diagnostics are
weighted by `n_p n_e ~ x_HII^2`, so `<T[NpNe]>` and `L(Hbeta)` are safe (G-M1c
measured -1.41 K and -0.13 K). **At risk is the shadow beyond the front and any
PDR application** -- which is precisely where M1 was expected to matter.

## 2. What the reference codes do

Read in the sources this session. The short answer: **not one of the five uses
Newton or secant on the ionization fixed point.** All five iterate the
ionization state by substitution, and reserve root finds for the scalar outer
variables.

| code | ionization state | `n_e` | `T_e` |
|---|---|---|---|
| Wood `ionize` | substitution, H and He each solved in closed form given the other; damping after 10 passes, hard stop at 20 (`h0find.f:38-80`) | H/He only, closed form | log secant (`tbal1.f:78-117`) |
| CMacIonize | substitution, inner fixed point (`IonizationStateCalculator.cpp:699-752`); metals one pass afterwards | H/He only, closed form | power-law secant on three trial temperatures (`TemperatureCalculator.cpp:730-753`, `:795`) |
| MOCASSIN | substitution by recursion, cap `maxIterateX = 25` (`update_mod.f90:1321`, `:84-85`, `:1608`), convergence tested on H0/He0/He+ only (`:1598-1604`) | summed, inside that recursion (`:1564-1577`) | bracketing secant / regula falsi, cap 25 (`:340`, `:399-455`) |
| TORUS | substitution with fixed under-relaxation 0.6 (`photoionAMR_mod.F90:6859-6904`, `:6853`, `:6914`) | summed (`ion_mod.F90:1418-1421`) | ZBRENT (`photoionAMR_mod.F90:6185-6241`), then under-relaxed 0.6 |
| Cloudy c25 | **no iteration**: an N x N matrix per element solved by LU with iterative refinement (`ion_solver.cpp:233`, `newton_step.cpp:346-400`) | **van Wijngaarden-Dekker-Brent root find** (`conv_eden_ioniz.cpp:44-150`, residual set only in `EdenError()`, `:271-276`) | Brent root find (`conv_temp_eden_ioniz.cpp:91-197`) |
| MoCHII today | substitution, damped 0.5 on `n_e`, cap 200 | inside the same substitution | bisection |

Two precedents matter for us.

**Cloudy solves the scalar MoCHII substitutes.** `n_e` is not a fixed-point
iterate there; it is the root of `EdenError()` found by a safeguarded bracketing
method. That is exactly the shape of option B below, and it is not exotic --
it is what the reference code does for the same variable.

**Wood already guards a near-cancellation by reformulating, not by iterating
harder.** `h0find.f:66-72` solves the quadratic for `h0` in closed form, but
when the discriminant term `t1 = 4 ch^2 (1 + AHe - AHe he0)/b^2` falls below
1e-3 it does **not** evaluate `(b - sqrt(b^2 - ...))/(2 ch)`; it switches to the
first-order expansion `h0 = ch (1 + AHe - AHe he0)/b`, because the difference of
two nearly equal numbers loses precision. That is the same disease as ours, met
with the same medicine as option C.

### 2.1 How charge exchange itself enters each iteration

The table above is about the outer solve.  The question that matters for M1 is
narrower: **with which iterate of the hydrogen state is the charge-exchange term
evaluated, and does the pair participate in the convergence test?**  Every one of
the five lags it -- none evaluates it implicitly -- but they lag it by very
different amounts, and only one of them checks the result.

**Wood `ionize` -- lagged by a whole Monte Carlo iteration.**  `h0find` solves
H and He with no charge-exchange term at all, and the metal blocks in
`ionize2.f` then run once, reading the converged `h0`, `he0` and `ne` of that
call (`:136-138` for O, `:83-85` for N).  There is no iteration on the metal
fractions and no feedback within the call, so the coupling closes only through
the next radiation pass.

**CMacIonize -- the same, and structurally so.**  Metals are computed after the
H/He state is frozen (`IonizationStateCalculator.cpp:116-185`, and the same
ordering in `TemperatureCalculator.cpp:270-345`), one pass, using that state.
Since the H/He solve carries no charge-transfer term either (no hit above
`:501`), the pair never appears on both sides in the same solve.

**TORUS -- lagged by one pass of a per-octal alternation.**  The charge-exchange
rates read the current `nHI`/`nHII` (`photoionAMR_mod.F90:6873`, `:6877`) inside
the loop that alternates the ionization and thermal solves (`:3909-3915`), and
the whole update is under-relaxed by 0.6 (`:6853`, `:6914`).  The under-relaxation
is applied to the ion fractions, not to the charge-exchange term specifically.

**MOCASSIN -- lagged by one pass of its own recursion, and switched off in cold
gas.**  The metal denominator reads the stored H0 fraction
`grid%ionDen(cellP, elementXref(1), 1)` (`update_mod.f90:1498-1499`), which the
previous recursion pass wrote, so the lag is one pass of `ionBalance` -- the
same structure as MoCHII's.  Two differences are worth carrying into any
comparison.  The `chex` table is repopulated inside the recursive subroutine on
every pass (`:1392-1451`), and, more consequentially, **the rate is hard-zeroed
outside 6000-50000 K**: `if (TeUsed < 6000. .or. TeUsed > 5.e4) chex(elem,ion,1) = 0.`
(`:1466`).  Below 6000 K MOCASSIN therefore has no charge exchange at all, where
MoCHII clamps the Kingdon & Ferland forms to the edge of their fit range and
keeps evaluating (`cx_rate`, forms 6 and 7).  A cold-gas or PDR comparison
against MOCASSIN is comparing two different physics inputs, not two solvers.
And no metal ion enters MOCASSIN's convergence test, which is on H0, He0 and
He+ only (`:1598-1604`), so a badly converged charge-exchange partner cannot
hold the cell open.

**Cloudy -- lagged inside one matrix, closed by the outer loop, and checked.**
Within a single LU the partner's density is a frozen coefficient, exactly as in
a lagged substitution.  What differs is everything around it.
`ChargTranEval()` is called once per pass of `ConvBase`
(`conv_base.cpp:502`) with the comment stating the reason -- "charge transfer
evaluation needs to be here so that same rate coefficient used for H ion and
other recombination" -- so both sides of every reaction are guaranteed to read
one number, which is the property MoCHII gets instead from the shared cache.
`iso_charge_transfer_update()` then re-accumulates hydrogen's totals over every
heavier element at the top of each ion solve (`ion_solver.cpp:68`).

Then Cloudy does the thing none of the others do: it forms the **net** rate of
each pair and refuses to call the zone converged if the two directions
disagree, `conv_base.cpp:578-583` and `:605-620`.  Two details of that test are
directly relevant to us.  The tolerance is adaptive and is written against the
cancellation --

```
ion_cmp = MAX2( 0.01*MIN2(RateIonizTot(nlo)*dl1, RateIonizTot(nhi)*dl2),
                1e-4*CharExcRecTo[nlo][nhi][0]*dl1*ul2 )
```

that is, the net must be small against the *smaller* of the two ionization
rates, not against either term of the difference.  And the check is applied to
**one pair only**: `lgCheckAll` is false, so the condition reduces to
`nlo == ipHYDROGEN && nhi == ipOXYGEN`.  Cloudy hard-codes the
hydrogen-oxygen pair as the one worth watching.

That is an independent confirmation of the analysis in section 1.  MoCHII's own
measurement found the same pair carrying over 99.99% of both sums at the
ionization front (2.835e-11 and 2.838e-11 s^-1, agreeing to four digits), and
identified that near-cancellation as the origin of the marginally stable mode.
The reference code arrived at the same place from the other direction: it is the
only pair whose bookkeeping it bothers to police.

**What MoCHII should take from this.**  Not the lag -- MoCHII's is already the
tightest of the five, since `R_in` and `R_out` are formed from the same cached
coefficients and the same chain walk as the electron sum, inside the `n_e` fixed
point rather than outside it.  What is missing is the second half: Cloudy pairs
that lag with a **net-flux test that knows about the cancellation**.  MoCHII's
`hydrogen_balance_residual` currently reports the residual of the whole balance;
adding the O-H pair net flux, compared against the smaller of the two rates in
Cloudy's manner, would make the diagnostic sensitive to precisely the failure
mode section 1 predicts, and is worth doing whichever of A, B or C is chosen.

## 3. Options

### Option A -- Aitken/Steffensen acceleration on the existing sequence

Take three consecutive substitution iterates and extrapolate:

```
x_acc = x0 - (x1 - x0)^2 / (x2 - 2 x1 + x0)      ! guard the denominator
```

applied every third pass, falling back to plain substitution when the
denominator is small or `x_acc` leaves `[0, 1]`.

- **For**: three lines; no derivatives; no new physics evaluation; turns a
  linearly convergent sequence with ratio -> 1 into a quadratically convergent
  one; degrades to the present behavior when it does not help, so the failure
  mode is "no worse".
- **Against**: relies on the error being dominated by a single mode. Here the
  slope analysis says it is (one eigenvalue near 1, set by the oxygen pair), but
  that is an argument, not a measurement, and it will be weaker where two
  elements contribute comparably. Does not remove the cancellation, only
  extrapolates through it, so the precision loss in `R_in - R_out` remains.

### Option B -- safeguarded secant on the scalar residual

Solve `R(x) = 0` instead of iterating `x <- F(x)`:

```
R(x) = x (B + R_in(x) + R_out(x)) - (A + R_in(x))
```

with the secant update, bracketed on `[0, 1]`, bisecting whenever the secant
step leaves the bracket or fails to reduce `|R|` -- i.e. a Brent-shaped method,
which is what Cloudy uses for `n_e`.

- **For**: convergence independent of the slope, so the slow mode disappears
  outright; superlinear (order ~1.618); **no derivatives**, and it reuses the
  `R_in`/`R_out` evaluation already in the loop, so no new physics code. It is
  the reference-code precedent.
- **Against**: substitution here is a self-map of `[0, 1]` and therefore
  **cannot** return an unphysical state; secant can, and this equation is known
  to admit unphysical roots -- G-M1a produced `x_HI = 6.47` and `578` when the
  denominator was wrong. Bracketing is therefore mandatory, not optional. Needs
  two starting iterates, so the first pass stays a substitution step. Kinks from
  `max(rrec, tinest)` and the `if (cch_cxi > 0)` branches are tolerable for
  secant but not free.

### Option C -- remove the cancellation analytically

Since the locked limit is known in closed form, solve for the **departure** from
it rather than for `x` itself: write `R_out = R_in (1-x)/x + D(x)`, where `D`
collects the terms that break the lock (`Gamma_O`, `n_e alpha_O`, and the other
elements), substitute, and let the large parts cancel on paper. What is left is
an equation for `x` whose stiffness is set by `D`, not by `R_in`.

- **For**: treats the cause. Both the slow convergence and the precision loss in
  the difference of two large terms disappear together, which no accelerator
  fixes. Wood's `sqrt` expansion is the same move.
- **Against**: needs a careful re-derivation, and the cancellation is only exact
  for one reaction pair at a time -- the algebra has to hold with 29 reactions
  and 11 elements, not just oxygen. Highest risk of a subtle error, and it
  changes the expression that G-M1a checks, so that gate must be re-derived
  alongside.

### Option D -- raise the cap, or add `|dx_HI|` to the exit test

- **Against, decisively**: 1285 passes for each of about 23 temperature trials
  for each leaf and each gas iteration is not affordable, and the exit test only
  *detects* the failure. These treat the symptom. Recorded so they are not
  proposed again as fixes; adding `|dx_HI|` to the exit test is still worth
  doing as a **diagnostic**, so a capped cell is visible rather than silent.

### Not an option -- special-casing oxygen

Solving the H-O pair as a 2x2 would work at HII40, where oxygen is over 99.99%
of both sums, and would break the generic registry the founding design is built
on (`docs/PLAN.md`). Rejected.

## 4. Recommendation

**A first, then B if A is not enough, and C only with a re-derived G-M1a.**

The order is by risk, not by elegance. A is nearly free and reversible; B is the
reference-code precedent and removes the slow mode but needs safeguarding
against a failure mode this equation demonstrably has; C is the right answer in
principle and the easiest to get subtly wrong.

**Measure before choosing.** The slope derived above is for the `x` map at fixed
`n_e`, but `x` sits inside the `n_e` fixed point with helium solved alongside;
what governs the real pass count is the spectral radius of the coupled map. If
the slow mode is partly in `n_e`, accelerating `x` alone moves the bottleneck
instead of removing it. **Instrument first: report the pass count and the
per-variable contraction ratio for a cell at `x_HII` near 0.02.**

## 5. Gates

- **G-M1f** (from `METAL_COUPLING_PLAN.md`). A cell at `x_HII` between 0.01 and
  0.05 must converge to the same state as repeated whole-cell solves, to
  <= 1e-10 relative, **inside** the pass cap.
- **G-M1g** (new, for any of A/B/C). `par%charge_exchange = 'metal_only'` must
  remain bit-identical to the pre-M1 binary, and the `two_way` HII40 and HII20
  results must not move by more than the run's own convergence tolerance
  (82 K and 69 K on `<T[NpNe]>`) -- the accelerator must change the path to the
  root, not the root.
- **G-M1h** (option B only). No leaf may leave the solve with `x_HI` outside
  `[0, 1]`, asserted in the solver rather than checked afterwards.

## 6. What this does not change

The physics of M1 is settled and gated; this is solely about how the fixed point
is walked. `'metal_only'` remains available and bit-identical, so a recorded
gate can always be reproduced while this is open.
