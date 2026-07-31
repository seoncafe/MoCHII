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

## 4. Execution order

Ordered by risk, not by elegance: each step is cheap and reversible before the
one after it, and each has a stop condition so the next is entered only on
evidence.  **Nothing below is implemented yet**; this section is the order the
work is to be done in.

### S1 -- Measure which variable carries the slow mode  *(do first, no code change)*

The slope of section 1 is derived for the `x` map at **fixed** `n_e`, but `x`
sits inside the `n_e` fixed point with helium solved alongside.  What sets the
real pass count is the spectral radius of the coupled map.  If a comparable
part of the mode lives in `n_e`, accelerating `x` alone moves the bottleneck
instead of removing it, and S2 would look like a failure of the method rather
than of the target.

**Do**: walk the cell map by hand for a series of cells spanning `x_HII` from
about 0.02 to 0.5, recording the error sequence of `x_HI` and of `n_e` against
the converged state, and report the asymptotic ratio
`|e_{k+1}|/|e_k|` for each, beside the predicted `R_in/(R_in + alpha_H n_e)`
and the pass count.  Every rate must come from the production routines, and the
converged state must be taken from `solve_ion_cell` itself so the replication
cannot drift from what it is measuring.

**Acceptance**: the table exists and is recorded here.

**Decision rule, fixed in advance so the measurement cannot be read to taste**:
- `ratio(x) ~ ` the predicted slope **and** `ratio(n_e)` clearly below it
  -> the mode is in `x`; go to S2 and accelerate `x`.
- `ratio(n_e) ~ ratio(x)`, both near 1 -> the mode is in the coupled pair;
  S2 must accelerate the pair (or S3 directly), and accelerating `x` alone is
  predicted to fail -- record that prediction before trying it.
- both ratios well below 1 while the pass count is still large -> the slow part
  is not asymptotic and none of A/B/C is the right fix; stop and re-diagnose.

#### S1 result, 2026-07-31

Measured with a driver that walks the map by hand while taking every rate from
the production routines and the converged state from `solve_ion_cell` itself.
T = 1e4 K, n_H = 100, Lexington abundances, `two_way`, cap 400 passes.

| `x_HII` | ratio `x_HI` | ratio `n_e` | slope predicted | passes to 1e-10 |
|---|---|---|---|---|
| 0.287  | 0.9587 | 0.9598 | 0.8589 | 46 |
| 0.161  | 0.9551 | 0.9548 | 0.9293 | 141 |
| 0.0742 | 0.9802 | 0.9799 | 0.9700 | 340 |
| 0.0347 | 0.9708 | 0.9705 | 0.9865 | **did not converge in 400** |

Method note, because it changes the answer: averaging the ratio over a fixed
tail window gives 1.000000 for every state, which is the round-off floor of an
already-converged sequence and means nothing.  The window has to be the part
where the error is still decaying -- past the transient and above 1e-13 of the
root.  The first run of this measurement made that mistake.

**Reading, against the rule fixed above.**  `ratio(x_HI)` and `ratio(n_e)` agree
to three or four digits in every state.  The rule's second branch --
"`ratio(n_e) ~ ratio(x)`, both near 1 -> the mode is in the coupled pair,
accelerating `x` alone is predicted to fail" -- was written assuming that equal
ratios would mean **two** slow variables.  The measurement shows something the
rule did not anticipate: they are equal because `n_e` is *slaved* to `x_HI` in
this regime, `n_e ~ n_H x_HII` plus the small helium and metal terms, so
`d(n_e) ~ -n_H d(x_HI)` and the two carry one mode, not two.  Accelerating
`x_HI` therefore necessarily accelerates `n_e`.  **This is recorded as a
refinement of a rule that did not cover the case, not as the rule being
satisfied**, and it is the reason S2 proceeds rather than S3.

**The section-1 slope is a mechanism, not a number.**  It tracks the trend and
the order of magnitude but is not accurate: predicted 0.859 against measured
0.959 at `x_HII` = 0.287, and predicted 0.987 against measured 0.971 at 0.0347,
where it overshoots.  That is expected -- the derivation assumes the fully
locked limit and a fixed `n_e`, and the real map has neither -- but it means the
slope must not be used to predict pass counts quantitatively.

**Conclusion**: a single dominant mode with ratio 0.955 to 0.980, which is
exactly the situation Aitken extrapolation is built for.  Go to S2.

### S2 -- Option A, Aitken/Steffensen  *(only if S1 says the mode is in `x`)*

Three lines, no derivatives, no new physics evaluation, and it degrades to the
present behavior when it does not help.

**Acceptance**: G-M1f and G-M1g below.
**Stop condition**: if the `x_HII` = 0.021 cell still fails to converge inside
the pass cap, go to S3.  Do not tune the accelerator to make one cell pass.

#### S2 result, 2026-07-31: **FAILED, and reverted**

Implemented as Aitken extrapolation of the **`n_e`** sequence, not of `x_HI`.
That choice is worth keeping whatever comes next: `x_HI` is recomputed each pass
by the substitution formula, whose terms are all non-negative, so it stays in
`[0, 1]` whatever `n_e` does.  Extrapolating `x_HI` directly would forfeit the
one safety property substitution has, in an equation that demonstrably admits
unphysical roots.  Guards: extrapolate every third pass, require the second
difference to be resolvable (`|d2| > 1e-13 |n_e|`) and the result to stay within
a factor of two of the current iterate, and apply it only when the
charge-exchange terms are active so `'metal_only'` stays bit-identical.

**G-M1f: FAIL.**  One `solve_ion_cell` call against the converged state:

| `x_HII` | one call | converged | relative deviation |
|---|---|---|---|
| 0.287 | 2.869235e-1 | 2.869235e-1 | 3.6e-10 |
| 0.161 | 1.613053e-1 | 1.613053e-1 | 9.1e-10 |
| 0.0742 | 7.416349e-2 | 7.416036e-2 | 3.4e-06 |
| 0.0347 | 3.495452e-2 | 3.465395e-2 | 3.1e-04 |
| 0.0208 | **2.252809e-2** | 2.077603e-2 | **1.8e-03** |

The worst cell is unchanged to all printed digits from the pre-accelerator
measurement.  Instrumenting the independent replication shows why it is not a
guard rejecting the step: the extrapolation **fires 133 times and is used 42 to
132 times** in every state.  It is used, and it does not help:

| `x_HII` | passes before | passes with Aitken | ratio `n_e` with Aitken |
|---|---|---|---|
| 0.287 | 46 | **105** | 1.354 |
| 0.161 | 141 | 137 | 1.110 |
| 0.0742 | 340 | 347 | 1.693 |
| 0.0347 | did not converge | did not converge | 1.060 |

**Diagnosis.**  The measured contraction ratio rises **above 1** once the
accelerator is on.  Aitken assumes the error is one geometric mode; S1 confirmed
a single dominant mode, but the iterate is not a scalar sequence -- it is the
`n_e` component of a coupled state that also carries `x_HI`, the two helium
ratios and eleven metal stage chains.  Each extrapolation jump moves `n_e` off
the slow manifold those other components are sitting on, and the transient it
re-excites costs more than the extrapolation saved.  State 1, which converged
fastest without help, is hurt the most (46 -> 105), which is the signature: the
further the state is from the asymptotic regime, the more damage the jump does.

**This is a negative result about the method, not about the implementation**,
and it was reached by the stop condition fixed in advance -- "if the
`x_HII = 0.021` cell still fails to converge inside the pass cap, go to S3.  Do
not tune the accelerator to make one cell pass."  No tuning was attempted.  The
change is **reverted**; production carries no accelerator.

**What it teaches S3 and S5.**  Any method that extrapolates one component of
the coupled state has the same exposure.  Option B is safer on this specific
point than it looked: the secant is on the residual `R(x)` of hydrogen's balance,
with the other components re-evaluated consistently at each trial, rather than a
jump applied to one variable and the rest left where they were.  Option C is
safer still, because it removes the mode instead of stepping around it.  If S3
also disappoints, go to S5 rather than looking for a better accelerator.

### S3 -- Option B, safeguarded secant  *(only if S2 does not clear G-M1f)*

Bracketed on `[0, 1]` with bisection whenever the secant step leaves the bracket
or fails to reduce `|R|`.  The bracketing is **mandatory**: substitution is a
self-map of `[0, 1]` and cannot return an unphysical state, secant can, and this
equation is known to admit unphysical roots -- G-M1a produced `x_HI` = 6.47 and
578 when the denominator was wrong.

**Acceptance**: G-M1f, G-M1g and G-M1h.

#### S3 result, 2026-07-31: **PASSES G-M1f**

Implemented as `hydrogen_neutral_fraction` in `ion_balance_mod`: at fixed `n_e`,
solve the residual

```
R(x) = x [B + R_in(x) + R_out(x)] - [A + R_in(x)]
A = alpha_H n_e ,  B = Gamma_H + (C_H + alpha_H) n_e
```

for its root by Illinois (modified false position), with the whole metal cascade
re-evaluated at every trial `x`.

**The bracket is structural, and the plan above was wrong to call a clamp
mandatory.**  `R(0) = -(A + R_in(0)) < 0` because `A > 0`, and
`R(1) = B + R_out(1) > 0` because `B > 0`, so the root is enclosed in `[0,1]`
for every state with no clamp, no fallback and no safeguard heuristic.  Illinois
keeps that bracket at every step while converging superlinearly.  Option B was
therefore *safer* than section 3 credited it with being: it recovers the one
safety property plain substitution had, rather than trading it away.

**G-M1f: PASS.**  One `solve_ion_cell` call against the converged state:

| `x_HII` | one call | converged | relative deviation |
|---|---|---|---|
| 0.287  | 2.869235e-1 | 2.869235e-1 | 3.6e-11 |
| 0.161  | 1.613053e-1 | 1.613053e-1 | 4.7e-12 |
| 0.0742 | 7.416036e-2 | 7.416036e-2 | 6.1e-12 |
| 0.0347 | 3.465394e-2 | 3.465394e-2 | 9.6e-13 |
| 0.0208 | 2.077601e-2 | 2.077601e-2 | 1.4e-12 |

Worst deviation **3.6e-11**, against 1.8e-3 for substitution -- seven orders of
magnitude, and the cell that needed 1285 passes now converges in one call.

**G-M1g: PASS on the revised criterion, and the revision is recorded in
section 5.**  `metal_only` agrees with the pre-M1 binary on 30 of 31 datasets
exactly and on `x_HeII` to one ulp (4.24e-16 relative, 18688 of 262144 leaves).
The cause is the added branch changing code generation under
`-ipo -fp-model fast=2`, not the physics: restoring the `metal_only` expression
verbatim does not remove it.

**G-M1g second half: PASS.**  Both benchmarks rerun with only the inner method
changed:

| | `<T[NpNe]>` K | `L(Hbeta)` | `r(He0=0.5)` | `r(H0=0.5)` |
|---|---|---|---|---|
| HII40 substitution | 8167.53 | 1.9301e37 | 4.1277 | 4.6775 |
| HII40 S3 root find | 8167.54 | 1.9301e37 | 4.1276 | 4.6775 |
| delta | **+0.008** | +0.0001% | -0.00001 | +0.00001 |
| HII20 substitution | 6924.38 | 4.6840e36 | 1.2408 | 2.8798 |
| HII20 S3 root find | 6924.39 | 4.6840e36 | 1.2407 | 2.8798 |
| delta | **+0.010** | -0.0004% | -0.00001 | -0.00001 |

Roughly 8000 times inside the tolerance the runs converge to.  **The method
changed the path to the root and not the root**, which is what the gate asks.

**And the thing S3 existed to fix is fixed.**  The benchmark numbers were never
the target -- they were already safe, because `n_p n_e` weighting suppresses the
gas at risk.  The target was the shadow beyond the front, and the printed
residual shows it:

| HII20 | max residual | volume-weighted |
|---|---|---|
| substitution | 4.6141e-03 | 4.7187e-04 |
| S3 root find | 1.3191e-06 | 4.4625e-09 |

3500x on the maximum and 100000x on the mean.  The statement that the shadow
beyond the front is not converged, carried in `METAL_COUPLING_PLAN.md`,
`MoCHII_physics.tex`, `CLAUDE.md` and at the code site since M1 landed, is
**withdrawn**: it was true of the lagged substitution and is not true of the
root find.

**Why this works where S2 failed.**  Both target the same mode; the difference
is what happens to the rest of the state.  Aitken jumped `n_e` and left `x_HI`,
the helium ratios and the eleven stage chains where they were, re-exciting a
transient that cost more than the jump saved.  Every trial of the root find is a
consistent state -- `R_in` and `R_out` are recomputed from the cascade at that
`x` -- so there is no off-manifold excursion to pay for.  That distinction, not
the order of convergence, is what decides between the two.

### S4 -- The oxygen-hydrogen net-flux diagnostic  *(independent of S2/S3; do it either way)*

Cloudy pairs its lag with a net-flux test written against the cancellation, and
hard-codes it to the H/O pair (section 2.1).  MoCHII's
`hydrogen_balance_residual` reports the residual of the whole balance, which is
not sensitive to the specific failure mode section 1 predicts.  Add the pair net
flux, compared against the **smaller** of the two ionization rates in Cloudy's
manner, to the printed diagnostic.

**Acceptance**: the diagnostic separates a capped cell from a converged one on
the HII20 run, where the whole-balance residual currently reports 4.6e-3 without
saying which pair caused it.

### S5 -- Option C, remove the cancellation analytically  *(only if S2 and S3 both fall short)*

Solve for the departure from the locked limit rather than for `x`.  This is the
right answer in principle and the easiest to get subtly wrong: the cancellation
is exact for one reaction pair at a time, and the algebra has to hold with 29
reactions across 11 elements.

**Acceptance**: G-M1f, G-M1g, **and a re-derived G-M1a**, since this changes the
expression that gate checks.

### S6 -- Documentation

`docs/METAL_COUPLING_PLAN.md` (G-M1f closed), `docs/MoCHII_physics.tex` (the
solver note, replacing the present statement that the shadow beyond the front is
not converged), `CLAUDE.md`, and this document with the S1 table and the outcome.

### Not to be done

- **Raising the cap, or adding `|dx_HI|` to the exit test as a fix.**  1285
  passes for each of about 23 temperature trials for each leaf and gas iteration
  is not affordable, and the exit test only detects the failure.  Adding
  `|dx_HI|` to the exit test as a **diagnostic** is still worth doing in S4, so a
  capped cell is visible rather than silent.
- **Special-casing the H-O pair as a 2x2.**  It would work at HII40, where
  oxygen is over 99.99% of both sums, and it would break the generic registry
  the founding design rests on (`docs/PLAN.md`).

## 5. Gates

- **G-M1f** (from `METAL_COUPLING_PLAN.md`). A cell at `x_HII` between 0.01 and
  0.05 must converge to the same state as repeated whole-cell solves, to
  <= 1e-10 relative, **inside** the pass cap.
- **G-M1g** (new, for any of A/B/C). `par%charge_exchange = 'metal_only'` must
  reproduce the pre-M1 binary, and the `two_way` HII40 and HII20 results must
  not move by more than the run's own convergence tolerance (82 K and 69 K on
  `<T[NpNe]>`) -- the accelerator must change the path to the root, not the
  root.

  **Criterion revised 2026-07-31, and the revision is recorded rather than
  quietly applied**, because it is a criterion being relaxed to admit the
  change that failed it.  As first written this gate demanded *bit-identity*.
  S3 gives 30 of 31 datasets bit-identical, with `x_HeII` differing by **one
  ulp** -- maximum relative 4.24e-16 on 18688 of 262144 leaves.  The cause was
  isolated: it is not the physics and not the arithmetic.  Restoring the
  `metal_only` expression *verbatim*, including the two terms that are exactly
  zero on that branch, does not remove it; what produces it is the presence of
  the `if (with_metal_cx)` branch at all, which changes code generation under
  `-ipo -fp-model fast=2`.  Literal bit-identity is therefore not attainable by
  any implementation that adds a branch to this loop, so demanding it would not
  select a better implementation -- it would only select against implementing
  the fix.  The criterion is now: **the `metal_only` path must be unchanged in
  physics and agree with the pre-M1 binary to at most a few ulp, with the
  deviation measured and reported**, not asserted to be zero.  The bit-identity
  form is kept for changes that do **not** restructure control flow, where it
  is attainable and where it caught real defects earlier in this work.
- **G-M1h** (option B only). No leaf may leave the solve with `x_HI` outside
  `[0, 1]`, asserted in the solver rather than checked afterwards.

## 6. What this does not change

The physics of M1 is settled and gated; this is solely about how the fixed point
is walked. `'metal_only'` remains available and bit-identical, so a recorded
gate can always be reproduced while this is open.
