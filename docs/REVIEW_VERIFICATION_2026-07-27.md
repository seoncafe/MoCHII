# Verification of MOCHII_ALGORITHM_PHYSICS_REVIEW.md, and the resulting fix plan

- Date: 2026-07-27
- Subject: `MOCHII_ALGORITHM_PHYSICS_REVIEW.md` (external review, same date)
- Method: every cited line range was opened and compared against the review's
  claim; the two build claims were reproduced by compiling; the test-script
  claim was checked by inspecting exit paths.
- Compilers used for the build checks: GNU Fortran 13.1.0, Intel MPI `mpif90`
  (which wraps `gfortran` on this machine), `mpiifort` (Intel oneAPI 2025.3).

### Revision history

**Rev. 1 (initial)** — verification only, no source file modified. Sections 1-3
record what the review got right, what its priorities should be, and what it
missed; section 5 is the resulting plan.

**Rev. 2 (2026-07-27, this document)** — five findings implemented and
measured. Section 4a is new and carries the measured effect of each; the plan
in section 5 is annotated with what is DONE and what is still open. Nothing
was committed.

| Finding | Status | Where measured |
|---|---|---|
| 3.4 C II/O I two-level denominator | **DONE** | §4a, +0.31% `<Te>` here, factor 3 in the high-T limit |
| 3.5 `metal_cooling_H` gated by `grain_pe` | **DONE** | §4a, -69.7% `<Te>` with `grain_pe = .false.` |
| 3.3 `laursen09_live` loses `opac_norm` | **DONE** | §4a, dust heating 79x |
| 6.10 working-directory-relative data paths | **DONE** | §4a, was destroying output silently |
| 6.8 rank-0-gated diagnostic before a collective abort | **DONE (6 of 55 sites)** | §4a |
| 3.1, 3.2, 3.6, 3.7, 3.8, 6.1, 6.2, 6.3, 6.4, 6.6, 6.7 | open | §5 |

**Rev. 3 (2026-07-28)** — impact of the Rev. 2 fixes measured against the
existing test suite (new section 4b) and the first test gates added or made
enforceable (new section 4c). The MOCASSIN HII20/HII40 gate metric does not
move; HII40 was already failing its 5% gate beforehand. Still nothing
committed. A separate performance plan,
[THERMAL_SOLVE_PARALLEL_PLAN.md](THERMAL_SOLVE_PARALLEL_PLAN.md), was written
after the timing measurements taken during this work; it is not a review
finding and is sequenced after this work.

**Rev. 4 (2026-07-28)** — work resumed from GitHub commit `4624f08` after an
incomplete edit of `hdf5io_mod.f90` was backed up outside the repository.
Plan items 6 and 8 are now implemented, and the case-sensitive half of item 9
is fixed. HDF5 operations preserve the first failure across cleanup, output
helpers share one status instead of resetting it, and rank-0 output routines
abort instead of printing success after an I/O failure. The HDF5-off branch
now contains ordinary Fortran procedures rather than CPP-generated stubs.
OpenMP flags are selected per compiler. A focused `tests/io/check_io.sh`
regression checks HDF5-off compilation, normal HDF5 output and readback,
failed-create propagation, first-error preservation, and uppercase
`.H5`/`.HDF5` detection.

Verification for Rev. 4:

| Check | Result |
|---|---|
| Intel debug full build, `HDF5=1` | **PASS** |
| Intel debug full build, `HDF5=0` | **PASS** |
| GNU default-line-length compile of HDF5-off backend | **PASS** |
| `tests/io/check_io.sh` | **PASS** |
| `tests/two_level/check_two_level.py` | **PASS** |
| GNU full link | Not claimed: the bundled SEDust archive and module files are Intel-built and must be rebuilt with GNU first |

Files changed in Rev. 2: `species_mod.f90`, `thermal_mod.f90`,
`physics_amr_mod.f90`, `grid_mod_amr.f90`, `gas_opacity_mod.f90`,
`utility.f90`, `define.f90`, `setup.f90`, `nebcont_mod.f90`,
`gaunt_vh14_mod.f90`, `sh95_mod.f90`, `cooling_mod.f90`, `dust_temp_mod.f90`.

## 1. Verdict

No factual claim in the review was found to be wrong. Every cited line range
matched within a line or two. The problems with the document are not accuracy
but **priority calibration** and **three missing items**: findings whose real
effect is a factor of 1.8 sit in the same priority band as findings whose real
effect is 0.4%.

### 1.1 Physics findings

| Review | Claim verified? | Effect size (this verification) | Priority as written | Priority after verification | Status |
|---|---|---|---|---|---|
| 3.4 C II/O I two-level denominator | Yes | Factor 1.80 (C II, 100 K) rising to 3; 1.6 (O I, high T). **Active in the density range the term was written for** | P0 | **P0, first** | **DONE** |
| 3.3 `laursen09_live` loses `opac_norm` | Yes | Unbounded opacity jump at the first refill | P0 | P0 | **DONE** |
| 3.5 `metal_cooling_H` gated by `grain_pe` | Yes | An independent coolant disappears with a heating switch | P1 | **P0 after measurement** | **DONE** |
| 3.2 AMR remap not conservative | Yes, and reachable | Coarsening only; default `refine_lbase = 5` will coarsen deeper input leaves | P0 | P0, coarsening only | open |
| 3.1 Final state and rates one iteration apart | Yes | Bounded by `gas_tol` when converged; unbounded only at the iteration cap | P0 | P0 unconverged / P1 converged | open |
| 3.7 Metal stage product overflow | Yes, as a latent risk | Not reachable with physical `n_e`; failure mode misstated (see 3.3 below) | P1 | P2 | open |
| 3.8 Single monotonic thermal root | Yes | Real limitation; the silent boundary clamp is the part with no diagnostic | P1 | P2 | open |
| 3.6 Diffuse packet count and spectrum | Yes | Rounding ~1e-6 relative; energy weighting 0.36% (H) and 1.4% (He I 2-photon); `eion_max` clip `~e^{-100}` | P1 | P2 | open |

### 1.2 Software findings

| Review | Verified how | Result |
|---|---|---|
| 6.2 GNU `HDF5=0` build broken | Compiled | **FIXED in Rev. 4.** The non-HDF5 backend now uses ordinary Fortran procedures and compiles under GNU's default free-form line limit |
| 6.3 `-qopenmp` on the GNU link | Compiled | **FIXED in Rev. 4.** Intel selects `-qopenmp`; `F90=mpif90` selects `-fopenmp` |
| 6.1 HDF5 status masked | Read and executed | **FIXED in Rev. 4.** Primary errors survive cleanup, nested writers propagate one status, and output routines check it before reporting success |
| 6.7 Test scripts do not fail the process | Read | Confirmed. `check_g1_stromgren.py`, `check_g2a.py`, `check_g4.py`, `check_pdr.py` print `GATE: PASS/FAIL` and exit 0 either way |
| 6.10 working-directory-relative data paths | **Reproduced, and worse than written** | The dense PDR case lost its entire output to this. See §4a. **FIXED in Rev. 2** |
| 6.10 `io_detect_format` case sensitivity | Read and executed | **FIXED in Rev. 4.** `.h5` and `.hdf5` matching is case-insensitive and covered by `tests/io/check_io.sh` |
| 6.4 Build portability | Read | Accurate. Adding that `OBJS` happens to be in topological order, so serial builds work and only `make -j` is unsafe |
| 6.6 Input validation | Read | **Overstated.** `setup.f90` already contains 47 validation checks with explicit `ERROR:` messages. The specific gaps named are real (see §3) but the summary reads as if validation were absent |
| 6.8 MPI error handling | **Reproduced. This verdict is revised** | Rev. 1 called it "standard practice recommendations, not defects". That was wrong: a rank-0-gated diagnostic in front of a collective abort produced a completely silent failure that also destroyed the output file. See §4a. **FIXED in Rev. 2 at the 6 file-open sites**; 49 setup-validation sites remain |
| 6.9 Global state and coupling | Read | Plausible; standard practice recommendation, not a defect |

## 2. The one finding whose priority must move

Section 3.4 is the highest-value item and the review understates it.

`src/species_mod.f90:502-513` computes

```text
Lambda = n_ion n_HI q_lu dE / (1 + n_HI q_ul / A)
```

The two-level steady state with `n_ion = n_l + n_u` gives

```text
n_u (A + n_HI q_ul) = n_l n_HI q_lu
Lambda = n_u A dE = n_ion n_HI q_lu dE / (1 + n_HI (q_ul + q_lu) / A)
```

so the upward rate is missing and the dense limit is overestimated by
`1 + q_lu/q_ul = 1 + (g_u/g_l) exp(-dE/kT)`:

- C II: 1.803 at 100 K, approaching 3 at high T (`g_u/g_l = 2`).
- O I: 1.062 at 100 K, approaching 1.6 at high T (`g_u/g_l = 0.6`).

The reason this is first rather than fourth: the H critical density is
`n_crit = A/q_ul = 3.05e3 cm^-3` for C II at 100 K, and the code comment at
`species_mod.f90:476-478` states the term was added for `nH ~ 5e3` PDR gas.
The error is therefore not latent — it is at or near full strength in exactly
the regime the term exists to model.

## 3. Corrections to the review text

1. **3.2** — state that the defect is coarsening only. Point sampling at the
   new leaf center is exact piecewise-constant prolongation under pure
   subdivision, so `sum(nH V)` and the rest are conserved there. The reachable
   trigger is `refine_lbase = 5` (default) against an input octree with leaves
   at level 6 or 7 away from the front.
2. **3.2** — the `species_resize` sub-point should be downgraded. Zeroing
   `egam`/`eheat` leaves the first new-grid metal opacity wrong, but
   `main.f90:201` sets `converged = .false.` after a refinement event, so the
   next transport rebuilds the rates. It is a one-iteration transient, not a
   persistent state error.
3. **3.7** — the failure mode is misstated. With `prod` overflowing to `Inf`,
   `s = Inf` gives `frac(1) = 1/Inf = 0`, and every later stage is
   `0 * finite = 0`. The result is **all stage fractions zero with a sum of
   zero**, which is silent, not `NaN`. A `NaN` needs `r(1) = Inf` as well,
   which needs `rrec` below `tinest` together with `rion > ~0.4 s^-1`. Silent
   zeroing is the more likely and the worse outcome, and it should be the
   headline of that finding.
4. **3.6** — add the measured sizes so the priority is defensible:
   relative luminosity error `0.5/N_diff` (order 1e-6 for any production
   packet count); mean packet energy low by 0.36% for the H ground continuum
   at 1e4 K (`<E>` = 14.460 eV sampled versus `<E^2>/<E>` = 14.512 eV) and by
   1.4% for the flat He I two-photon branch (17.109 versus 17.354 eV); the
   `eion_max` clip is suppressed by `exp(-100)` at the 100 eV default.
5. **6.3** — record that the Makefile's compiler detection selects `mpiifort`
   on this machine, so the broken flag is on the `F90=mpif90` path only.
6. **6.6** — replace "lack complete checked parsing" with the enumerated gaps,
   and state that 47 checks already exist. The real gaps are: `open`/`read` of
   the namelist without `iostat`/`iomsg` (`setup.f90:34-36`); `conv_crit` not
   validated against an enumeration, so anything other than `'vol'` silently
   selects the cell criterion at `main.f90:167`; no `refine_lbase <=
   refine_lmax` check; `no_photons` declared `real(wp)` and truncated into an
   `int64`; every rank reading the namelist independently.

## 4. Items the review missed

1. **`par%taumax` and `par%tauhomo` are overwritten with the realized values**
   at `grid_mod_amr.f90:335-336`, after normalization. Any later attempt to
   recover `opac_norm` from `par` is therefore circular. This constrains how
   3.3 can be fixed and the review's recommendation does not mention it.
2. **The electron and hydrogen colliders are treated as two independent
   two-level systems**, each assuming the full `n_ion` in the lower level:
   `metal_cooling` uses the electron n-level table, `metal_cooling_H` adds a
   separate saturated two-level term. Where both colliders matter the ion
   population is counted twice. The correct form has one denominator,
   `1 + [n_HI(q_ul^H + q_lu^H) + n_e(q_ul^e + q_lu^e)]/A`. In the PDR
   (`n_HI >> n_e`) the split is harmless, which is why it has not shown up,
   but it is the deeper form of finding 3.4 and the same fix addresses both.
3. **The diffuse ground continuum drops the Milne weighting.** The true
   photon-number spectrum is `nu^2 sigma(nu) exp(-h(nu-nu_0)/kT)`; the code
   samples `exp(-(E-E_th)/kT)`, which assumes `nu^2 sigma(nu)` is constant
   across the emitted range. The module header at `diffuse_mod.f90:13-14`
   calls it "the near-threshold Milne shape" without stating that the
   `nu^2 sigma` factor was dropped or over what range that holds. Under the
   project rule for deliberate approximations, the domain of validity belongs
   at the code site.

## 4a. Measured effect of the fixes that were applied

Items 1, 2, 3 and the two infrastructure findings 6.8/6.10 were implemented and
exercised. Test case: `tests/pdr/pdr_dense.in` (Cartesian 48^3, `nH = 5e3`,
Draine ISRF, no ionizing source, `use_metals`, `solve_te`, `te_min = 10 K`,
30 iterations, 8 MPI ranks), plus `tests/d_dusty/d_dust_l09.in` with
`taumax = 5` added for the dust normalization.

**Item 1, C II / O I two-level denominator.** Analytic limits first: after the
fix the dense limit equals `A * n_u/n_ion` with `n_u/n_ion` the Boltzmann
population to 6 digits, and the low-density limit is unchanged. Before the fix
the same limit exceeded the Boltzmann value by `1 + (g_u/g_l)e^{-dE/kT}` —
for C II at 10^4 K that is a factor 2.98, i.e. an implied upper-level
occupation of 198% of the ion population, which is impossible. That is the
decisive evidence.

The effect on this particular test is small, and smaller than a naive estimate
suggests. **The correction factor is strongly temperature dependent and this
PDR is cold**: the gas sits at 28-37 K, where `q_lu/q_ul = 0.077` and the
saturated cooling changes by only 4.4%, against 50% at 100 K and a factor 3 in
the high-temperature limit. Measured over the 42168 gas leaves:

```text
volume-weighted <Te>   27.586 -> 27.672 K   (+0.31%)
max |dTe|/Te           +1.41%
leaves with |dTe|/Te > 0.1% / 1%      69.0% / 9.0%
```

The sign is right (less cooling, warmer gas). A warmer PDR would show much
more; this case understates the fix.

**Item 2, `metal_cooling_H` no longer gated by `grain_pe`.** With
`grain_pe = .true.` the term merely moved, and the run reproduces the
item-1-only run (identical convergence metrics at iteration 30). With
`grain_pe = .false.` the old code never called the H channel at all, so the
pair isolates the fix exactly:

```text
volume-weighted <Te>   48.729 -> 14.788 K   (-69.7%)
dTe/Te                 -27.8% to -84.7%, every gas leaf affected
```

Metal photoheating was running without its dominant coolant. The coldest cells
reach 11.0 K against a `te_min` of 10 K, so they are near but not at the floor.

**Item 3, live-dust normalization.** With `laursen09_live` and `taumax = 5`
(`te_fixed`, metals off, so items 1 and 2 are inert):

```text
                   before        after      ratio
Heat_dust sum    6.69e-17     5.32e-15      79.5
dust IR L        4.35e+25     2.14e+27      49.2
J_nu sum         8.86e-12     1.02e-12      0.115
<x_HI>              0.366        0.514
```

Before the fix the requested optical depth vanished at the first refill.
Afterward the attenuation, dust heating and neutral fraction are all consistent
with the requested dust scale.

**Findings 6.8 and 6.10.** These were fixed because they blocked the work
above, and the failure they produce is worse than the review states. Running
the dense PDR case with the binary outside the source tree, `gammaHI.dat` was
looked up relative to the current working directory, failed, and non-zero ranks
called `MPI_ABORT` while rank 0 was still inside `gas_rates_write`. The
diagnostic was gated on rank 0, so **nothing at all was printed**, and the
partially written rates file was destroyed — after 30 iterations of transport
had already been paid for. A silent loss of output after a full run is a more
severe outcome than "behavior depends on the working directory".

Changes: `par%data_dir` resolved at setup to `<directory of MoCHII.x>/data`
(explicitly settable), used by `nebcont_mod` and `gaunt_vh14_mod`;
`par%atomic_dir` defaults to `<data_dir>/atomic` instead of the working-
directory-relative `data/atomic`; and `utility::fatal_error` writes a
rank-qualified message from the rank that actually failed before aborting,
now used at the six "cannot open" sites. For the normal layout (binary in the
source root, run from `tests/<case>/`) the resolved paths are unchanged; only
the dependence on the working directory is gone.

Not done, reported only: 49 of the 55 `MPI_ABORT` sites remain rank-0-gated.
Those are setup-time input validation where every rank fails identically, so
rank 0 does print; converting them is mechanical follow-up work.

Verification cost note: `define.f90` and `utility.f90` both changed, and the
incremental `make MoCHII.x` linked without complaint despite stale `.mod`
files elsewhere. A full rebuild was required. That is finding 6.4 biting in
practice.

## 4b. Impact of the applied fixes on the existing test suite

Item 2 turns on `metal_cooling_H` in every case that has `use_metals` without
`grain_pe` — **38 of the shipped test inputs**, including the MOCASSIN HII20
and HII40 benchmarks, which are this code's external calibration. Those
benchmarks were therefore run before and after, at 8 MPI ranks, with
`tests/g2_hii/hii20.in` and `hii40.in` unchanged apart from absolute paths.

### The gate metric does not move

`check_g2_hii.py` gates on the ion-weighted `Te(H II)` against the MOCASSIN
baseline:

| case | before | after | verdict |
|---|---|---|---|
| HII20 | 6922.4 K (dev -4.50%) | 6922.1 K (dev -4.51%) | PASS -> PASS |
| HII40 | 7929.0 K (dev -5.65%) | 7928.8 K (dev -5.65%) | FAIL -> FAIL |

The shift is -0.004%. **HII40 already failed the 5% gate before any of this
work**, at the same -5.65%; that is a pre-existing state, not a regression
introduced here, and it is recorded so it is not later attributed to these
fixes.

### The neutral-zone metric moves as the physics requires

`Te(H I)` is weighted by `n_HI n_e V`, so it samples exactly the region where
the H-impact channel acts:

| case | before | after | absolute |
|---|---|---|---|
| HII20 | dev +0.40% | dev **-0.24%** | 7074.8 -> 7029.7 K |
| HII40 | dev -9.30% | dev -9.69% | 8149.2 -> 8114.3 K |

HII20 moves closer to MOCASSIN, HII40 slightly further. Individual neutral
leaves cool by up to 26%, the median leaf does not move, and the change is
**one-sided**: `max dTe = 0`, no leaf gets warmer. That is the correct
signature — with `grain_pe` off the old code called the H channel nowhere, so
only the added coolant is at work.

Mean ionic fractions move in the fourth decimal at most (`C II` 0.9645 ->
0.9646, `S II` 0.6017 -> 0.6015, and similar for N, O, Ne, He). No benchmark
conclusion depends on them.

**Verdict: safe for the external benchmarks.** `metal_cooling_H` scales with
`n_HI` and vanishes in ionized gas, which is what the numbers show.

### Side finding: the benchmark inputs do not converge

All four runs stopped at the 100-iteration cap with
`max |dx_HII|` = 3.5e-2 (HII20) and 4.4e-2 (HII40) against `gas_tol = 5e-3`.
`hii20.in` and `hii40.in` set no `conv_crit`, so the default `'cell'`
(cell-maximum) applies, which these cases do not reach. The `_bench` variants
set `conv_crit = 'vol'` with `gas_tol_te = 1e-2` and converge in 27 and 38
iterations.

This predates the present work, but it means **the MOCASSIN comparison is
being made on a state that is not at equilibrium**, and finding 3.1 applies to
it in full: at the iteration cap with `require_convergence = .false.` the
written gas state and the written rates are one iteration apart with no bound
on the mismatch. It is an additional argument for the final consistency pass
(plan item 5), and separately the gate inputs should either adopt the `_bench`
convergence settings or raise the cap.

## 4c. Test gates added and converted

- **New `tests/two_level/check_two_level.py`** — gate for finding 3.4 that
  needs no Monte Carlo run. It reads the denominator form out of
  `species_mod.f90`, then checks both density limits at six temperatures
  (low density `Lambda/(n q_lu) = 1`; LTE `Lambda/(A x_u^Boltzmann) = 1`) and
  rejects an implied upper-level occupation above the ion population.
  Confirmed to catch the defect: reinstating the old C II denominator makes it
  exit 1.
- **Exit codes added** to `check_g4.py`, `check_g2a.py`, `check_peel.py` and
  `check_g0_gamma.py`; all four verified to exit 0 on their current data.
  `check_g4.py` additionally treats a missing input file as a failure — it
  previously printed "skipped" and exited 0 with the gate never evaluated.
- **Still open**: 15 of 25 `check_*.py` scripts have no process exit. Six of
  them (`check_pdr.py`, `check_g1_stromgren.py`, `check_lines.py`,
  `check_fuv.py`, `check_sedust.py`, `check_v3.py`) print numbers with **no
  pass/fail criterion at all**. Defining those criteria is a physics judgment
  and was deliberately not done unilaterally.

## 5. Fix plan

The review's own section 8 puts build and I/O first. This plan puts the active
physics errors first, because physical correctness is the acceptance criterion
and items 1 and 2 change current PDR results.

### Decisions taken

- **`laursen09_live` combined with `taumax`/`tauhomo` keeps the initial
  normalization and then evolves freely.** `opac_norm` is computed once, from
  the initial gas state, and is held as persistent model state; every later
  `gas_opacity_fill` multiplies by that stored constant. It is never
  recomputed. `taumax`/`tauhomo` therefore act as a one-time calibration of
  the dust scale — replacing the absolute `DGR * cext_dust` normalization —
  and after that the optical depth moves only through the ionization-dependent
  survival factor `f_ion_dust`/`f_ion_pah`. There is deliberately **no**
  re-normalization at each iteration: pinning the realized optical depth to
  the target every iteration would defeat the live dust model.

  The semantics to document at the code site: `taumax`/`tauhomo` is the
  optical depth of the **initial** configuration (whatever `xHI_init` or the
  input `xHI` column gives), not of the converged one. The realized value will
  drift away from the requested one as the gas ionizes, and that drift is the
  physical result, not an error.

  This differs from the current code in one respect only, but a decisive one:
  the constant must survive the refill. Dropping it, as
  `gas_opacity_mod.f90:201-213` does now, makes the initial setting
  meaningless after a single iteration and produces a discontinuous global
  opacity jump.

### Order

1. **[DONE, Rev. 2] C II/O I denominator** (`species_mod.f90:502-513`). Use
   `1 + nHI*(qul + qlu)/A`. Two lines. Limits checked analytically (§4a); the
   low-density scaling is unchanged and the dense limit now equals the
   Boltzmann population. **Still open: a regression test in `tests/` asserting
   both limits — the check was run once by hand, not committed as a gate.**
2. **[DONE, Rev. 2] Move `metal_cooling_H` out of the `grain_pe` block** into
   the `use_metals` cooling path. Measured in §4a: with `grain_pe = .false.`
   the coolant was absent entirely and the gas ran 70% too warm. **Still open:
   a committed test that toggling `grain_pe` changes the photoelectric heating
   but leaves the H-impact cooling in place.**
3. **[DONE, Rev. 2] Persist the live-dust normalization.** Store `opac_norm` as module state
   (it is currently a local at `grid_mod_amr.f90:41`), broadcast it, and apply
   it in the `laursen09_live` branch of `gas_opacity_fill`. Compute it once,
   at grid setup, and never again. Because `par%taumax`/`par%tauhomo` are
   overwritten with the realized values at `grid_mod_amr.f90:335-336`, the
   factor cannot be recovered from `par` afterward and must be kept in its own
   variable — see §4.1. Report the realized pole or volume-average optical
   depth each iteration so the drift away from the initial target is visible.
   Test: run `laursen09_live` with `taumax` and again with `tauhomo`, and
   require that the first refill changes the opacity by the survival factor
   alone, with no jump of order `1/opac_norm`.

   **[DONE, Rev. 2]** implemented as `physics_amr_mod::dust_opac_norm`, set in
   `grid_mod_amr` after the broadcast and applied in the `laursen09_live`
   branch of `gas_opacity_fill`. Verified with `taumax = 5` (§4a). **Still
   open: the `tauhomo` variant, the per-iteration realized-tau report, and a
   committed test.**
4. **AMR coarsening**: short term, refuse or warn when `refine_lbase` is
   coarser than the deepest input leaf, and validate `refine_lbase <=
   refine_lmax`. Longer term, remap extensive quantities through old/new leaf
   volume overlaps and assert conservation of `sum(nH V)`, each ion stage
   count, `sum(ne V)`, `sum(rhokap V)`, and thermal content before and after.
5. **Unconditional final transport-and-rates pass** after the nonlinear loop
   in `main.f90`, with no gas solve, so the written rates, dust heating, and
   dust temperature belong to the written gas state. Cost is one extra
   transport, that is `1/niter` of the run.
6. **HDF5 status preservation**: keep the first error in a dedicated variable,
   never let a cleanup call overwrite it, and return a checked status from
   every output routine so a failed write cannot print a success line.

   **[DONE, Rev. 4]** `hdf5io_mod` sequences operations, records the first
   backend failure, and uses separate cleanup statuses. Nested output helpers
   no longer reset status, and every HDF5-producing rank-0 output path checks
   the final status before printing success. Covered by `tests/io/check_io.sh`.
7. **Make the regression suite enforceable**: give every `check_*.py` a
   nonzero exit on a failed gate or a missing product, and add one top-level
   runner. This has to land early — without it the changes above cannot be
   validated safely.
8. **Build**: replace the preprocessor-generated stubs in the non-HDF5 branch
   of `hdf5io_mod.f90` with ordinary procedures, and select the OpenMP flag by
   compiler.

   **[DONE, Rev. 4]** Both Intel HDF5-on and HDF5-off debug builds link. The
   HDF5-off module also compiles with GNU under its default line-length limit,
   and the GNU Makefile path now uses `-fopenmp`. A full GNU application link
   still requires rebuilding the bundled Intel SEDust archive and `.mod`
   files with GNU.
9. **Low-cost remainder**: give each diffuse packet `diffuse_lum/diffuse_nphot`
   instead of the stellar `Lpacket`; make `io_detect_format` case-insensitive;
   report the count of leaves pinned at `te_min`/`te_max` by the thermal
   solver; add `iostat`/`iomsg` to the namelist read; validate the `conv_crit`
   enumeration.

   **[PARTIAL, Rev. 4]** `io_detect_format` is now case-insensitive and tested.
   The other low-cost items remain open.
10. **Accuracy improvements with small current effect**: energy-luminosity CDF
    sampling inside each diffuse channel; log-sum-exp metal stage fractions;
    multi-bracket thermal root scan; H colliders folded into the n-level
    solver so the two-level and n-level paths cannot diverge (this also
    removes the double counting in §4.2).
11. **Structural work**, only once 7 is green: validated central parsing, MPI
    error checking, build-system replacement, splitting `define.f90`.

### Items pulled forward in Rev. 2 (out of the original order)

12. **[DONE, Rev. 2] Executable-relative data root** (was part of item 10).
    `par%data_dir`, resolved at setup, replaces the working-directory-relative
    `'../../data/x'` / `'data/x'` lookups in `nebcont_mod` and
    `gaunt_vh14_mod`; `par%atomic_dir` now defaults to `<data_dir>/atomic`.
    Pulled forward because it was silently destroying run output (§4a).
13. **[DONE, Rev. 2, partial] Rank-qualified fatal abort** (was part of item
    10). `utility::fatal_error` prints from the rank that failed, then aborts.
    Applied at the six "cannot open" sites. **Still open: the other 49
    `MPI_ABORT` sites**, all setup-time input validation where every rank
    fails identically, so rank 0 does print. Also still open: `mpi_check` on
    MPI return codes, one collective fatal path, and `mpi_f08`.
