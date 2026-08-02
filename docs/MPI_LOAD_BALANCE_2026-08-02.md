# MPI load balance of the ionizing-band transport

Date: 2026-08-02.  Status: measured, then implemented and re-measured.  The
dynamic schedule now exists (`src/photon_schedule_mod.f90`) and is off by
default (`par%use_master_slave = .false.`); see "What was implemented" below.

## Why this was measured

Until the change described below, MoCHII split the photon loop **cyclically**
over MPI ranks, in the driver itself:

```fortran
do ip = mpar%p_rank+1, par%nphotons, mpar%nproc      ! src/main.f90, before the change
```

so every rank carries the same photon *count* to within one packet.  Whether it
carries the same *work* is a different question, and a global reduction
(`jtally_ion_reduce`) follows immediately, so the pass costs the time of its
**slowest** rank.

MoCafe v2.00 answers this with a master-slave scheduler (`run_master_slave`,
`use_master_slave = .true.` by default): rank 0 hands out batches of
`num_send_at_once` photons and receives with `MPI_ANY_SOURCE`, so a rank that
finishes early immediately gets more.  When this was measured MoCHII had not
ported it: `use_master_slave`, `num_send_at_once` and `use_reduced_wgt` were all
declared and read nowhere (see `MOCAFE_PARAMETER_DIVERGENCE.md`).  The first two
are read now; `use_reduced_wgt` is still read nowhere.

The question is therefore not whether dynamic scheduling is good in general — it
plainly is for some transport problems — but **how much idle time there is in
MoCHII's ionizing band to give back**.

## Method

A temporary instrumentation in `src/main.f90` — where the transport pass lived
at the time; it is `src/ionizing_field_mod.f90` now — timed each rank's photon
loop (stellar + diffuse) with `MPI_Wtime`, gathered the 70 values to rank 0 with
`MPI_GATHER`, and printed min / q1 / median / q3 / max, the mean, `max/mean`,
and

```
idle at the reduction = (max - mean) / max
```

which is the fraction of the pass's wall time spent waiting for the slowest
rank — the upper bound on what any scheduler could recover.  The stellar loop
was timed separately from the pass total so the diffuse packets could be
separated out.

The patch is kept at `scratchpad/load_instrument.patch` for the session; it is
54 added lines and nothing else, and it was reverted before anything below was
written, so the numbers here and the restructuring below are independent.

Runs: **70 MPI ranks** on a 72-core node, `4.0e7` photons per pass, the L6
grids, `gas_niter` cut to 3 so each case ends after three iterations plus the
final consistency pass.  Cutting the iteration count does not change the work
per pass, which is what is being measured.

| case | input | what it isolates |
|---|---|---|
| A_pure | `tests/g1_stromgren/g1_stromgren.in` | path-length variance alone: no dust, no metals, no diffuse field |
| B_bench | `tests/g2_hii/hii40_bench.in` | production shape: metals + diffuse field + thermal balance |
| C_scat | `tests/d_dusty/d_dust_scat.in` | `ion_add_dust` + `ion_dust_scatter`, the only path with a scattering loop |

## Result

| case | pass | mean [s] | min | med | max | max/mean | **idle** |
|---|---|---|---|---|---|---|---|
| A_pure | iter 1 | 1.385 | 1.078 | 1.393 | 1.586 | 1.145 | 12.7% |
| | iter 2 | 1.583 | 1.024 | 1.594 | 1.943 | 1.228 | 18.5% |
| | iter 3 | 1.776 | 1.169 | 1.788 | 2.141 | 1.206 | 17.1% |
| | final | 1.960 | 1.278 | 1.974 | 2.219 | 1.132 | 11.7% |
| B_bench | iter 1 | 40.140 | 25.509 | 38.078 | 47.827 | 1.192 | 16.1% |
| | iter 2 | 40.482 | 25.921 | 38.272 | 48.583 | 1.200 | 16.7% |
| | iter 3 | 39.607 | 25.069 | 37.409 | 47.576 | 1.201 | 16.8% |
| | final | 39.145 | 24.762 | 36.998 | 47.001 | 1.201 | 16.7% |
| C_scat | iter 1 | 18.650 | 12.176 | 17.428 | 21.985 | 1.179 | 15.2% |
| | iter 2 | 15.799 | 10.366 | 14.960 | 18.364 | 1.162 | 14.0% |
| | iter 3 | 10.728 | 7.293 | 10.594 | 11.639 | 1.085 | 7.8% |
| | final | 9.443 | 6.471 | 9.447 | 9.947 | 1.053 | 5.1% |

In B_bench the stellar loop is 71-73% of the pass; the remaining quarter is the
diffuse packets.  A and C have no diffuse field, so the pass is the stellar loop.

## What the distribution looks like, and what causes it

A second run instrumented the same loop with **work counters on every rank** —
the photon count and `direct_segments`, the number of path segments the rank
actually scored — gathered alongside the time.  That separates the two possible
causes completely.  `hii40_bench`, 70 ranks, `4.0e7` photons:

| quantity | spread across the 70 ranks |
|---|---|
| photons per rank | 571428 - 571429 (**0.000%**, equal by construction) |
| path segments per rank | 33 253 313 - 33 310 451 (**0.17%**) |
| wall time per rank | 17.1 - 35.8 s (**109%**) |

**Every rank does the same work and takes up to twice as long to do it.**  The
photon assignment is not the cause; the execution environment is.

The node explains it.  `lscpu` reports 72 logical CPUs as **36 physical cores x
2 hyperthreads** (2 sockets x 18 cores, 4 NUMA nodes with the CPU numbering
interleaved), and `/sys/.../thread_siblings_list` pairs CPU *c* with CPU *c*+36.
Running 70 ranks therefore packs most ranks two-to-a-core.  Re-running the same
case with **36 ranks** — one per physical core — collapses the spread:

| ranks | min | median | max | max/min |
|---|---|---|---|---|
| 70 | 17.1 | 26.1 | 35.8 | **2.09** |
| 36 | 28.6 | 29.1 | 34.9 | **1.22** |

so hyperthread contention accounts for most of the 70-rank spread, and about a
20% residual survives at 36 ranks.  (A prediction made before the run — that
ranks 34 and 35 would be the fast ones, being the only pair whose sibling CPUs
70 and 71 are unused — was **wrong**: the fast ranks are 53 and 18, the slow
ones 4, 5, 8, 17.  The rank-to-CPU pinning is not the linear map that
prediction assumed, and the actual pinning was not determined.)

**Hyperthreading buys nothing here.**  70 ranks finish the pass in 35.79 s;
36 ranks finish the same pass in 34.86 s.  Twice the ranks, slightly *worse*
wall time.

## What a dynamic scheduler was predicted to recover

This section is the prediction made **before** implementing anything; the
measured outcome is two sections below.

Cause matters less than it first appears.  A scheduler cannot create physical
cores, but it does not need to: it removes the *waiting*, by letting whichever
ranks run fast take more photons.  With work handed out on demand, all ranks
finish together and the pass takes the harmonic mean of the individual rank
times instead of the maximum:

```
static      T = max_i t_i
master-slave T = n / sum_i (1/t_i)        (i over the workers)
```

The measured rank times give:

| configuration | static wall | master-slave (rank 0 serving) | gain |
|---|---|---|---|
| **70 ranks** | 35.79 s | **27.82 s** | **1.287x** |
| 36 ranks | 34.86 s | 31.00 s | 1.125x |

Losing rank 0 to the master costs 0.3 s of the 8 s recovered at 70 ranks — one
rank in 70, as expected.

Ranking the four combinations of this one pass:

| | wall time |
|---|---|
| 70 ranks, static (what MoCHII does today) | 35.79 s |
| 36 ranks, static | 34.86 s |
| 36 ranks + master-slave | 31.00 s |
| **70 ranks + master-slave** | **27.82 s** |

Note the interaction: master-slave makes hyperthreading useful again.  Under
static splitting a contended thread holds up the whole pass; under dynamic
splitting it is simply a slower worker that takes a smaller share, and its
contribution is additive.  That is why 70+MS beats 36+MS.

The harmonic mean is an **upper bound**: it assumes zero scheduling overhead and
infinitely divisible work.  Batch granularity is not the limit here (4e7 photons
in batches of `num_send_at_once` = 1e4 is 4000 batches over ~69 workers, ~58
each), but the communication round trip of each batch is real and is paid on
every pass, of which a model runs 3-100.  The measurement below lands below this
bound, as it should.

## Why this differs from Lyman-alpha and from MoCafe

The measurement above should not be read as an argument against master-slave in
general.  **Resonant line transfer is the opposite extreme**: in Lyman-alpha RT
a photon in optically thick gas scatters an enormous number of times before it
escapes, the cost distribution acquires a very long tail, and a handful of
photons can dominate a whole rank's wall time.  In that regime master-slave is
reported to be worth up to a **factor of 10** — not a 15% correction but the
difference between a run finishing and not.  MoCafe's dust RT sits between the
two: photons are scattered rather than absorbed, so its tail is real, which is
why its default is `use_master_slave = .true.`.

MoCHII's ionizing band has no resonance and no long-lived packets, so it sits at
the benign end.  That is a property of the physics being transported, not of the
implementation, and it is the reason the same scheduler choice comes out
differently in the two codes.

## What was implemented

The decision the measurement pointed to — implement master-slave — was carried
out on the same day, in two steps.

**Step 1: one transport pass, one loop.**  `src/main.f90` carried the transport
pass twice, once in the `contains` routine `ionizing_field_and_rates` and once
as a copy for the peel-off imaging pass, so any change to the loop had to be
made in both.  Both are now the single
`ionizing_field_and_rates(label, with_peel)` in the new
`src/ionizing_field_mod.f90`; `main.f90` drops from 471 to 340 lines.  Behavior
is unchanged, and the `peel_cube` regression reproduces bit for bit.

**Step 2: the schedule behind an interface.**  `src/photon_schedule_mod.f90`
holds both schedules behind three calls, so the loop body is written once:

```fortran
call photon_schedule_begin(ntotal [, report_stride])
do while (photon_schedule_next(ip))
   ... transport packet ip ...
end do
call photon_schedule_end()
```

- `par%use_master_slave = .false.` (the MoCHII default; MoCafe defaults it
  `.true.`) walks `ip = p_rank+1, p_rank+1+nproc, ...` — the same indices in the
  same order as the loop it replaced, so results are unchanged bit for bit.
- `.true.` makes rank 0 the master: it serves batches of
  `par%num_send_at_once` (default 10000) consecutive indices to whichever rank
  asks next and **transports nothing itself**.  A rank that runs slow takes a
  smaller share instead of holding up the reduction that closes the pass.

MoCafe keeps the two schedules as two whole routines with the transport loop
inlined in each (`run_master_slave` and `run_equal_number` in
`src/run_simulation_mod.f90`).  MoCHII separates the schedule from the loop body
instead, which is why one loop covers both cases and the peel-off pass gets the
scheduler for free.

**Scope.**  The interface is used by the stellar loop **and** the diffuse LyC
loop, and therefore also by the peel-off imaging pass, which runs the same
routine with `with_peel`.  It is **not** used by the dust emission infrared
transport (`src/dust_emis_transport_mod.f90:269`) or by the SEDust emission
solve (`src/sedust_mod.f90:221`).  Those two distribute **leaves** cyclically,
not photons, and allocate packets to a leaf by stratifying on its share of the
cumulative emission, so the packet index is not a global quantity a master could
hand out; adapting them needs a different interface, not this one.

A run reports the choice in the startup summary
(`Master-slave scheduling   :`, `src/setup.f90:615`).

## What it actually recovered

`tests/g2_hii/hii40_bench.in` with `gas_niter = 1`, 4e7 photons, **70 ranks**,
six runs alternating between the two schedules:

| | min | median | max |
|---|---|---|---|
| transport pass, static | 45.5 s | **52.1 s** | 58.3 s |
| transport pass, dynamic | 43.1 s | **43.4 s** | 43.7 s |
| total wall clock, static | 159.0 s | **174.9 s** | 178.3 s |
| total wall clock, dynamic | 151.5 s | **152.2 s** | 154.7 s |

**The transport pass runs 1.20x faster (52.1 → 43.4 s, 16.7% less time) and the
whole run 1.15x (174.9 → 152.2 s, 13.0% less time).**  Both figures below are
speed ratios, not time reductions; the two differ and it is worth keeping them
apart.  1.20x is under the 1.287x predicted above, as it should be: the harmonic
mean assumes zero communication cost, so it is an upper bound and not a target.

The **variance** is the more striking result.  The static schedule ran between
45.5 and 58.3 s for identical work — a 28% swing driven by whatever else the
node was doing — while the dynamic schedule ran 43.1 to 43.7 s, **1.5%**.  A
schedule that hands out work on demand absorbs the contention instead of
sampling it, so a timing measured under it means something.

## Correctness of the dynamic path

With `launch_sequence = 'sobol'` the two schedules transport the *same* packet
set (the launch is fixed by the global index) and differ only in the order the
contributions are summed, so they can be compared directly.  Static against
dynamic at 16 ranks:

- the diffuse packet count agrees to the last digit in all three passes;
- `rates.h5`: largest relative difference **1.16e-12**;
- `<base>_lines.txt` and `<base>_nebcont.txt` are textually identical, including
  L(Hβ) = 2.262513E+37 on both sides.

Four deadlock boundaries were run and all terminated normally: `np = 1` (no
worker exists, so the request degrades to the static schedule and rank 0
transports everything), `no_photons = 5` (fewer batches than workers),
`num_send_at_once = 1`, and a pass with zero diffuse photons.

## Why the default stays `.false.`

With the default `launch_sequence = 'random'` each rank draws from its own
stream, so which rank processes a photon selects the realization: statistically
equivalent, not numerically identical.  Every recorded gate would move.  The
switch is therefore opt-in, and a run that wants both the speed and a comparable
number should set `launch_sequence = 'sobol'` with it.

## Bit-level reproducibility is not free even without this switch

This matters here because the project verifies regressions by bit identity.

**Confirmed by measurement:** two runs of the *same binary* on the *same input*
at the *same rank count*, with no code change of any kind between them, are
already not bit-identical.  At 16 ranks on `tests/cooling_ne/id_new.in`, only
**5 of the 31 datasets** in `rates.h5` matched bit for bit; `J_nu` differed by
**5.7e-14** in relative terms.  (`<base>_lines.txt` was identical.)

**Reported, not re-verified here:** the cause was diagnosed as Intel MPI
choosing its collective algorithms adaptively from run to run, and pinning them
with `I_MPI_ADJUST_ALLREDUCE=1` (plus the `REDUCE` / `BCAST` / `GATHER`
equivalents) was reported to restore bit identity.  That remedy has **not** been
re-confirmed, and this document does not claim it works.

The practical consequence stands either way: a bit-identity check is a statement
about a *fixed collective algorithm selection*, not about the code alone.  A
comparison that fails at the 1e-14 level has not necessarily found a defect.

## A separate, cheaper finding

Hyperthreading gives this workload nothing under static splitting (36 ranks beat
70).  Anyone running MoCHII on the static schedule can have most of the 70-rank
penalty back by using one rank per physical core, with no code change at all.

## What was not measured

- Only the ionizing-band transport loop was instrumented for the rank-time
  distribution.  The thermal solve, the dust emission solve and I/O were not,
  and in B_bench the pass is only part of the iteration.  (The total wall clock
  in the comparison above does include them.)
- One node, 70 ranks, one grid resolution (L6).  Nothing here says how the
  imbalance scales with rank count, and the interesting regime for a scheduler
  is many more ranks than photons per rank.
- The `hii20` benchmark and the plane-parallel slab geometry were not run.
- The dynamic schedule was not timed under `launch_sequence = 'random'`; the
  timing runs and the agreement runs used different rank counts (70 and 16).
- The `I_MPI_ADJUST_*` remedy above was not re-tested.
