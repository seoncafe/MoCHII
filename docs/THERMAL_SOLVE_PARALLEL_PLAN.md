# Plan: distribute the thermal solve across the node-local ranks

- Date: 2026-07-28
- Target: `gas_thermal_update` in `src/thermal_mod.f90`
- Status: implemented and validated locally

## 1. Why

`gas_thermal_update` runs on one rank per node. Every other rank on that node
idles through it. Measured on this machine, that idle time is now the larger
part of a run.

| case | leaves | photons | iterations | transport | everything else |
|---|---|---|---|---|---|
| `tests/pdr/pdr_dense.in` | 110592 | 4e6 | 30 | 0.76 min (4.6%) | **15.56 min (95.4%)** |
| `tests/g2_hii/hii40.in` | AMR L6 | 4e7 | 84 | 71.5 min (56.3%) | **55.6 min (43.7%)** |
| `tests/g2_hii/hii20.in` | AMR L6 | 4e7 | 83 | 75.6 min (59.4%) | **51.7 min (40.6%)** |

All three ran with 8 MPI ranks on one node. The transport column is already
distributed over all 8; the remainder is almost entirely the thermal solve and
is executed by a single rank. The fraction rises as the leaf count grows
relative to the photon count, and rises further with `use_metals`, because each
trial temperature then re-solves the metal cascade.

Amdahl with 8 ranks per node, taking the measured fractions:

```text
pdr_dense : 1/(0.046 + 0.954/8) = 6.1x  on total run time
hii40     : 1/(0.563 + 0.437/8) = 1.6x
```

## 2. Why it is serial today

This is a deliberate choice, not a constraint. `src/thermal_mod.f90:140-144`:

```text
only the node-local rank 0 solves; jt_ion is ALLREDUCEd over
MPI_COMM_WORLD so every node's rank 0 is identical.  The other node
ranks skip the (expensive) bisection-with-metals and receive the
state through shared memory + the metrics by broadcast.
```

The gas state arrays live in node-shared memory, so having one writer avoids
both duplicated work and any question about who writes what. The cost is that
the other ranks on the node do nothing at all during the solve.

Note the redundancy that already exists across nodes: **every** node's rank 0
solves the whole grid. This plan does not change that. It fills in the idle
ranks within each node, which is where the ratio is worst.

## 3. The loop is independent

Verified by reading `gas_thermal_update` (`thermal_mod.f90:152-215`) and
everything it calls:

1. **Each leaf touches only its own index.** `net_rate(il, ...)` reads
   `gas_nH(il)`, `gamma_HI/HeI/HeII(il)`, `sec_*(il)`, `amr_grid%rhokap(il)`,
   `g0_fuv(il)`, and through `species_fractions` the block `egam(ie)%g(it,il)`.
   `solve_ion_cell`, `charge_neutrality_ne`, `metal_cooling`,
   `metal_cooling_H`, `metal_heating` and `metal_freefree` are all evaluated
   at one leaf.
2. **No loop-carried state.** `xHI`, `xHeI`, `xHeII` are re-seeded from the
   leaf at line 162 on every pass, so the bisection of one leaf cannot inherit
   anything from the previous one.
3. **Reads and writes are already separated.** The solve fills the local
   arrays `xHI_new`, ..., `te_new`; the shared arrays are written in a second
   loop at lines 220-226. The design is already in the shape parallel
   execution needs.
4. **The only cross-leaf quantities are reductions**: `max_dx`, `max_dte`
   (maxima) and `sum_dxv`, `sum_xv`, `sum_dtev`, `sum_tev` (sums).

Point 3 also means the ranks may write the shared arrays in place at their own
leaves without a hazard: within one leaf, every read of the old value
(lines 162, 187-189, 208-214) precedes the write.

## 4. Design

Distribute over `mpar%hostcomm`, the communicator that already backs the
shared-memory windows.

1. Remove the `if (mpar%h_rank == 0)` guard around the solve.
2. Stride the solve loop and the write-back loop identically:

   ```fortran
   do il = mpar%h_rank + 1, gas_nleaf, h_size
   ```

3. Replace the four `MPI_BCAST` calls with reductions over `hostcomm`:
   `MPI_ALLREDUCE(..., MPI_MAX, ...)` for `max_dx` and `max_dte`,
   `MPI_ALLREDUCE(..., MPI_SUM, ...)` for the four accumulators. Form
   `dx_vol` and `dte_vol` from the reduced sums, not before.
4. Publish the five gas-state shared windows with their existing fence
   epochs, then synchronize the node-local writers.
5. Broadcast rank 0's four reduced metrics over `MPI_COMM_WORLD`. This is
   required because different nodes may have different host sizes or
   reduction orders; without it, ranks could make different convergence
   decisions and later deadlock in a world collective.

`mpar` stores `hostcomm` and `h_rank` but not the host size, so either add an
`h_size` field alongside them in `src/define.f90` and set it where `hostcomm`
is built, or call `MPI_COMM_SIZE(mpar%hostcomm, ...)` once and cache it in the
module. Adding the field is preferable — the size is a property of the MPI
context and other routines will want it.

The local arrays `xHI_new`, ... stay full length. At 110592 leaves that is
5 x 110592 x 8 B = 4.4 MB on each rank, which is not worth the complication of
packed indexing.

## 5. What is guaranteed and what changes

**The gas state is bit-identical.** Every leaf's arithmetic is untouched and
each leaf is computed by exactly one rank, so `gas_xHI`, `gas_xHeI`,
`gas_xHeII`, `gas_ne` and `gas_Te` reproduce the serial result exactly. This
is a strong property and the validation in section 7 should assert it rather
than a tolerance.

**`max_dx` and `max_dte` are exact.** `MPI_MAX` over a partitioned set gives
the same value as a serial maximum regardless of order.

**`dx_vol` and `dte_vol` change in the last bits.** The four sums are
accumulated in a different order, so their rounding differs. These two numbers
only gate the convergence test, and a run that sits near the tolerance could
in principle stop one iteration earlier or later. That is acceptable, but it
must be stated rather than discovered: the change is not bit-reproducible in
the convergence *metrics*, only in the *state*.

If exact reproduction of the metrics is ever wanted, accumulate them with a
compensated sum over a fixed leaf ordering after the reduction of the state.
That is not proposed here.

## 6. Load balance and the choice of MPI over OpenMP

**Stride, not contiguous blocks.** The bisection cost varies strongly from one
leaf to the next: a leaf whose net heating does not bracket a root costs two
`net_rate` calls, while a bracketed leaf costs up to 62. Those two populations
are spatially clustered (the ionized interior clamps, the front brackets), so
a contiguous partition would be badly unbalanced while a stride mixes them.

**MPI rather than OpenMP.** `species_mod` caches the transition coefficients
of the current leaf and temperature in the module variables `cch_gam`,
`cch_ci`, `cch_rec`, `cch_cxi`, `cch_cxr`, filled by `species_ne_prepare` and
consumed by `species_ne_cached`. Separate MPI processes each hold their own
copy, so this is safe as written. OpenMP threads would share those arrays and
race; making it thread-safe means `threadprivate` or threading the cache
through the call chain. The Makefile links `-qopenmp` today but the source
contains no directives, so nothing is lost by staying with MPI.

## 7. Validation

1. **Bit-identical state.** Run `tests/pdr/pdr_dense.in` before and after with
   the same rank count and require `gas_xHI`, `gas_xHeI`, `gas_xHeII`,
   `gas_ne`, `gas_Te` to match exactly, not within a tolerance.
2. **Rank-count independence.** Run at 1, 2, 4 and 8 ranks on one node and
   require the same state. This is the test that catches an ownership mistake
   in the stride, which a single rank count cannot.
3. **Convergence metrics.** Record `dx_vol` and `dte_vol` at every iteration
   and confirm they agree with the serial values to roundoff and that the
   iteration at which the run converges is unchanged for a case that is not
   sitting on the tolerance.
4. **Timing.** Report the transport and non-transport split from section 1 for
   the same cases and confirm the predicted scaling.
5. **Multi-node.** Run on two nodes and confirm the state still matches, since
   the cross-node redundancy is unchanged but the intra-node partition is new.

Test 2 is the important one. Tests 1 and 2 should become a committed gate, not
a one-time check.

### Local implementation result (2026-07-28)

The optimized Intel build was run on one node using the 262144-leaf HII20
smoke model, 100000 Sobol photons, one thermal iteration, and metals enabled.
The Sobol launch was used so the three-rank transport field was identical
between the serial-thermal baseline and the parallel-thermal build.

| build | MPI ranks on node | thermal solve |
|---|---:|---:|
| Git HEAD baseline (serial thermal loop) | 3 | about 25.8 s |
| node-local MPI thermal loop | 3 | 8.862 s |

This is a 2.91x thermal-solve speedup with three ranks. The five persisted gas
state datasets (`x_HI`, `x_HeI`, `x_HeII`, `n_e`, and `T_e`) were compared
with `h5diff` and matched exactly. The same optimized build also completed at
one rank, and its printed convergence metrics agreed with the three-rank run.
The reusable exact-state comparison is
`tests/thermal_parallel/compare_state.sh`.

## 8. Risks

- **Ownership mistakes in the stride** are the main failure mode: a leaf
  solved by nobody keeps its previous value and a leaf solved by two ranks is
  written twice with the same value, so neither shows up as a crash. Test 2
  above is what detects this.
- **The write-back loop must use the same stride as the solve loop.** If they
  diverge, a rank writes leaves it did not solve, reading uninitialized
  entries of its local arrays.
- **`h_size` must come from `hostcomm`**, not from `mpar%nproc`. On a
  multi-node run they differ and using the wrong one silently drops leaves.

## 9. Out of scope

- Cross-node distribution of the solve. Every node's rank 0 currently solves
  the whole grid; removing that redundancy means synchronizing the shared
  arrays between nodes, which the present design deliberately avoids.
- The thermal solver's algorithm. Finding 3.8 of
  [REVIEW_VERIFICATION_2026-07-27.md](REVIEW_VERIFICATION_2026-07-27.md)
  (single monotonic root, silent boundary clamping) is a separate change. It
  should land either well before or well after this one — changing the root
  finder and its parallel structure together would make a regression hard to
  attribute.
- OpenMP inside each rank, which would need the `species_mod` cache made
  thread-safe first.

## 10. Sequencing

This performance change was implemented after the preceding I/O and build
checkpoint had been validated, committed, and pushed. It does not alter the
thermal root-finding algorithm or its heating/cooling physics.
