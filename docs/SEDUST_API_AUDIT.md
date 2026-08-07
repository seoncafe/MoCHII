# SEDust v1.00 library API — audit from the MoCHII host

Audited 2026-08-07 against SEDust v1.00, the scalar branch, at
`/home/kiseon/MoCafe/Grain/SEDust_v1.00`. MoCHII's copy under `SEDust/` is
byte-identical to that tree for all 27 library sources and `sed/build_lib.sh`,
so the code this file describes is the code MoCHII links.

**Nothing in SEDust was changed from here.** The audit found five places where
the host had to work around the library. All five are now resolved in the tree
above, four of them by a data-root and status rework in SEDust and the fifth by
rewriting the optics products; this file records what each was, how it was
resolved, and what MoCHII dropped as a result.

## What was measured, and how

A probe program linked against `SEDust/sed/lib/libsedust.a` calls

```fortran
call build_dust(m, model, data_dir, 200, 2.7_wp, 5.0e3_wp, &
                include_euv = .true., status = st)
```

for each of the three models MoCHII can run, then `dust_extinction`. Every file
the process opened was recorded with `strace -e trace=openat`, so the list below
is what the library actually reads, not what its comments say it reads. The
data-root behavior was re-measured with SEDust's own `check_data_root.f90`,
linked against MoCHII's archive and run from a directory that is not a SEDust
`sed/`.

## The normal route works

All three models return `status = 0`, and `dust_extinction` returns `status = 0`:

| model | NLAM | lambda(1) [um] | lambda(NLAM) [um] | C_ext(lambda_1) [cm^2/H] |
|---|---:|---|---|---|
| astrodust | 1762 | 1.0000000E-04 | 3.9810000E+04 | 1.3249854E-21 |
| dl07 | 1823 | 6.2053992E-05 | 3.9810000E+04 | 4.6046705E-25 |
| zubko | 1201 | 1.0000000E-03 | 1.0000000E+04 | 6.9773439E-23 |

Files opened, by model:

| model | opened |
|---|---|
| astrodust | `astrodust/sedust_astrodust.h5`, `dielectric/index_CpaD03`, `dielectric/index_CpeD03`, `dielectric/q_D16graphite.dat`, `release/size_distribution.dat` |
| dl07 | `dl07/sedust_dl07.h5`, `dielectric/index_CpaD03`, `dielectric/index_CpeD03`, `dielectric/index_silD03`, `release/size_distribution.dat` |
| zubko | `zubko/sedust_zubko.h5`, `zubko/Graphitic_Calorimetry_1000.dat`, `zubko/Silicate_Calorimetry_1000.dat`, `zubko/ZDA_BARE_GR_S_Config.dat` |

Every one of those paths is resolved from the `data_dir` argument. No `q_*.dat`
and no `kext_*.dat` is opened by any model.

## Finding 1 (resolved) — `data_dir` resolved the model directories but not the shared ones, and the failure was a runtime abort rather than `status`

`build_dust` resolved `data_dir/<model>/` and
`data_dir/release/size_distribution.dat` from its argument. The dielectric
functions did not come from `data_dir` at all: `q_graphite.f90`,
`q_silicate.f90`, `q_astrodust.f90`, `pah_ld01.f90`, `pah_ioniz.f90` and the two
D16 graphite modules each opened `../data/dielectric/...`, a compile-time path
relative to the working directory, and opened it with no `iostat`. A host that
pointed `data_dir` anywhere else separated the model from the optical constants
its optics were computed on, and the Fortran runtime terminated the image
(`forrtl: severe (29): file not found`) before the builder could set a code.
For an MPI host that is the worst shape a failure can take: one rank dying
uncollected tears the job down before any rank can report.

**Resolved.** `sed_paths.f90` holds the data root; every one of those opens now
goes through `sed_data_path` and carries `iostat`, and `build_dust` sets the
root to its `data_dir` argument for the length of the build and restores it at
every exit. Measured with `check_data_root.x` from a scratch directory: all
three models build with `status = 0` from an absolute `data_dir`, and a
directory that does not exist comes back as `status` (1 / 1 / 9) rather than as
an abort — `PASS`.

MoCHII dropped `chdir(par%sed_workdir)` and the parameter with it.

## Finding 2 (resolved) — `lam_min` extended for two models and truncated for the other two

For astrodust and DL07, `lam_min` prepended log-spaced points below the stored
axis. For zubko and `from_files` it did the opposite: an explicit `lam_min` won
over `include_euv`, and since those models' grids are their own optics tables,
the request became a cut. A host with one `build_dust` call and one physically
motivated floor passed that floor for every model and silently truncated zubko,
from 1.0e-3 um to whatever it asked for, with nothing in `status` to say so.

**Resolved.** `lam_min` is now a coverage requirement rather than a grid: a
model that cannot reach the asked-for wavelength refuses to build and says so
(`status = 4` in `build_dust`'s vocabulary). Nothing is cut.

MoCHII still passes no `lam_min`, for a different reason than before — refusing
is not what a host wants from a floor its own band never reaches, and with
`include_euv = .true.` the astrodust axis already reaches 1.0e-4 um
(12398.4 eV). Whether the grid spans the transported band is checked with the
numbers in `gas_opacity_mod` at setup, for every model alike.

## Finding 3 (resolved) — `status` codes meant different things model by model

The four builders each numbered their own stages, so 5 was an unreadable
extinction table for astrodust and DL07, 6 for zubko, 9 for `from_files`, and 6
meant EUV spheroid optics for astrodust. A host using the single entry point
still had to branch on the model name to read a code.

**Resolved.** `build_dust` maps every builder's code onto one vocabulary:

| code | meaning |
|---:|---|
| 0 | built |
| 1 | optics table (Q table, or a component's optics) |
| 2 | size distribution |
| 3 | dielectric function |
| 4 | `lam_min` not coverable by this model |
| 5 | extinction table |
| 6 | calorimetry |
| 7 | grid inconsistent between components |
| 8 | EUV spheroid optics unavailable (no T-matrix registered) |
| 9 | model definition (the ZDA config, or the descriptor) |
| 90 | model name not one of the four |
| 91 | `from_files` without a descriptor |
| 92 | `zubko_optics` not one of `zda` / `mie_d03` |

and an optional `message` returns the same thing in words, for a host that has
to print one line before a collective abort.

MoCHII dropped its `st_kext = 5 / 6` mapping — the extinction table is 5 for
every model — and prints `message` on the codes it has nothing more specific to
say about.

## Finding 4 (resolved) — the DL07 text route ignored `data_dir`, and its message named a directory it did not open

`stored_q_on_model_grid` tried the HDF5 product first and then the text tables
at a hard-coded `'../data/dl07/'`, so a DL07 host that set `data_dir` elsewhere
lost the text route without being told and the build solved the optics from the
dielectric functions instead — a different number. The message it printed named
`../data/dl07/` even on the HDF5 route, where no file under that directory other
than `sedust_dl07.h5` was opened.

**Resolved.** The hard-coded directory is gone, and the message now names the
file that was read:

```
 sed_init_dl07: optics read from <data_dir>/dl07/sedust_dl07.h5
```

## Finding 5 (resolved) — `/kext` and `/qtable` of the same product were not computed from the same numbers

`calc_kext.x`, which writes `/kext`, called `build_dl07` and `build_zubko`
directly. `stored_q_h5` was set only inside `build_dust`, so those two builds
fell to the text `q_*.dat` optics — written with 7 significant digits — while
`/qtable` of the same file carried full double precision. The curve a host was
served was therefore the size integral of a 7-digit copy of the optics its model
was built on. Astrodust escaped it: its `/qtable` is itself the 7-digit T-matrix
table, so both halves of that product held one set of numbers.

Measured then, by rebuilding each model with `build_dust` and comparing
`dust_extinction` (the stored `/kext`) with `size_integrated_extinction` (the
same integral over the model as built from `/qtable`) — worst single wavelength,
relative on the cross sections, absolute on (albedo, albedo\*&lt;cos&gt;):

| model | lambda >= 0.0912 um | ionizing band | C_ext(0.55 um) |
|---|---|---|---|
| astrodust | 0.0 / 0.0 | 3.3e-16 / 1.1e-16 | 0.0 |
| dl07 | 1.05e-7 / 3.1e-8 | 9.6e-8 / 3.0e-8 | 1.4e-9 |
| zubko | 8.5e-8 / 3.2e-8 | 5.3e-8 / 1.9e-8 | 1.4e-8 |

**Resolved**, by rewriting the products so that both halves come from one set of
numbers. The same measurement now:

| model | lambda >= 0.0912 um | ionizing band | band s_abs | C_ext(0.55 um) |
|---|---|---|---|---|
| astrodust | 0.0 / 0.0 | 3.3e-16 / 1.1e-16 | 0.0 | 0.0 |
| dl07 | 2.7e-15 / 6.7e-16 | 4.2e-16 / 3.3e-16 | 0.0 | 0.0 |
| zubko | 5.8e-16 / 2.2e-16 | 6.3e-16 / 2.2e-16 | 6.9e-15 | 7.2e-15 |

That is rounding on a double. `tests/grain_kext/check_kext_table.f90` therefore
holds all three models to one tolerance, `TOL_KEXT_NODE = 1e-12`, in place of
the 1e-9 / 1e-6 split the 7-digit route had needed.

## What MoCHII carries now

| finding | what MoCHII used to do | what it does |
|---|---|---|
| 1 `data_dir` half-honored, abort not `status` | `chdir(par%sed_workdir)` around the build, restoring the working directory | nothing — `par%sed_workdir` is gone, and `par%sed_data_dir` is the whole of it |
| 2 `lam_min` truncated zubko | passed no `lam_min`, to avoid the cut | passes no `lam_min`, because its own band never reaches a floor any model would refuse; the coverage check stays in `gas_opacity_mod` |
| 3 status vocabulary not shared | mapped the extinction-table code by model (5 / 5 / 6) | reads 5, and prints SEDust's `message` on anything else |
| 4 DL07 text path and its message | read the HDF5 products only | unchanged — the text route is no longer reachable by accident either |
| 5 `/kext` and `/qtable` from different numbers | `tests/grain_kext` held astrodust to 1e-9 and DL07 / Zubko to 1e-6 | one tolerance, 1e-12, for all three |
