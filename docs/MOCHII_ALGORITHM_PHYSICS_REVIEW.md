# MoCHII Source Review: Algorithms, Physics, and Software Engineering

Review date: 2026-07-27

## 1. Scope and method

This review is based on direct inspection of the implementation in
`src/*.f90`, the build configuration, and the executable regression scripts.
It does not rely on the descriptive documents as evidence for how the code
works.

The physics and algorithm review covered:

- photon generation, frequency sampling, AMR ray walking, and path-length
  estimators;
- H/He photoionization and recombination equilibrium;
- electron-density closure, metal ion stages, line cooling, and PDR terms;
- thermal root finding and the nonlinear iteration/output sequence;
- AMR re-refinement and state remapping;
- diffuse recombination radiation and live dust-opacity feedback.

The software-engineering review covered:

- build reproducibility and compiler portability;
- I/O error propagation and output integrity;
- input and atomic-data validation;
- MPI error handling, shared-memory ownership, and scaling limits;
- regression-test enforceability;
- global state, module coupling, and maintainability.

Judgments are classified as:

- **Confirmed defect**: follows directly from the executed code path or
  equation.
- **Conditional confirmed defect**: definite when the stated option or regime
  is used.
- **Numerical risk**: the implementation permits a failure, but the inspected
  sources alone do not prove that existing production inputs trigger it.
- **Model limitation**: a boundary of the implemented physics rather than an
  implementation error.

Atomic coefficients were not validated against external literature in this
review. No production source file was modified.

## 2. Executive summary

### 2.1 Physics and algorithm findings

| Priority | Classification | Finding | Main consequence |
|---|---|---|---|
| P0 | Confirmed defect | The final gas update is not followed by an unconditional radiation/rate consistency pass | Output gas state can be paired with rates and dust heating from the preceding state |
| P0 | Confirmed defect | AMR re-refinement remaps state by sampling one old cell at each new leaf center | H mass, ion counts, electron count, thermal content, and dust optical depth are not conserved under coarsening |
| P0 | Conditional confirmed defect | `laursen09_live` loses the initial `taumax`/`tauhomo` normalization during opacity refill | Dust optical depth can jump after the first gas update |
| P0 | Confirmed physics defect | The H-collider C II/O I two-level saturation denominator omits the upward collision rate | Dense-gas cooling is overestimated |
| P1 | Confirmed coupling defect | H-collider metal cooling is controlled by the `grain_pe` heating switch | Disabling photoelectric heating also disables an independent cooling process |
| P1 | Confirmed sampling defect | Diffuse equal-energy packet count and intra-channel spectrum do not preserve the intended energy distribution | Weak diffuse luminosity can disappear and the spectrum is biased low |
| P1 | Numerical risk | Metal-stage fractions use an unscaled product chain | Extreme rate ratios can produce overflow and NaNs |
| P1 | Algorithmic limitation | The thermal solver assumes one monotonic root and silently accepts temperature bounds | Multiple or hidden thermal equilibria can be missed |

### 2.2 Software-engineering findings

| Priority | Classification | Finding | Main consequence |
|---|---|---|---|
| P0 | Confirmed defect | HDF5 errors can be overwritten by later cleanup calls, and output callers commonly ignore final status | Partial or failed output may be reported as successfully written |
| P0 | Confirmed build defect | The GNU `HDF5=0` build does not compile because preprocessor-generated stubs exceed the free-form line limit | The advertised optional-HDF5 configuration is broken |
| P0 | Confirmed build defect | `-qopenmp` is linked for every compiler, including GNU `mpif90` | The GNU link command fails |
| P1 | Build-system risk | No Fortran module dependencies are declared; dependencies and external-library paths are hard-coded | Parallel builds are unsafe and builds are machine-specific |
| P1 | Validation defect | Namelist, atomic-data, and several numerical inputs lack complete checked parsing | Typos or malformed tables can cause silent mode changes, truncation, or out-of-bounds access |
| P1 | Test-infrastructure defect | Several scripts print `PASS`/`FAIL` but do not return a failing process status | A CI runner could accept a failed scientific regression |
| P1 | MPI robustness risk | Most MPI return codes are ignored and fatal paths mix `STOP`, `ERROR STOP`, `MPI_FINALIZE`, and `MPI_ABORT` | Rank-local failures can be hard to diagnose and may leave other ranks blocked |
| P2 | Maintainability risk | A large public global module and broad imports couple most of the program | Local reasoning, unit testing, and safe refactoring are unnecessarily difficult |

The first four physics defects and the first three software defects should be
addressed before adding major new features.

## 3. Detailed physics and algorithm findings

### 3.1 [P0, confirmed] Final gas state and radiation rates can be one iteration apart

Evidence:

- `src/main.f90:116-165` transports packets and computes `jt_ion`, photo-rates,
  and heating using state `x_k`.
- `src/main.f90:168-194` solves for state `x_(k+1)` and refills opacity.
- `src/main.f90:203-208` exits immediately on convergence or at the iteration
  limit, without transporting through the new opacity.
- `src/main.f90:283-292` then writes the current gas state and uses the
  previously calculated rates and `heat_dust`.

For a normal run with `ion_peel = .false.`, the final products therefore mix:

```text
gas fractions, ne, Te, line emissivities : x_(k+1)
J, Gamma, gas heating, dust heating      : radiation field at x_k
```

At convergence, the mismatch may be limited by the convergence tolerance. If
the iteration cap is reached and `require_convergence = .false.`, its magnitude
is not bounded.

When `ion_peel = .true.`, the extra imaging transport at
`src/main.f90:232-279` rebuilds the tally and rates using the final state. That
option-dependent correction is not a suitable guarantee of scientific state
consistency.

Recommended change:

1. Always run a final transport-and-rates pass, without another gas solve, after
   the nonlinear loop.
2. Rebuild diffuse emission from the final `ne`, ion fractions, and `Te`.
3. Compute output rates and dust temperature only after this pass.
4. Record that the final consistency pass was completed.

### 3.2 [P0, confirmed] AMR re-refinement is not conservative

`src/amr_refine_mod.f90:99-113` builds a new leaf list and assigns each new
leaf the state of the old leaf containing the new leaf center.

This is valid as piecewise-constant prolongation when one old leaf is only
subdivided. It is not valid when heterogeneous old leaves are coarsened into a
larger new leaf: the state of one center-containing cell replaces the
volume-weighted state of the entire region.

The following global quantities are not conserved:

```text
total H number             = sum(nH V)
number in each ion stage   = sum(nH x_i V)
electron number            = sum(ne V)
dust extinction integral   = sum(rhokap V)
thermal content            = sum(u V)
```

This is most consequential at the density and ionization discontinuities that
motivate re-refinement. In addition, `species_resize` at
`src/species_mod.f90:164-175` zeros the metal photo-rate arrays, so the first
new-grid metal opacity is not based on the pre-remap radiation state.

Recommended change:

- Compute old/new leaf volume overlaps and remap extensive quantities.
- Remap H and He stage numbers, electron number, dust extinction integral, and
  an explicitly defined thermal energy.
- Recover intensive quantities only after dividing by the new volume.
- Remap metal rates or require a new-grid radiation pass before metal opacity
  is considered consistent.
- Print and assert global conservation errors before and after rebuilding.
- As a short-term safety option, prohibit coarsening and allow refinement only.

### 3.3 [P0, conditional confirmed] `laursen09_live` loses dust normalization

The initial grid path computes the live dust opacity at
`src/grid_mod_amr.f90:234-280`. It then derives `opac_norm` from `taumax` or
`tauhomo` and multiplies every leaf opacity at
`src/grid_mod_amr.f90:284-308`.

After a gas update, `src/gas_opacity_mod.f90:203-213` recomputes
`laursen09_live` opacity from metallicity, density, ion fraction, cross
section, and DGR, but does not apply the original `opac_norm`. The normalization
factor is not retained as persistent model state.

Therefore, when `laursen09_live` is combined with `taumax > 0` or
`tauhomo > 0`, the requested dust scale disappears at the first opacity refill.
This causes a global opacity jump unrelated to physical dust survival.

Recommended change:

- Store and broadcast a persistent `dust_opac_norm`.
- Apply it in every live-opacity refill.
- Define the semantics as normalization of a neutral reference opacity, with
  the survival factor changing only the relative opacity afterward.
- Report the realized pole or mean optical depth each iteration.

### 3.4 [P0, confirmed] Incorrect dense-limit C II/O I two-level cooling

`src/species_mod.f90:485-513` calculates H-impact C II 158 micron and O I
63 micron cooling as:

```text
Lambda_code = n_ion n_H q_lu DeltaE / (1 + n_H q_ul/A)
```

For a two-level steady state,

```text
n_u(A + n_H q_ul) = n_l n_H q_lu
n_ion = n_l + n_u
```

the correct expression is:

```text
Lambda = n_ion n_H q_lu DeltaE
         / (1 + n_H(q_ul + q_lu)/A)
```

The upward rate `q_lu` is missing from the code denominator. The low-density
limit is correct, but the LTE/high-density limit is not.

The high-density overestimate is `1 + q_lu/q_ul`:

- C II: about 1.80 at 100 K and approaching 3 at high temperature.
- O I: approaching 1.6 at high temperature.

Recommended change:

- Use `1 + nHI*(qul + qlu)/A`.
- Preferably add H colliders to the generic n-level solver so that the manual
  two-level expression and the electron-collider solver do not diverge.
- Test both the low-density scaling and the high-density Boltzmann limit.

### 3.5 [P1, confirmed] H-impact cooling is incorrectly gated by `grain_pe`

`src/thermal_mod.f90:78-92` calls `metal_cooling_H` only inside
`if (par%grain_pe)`.

H collisions with C II and O I are independent of whether grain
photoelectric heating is enabled. With the present coupling, changing a
heating switch also removes an unrelated cooling channel.

Recommended change:

- Move `metal_cooling_H` to the general `use_metals` cooling path.
- If backward-compatible control is needed, add a separate
  `metal_H_cooling` option.

### 3.6 [P1, confirmed] Diffuse packet luminosity and spectrum are biased

#### Packet-count rounding

`src/diffuse_mod.f90:100-109` sets:

```text
N_diff = nint(L_diff/L_packet)
```

Every emitted diffuse packet then carries the stellar `L_packet` at
`src/diffuse_mod.f90:187-189` and `:273-275`. The emitted diffuse luminosity is
therefore `N_diff L_packet`, not generally `L_diff`.

- Absolute error can reach `0.5 L_packet`.
- If `L_diff < 0.5 L_packet`, the diffuse field vanishes.
- Relative error is unbounded for a weak diffuse component.

#### Intra-channel energy sampling

Channels are selected by energy luminosity, which is appropriate for
equal-energy packets. Inside a ground recombination continuum, however,
`src/diffuse_mod.f90:149-167` and `:241-258` sample the adopted photon-number
model:

```text
p_N(E) proportional to exp[-(E-E_th)/kT]
```

Within that same exponential approximation, an equal-energy packet must be
sampled from:

```text
p_E(E) proportional to E exp[-(E-E_th)/kT]
```

The present implementation therefore shifts packet energy toward lower
frequencies. The flat He I two-photon photon-number spectrum has the same
energy-weighting issue.

Samples above `eion_max` are also clipped to `0.999 eion_max`, creating an
artificial pile-up in the last bin while the luminosity calculation still uses
the unclipped mean energy.

Recommended change:

1. Give each diffuse packet `L_diff/N_diff`, or use unbiased stochastic
   rounding with an explicit weak-field policy.
2. Sample the energy-luminosity CDF inside every channel.
3. Normalize a consistently truncated spectrum and account separately for
   out-of-band luminosity.

### 3.7 [P1, numerical risk] Metal-stage ratio products can overflow

`src/species_mod.f90:563-593` and the cached path at `:416-442` calculate
successive products of:

```text
r_i = R_ion(i)/R_rec(i)
```

in linear floating-point space. In low-density, strongly irradiated cells,
these products can overflow. Normalization can then produce `Inf`, zero, and
eventually `0*Inf = NaN`.

Recommended change:

- Construct stage weights in log space and normalize with log-sum-exp.
- Handle exact zero and infinite rate limits explicitly.
- Assert that fractions are finite, non-negative, and sum to one.

### 3.8 [P1, algorithmic limitation] Thermal root finding assumes one monotonic root

`src/thermal_mod.f90:161-184` evaluates net heating at `te_min` and `te_max`.
It clamps to a boundary when the endpoint signs do not provide the expected
bracket; otherwise it performs one log-temperature bisection.

Metal line cooling, collisional ionization, charge exchange, and grain
charging can produce a non-monotonic net heating curve. Equal endpoint signs
do not prove that there is no interior root, and multiple roots may contain
both stable and unstable branches.

Recommended change:

- Scan a log-temperature grid for all sign changes.
- Refine each bracket and select a thermally stable root near the previous
  iteration's temperature.
- Record `root_found`, `clipped_low`, `clipped_high`, and `multiple_roots`.
- Do not declare a run converged if a significant volume is boundary-clipped.

## 4. Physics and algorithm components that appear sound

The following central components are internally consistent under the reviewed
assumptions.

1. **Path-length estimator conversion**

   `src/gas_rates_mod.f90:116-131` converts the code-length tally using
   `1/(V_code distance2cm^2)` to obtain `4 pi J_nu dnu`. Multiplication by
   `sigma/(h nu)` and by the excess-energy fraction produces the correct units
   for photoionization and heating.

2. **Analytic direct attenuation tally**

   The `(exp(-tau_in)-exp(-tau_out))/kappa` integration in
   `src/raytrace_amr.f90:225-265` is the exact cell integral of the attenuated
   path length and has the correct zero-opacity limit. The forced-first
   interaction, albedo weighting, and Russian roulette form a structurally
   unbiased scattered estimator.

3. **H/He equilibrium and charge closure**

   The H and He stage ratios, electron-density closure, threshold treatment,
   and excess-energy heating are mutually consistent.

4. **Stellar equal-energy source CDF**

   Stellar and external source components and frequency bins are selected by
   energy luminosity, matching equal-energy packets. Threshold-aligned band
   edges also reduce group errors at ionization edges.

5. **Slab/external angular normalization**

   Lambert sampling with `mu=sqrt(U)`, face-area weighting, `pi I A`
   normalization, and the escaping-intensity conversion are consistent.

## 5. Physics model limitations

These are interpretation limits, not necessarily implementation defects.

- The neutral/PDR cooling set is concentrated on H-impact C II 158 micron and
  O I 63 micron. Cold dense gas requiring C I, additional O I levels,
  molecular cooling, cosmic-ray heating, or gas-grain thermal exchange is not
  a complete thermal model.
- Electron n-level tables may be clamped at a table boundary below their
  supported temperature range. This should be reported separately from a
  genuinely evaluated low-temperature rate.
- Finite frequency groups still require `nnu` convergence tests, especially
  for hard spectra and metal edges.
- Random packet histories change the nonlinear iteration noise unless common
  random numbers or a packet-keyed counter generator are used.
- The Sobol implementation at `src/qmc_mod.f90:186-205` repeats indices modulo
  `2^32`. Enforce `nphotons <= 2^32` or extend the direction numbers.

## 6. Detailed software-engineering findings

### 6.1 [P0, confirmed] I/O failures can be masked and reported as success

The HDF5 dataset writers reuse one `ierr` for creation, writing, dataset close,
and dataspace close:

- `src/hdf5io_mod.f90:335-346`
- `src/hdf5io_mod.f90:365-376`
- `src/hdf5io_mod.f90:389-400`

If dataset creation or writing fails, a later successful close can overwrite
the nonzero error code. The final test then sees zero. `hdf5_open_new` has a
similar pattern at `src/hdf5io_mod.f90:422-435`, where the property-list close
can overwrite the file-creation error.

Read paths show the same risk, for example
`src/hdf5io_mod.f90:1632-1636`: the result of `h5dread_f` is overwritten by
`h5dclose_f`.

At the call sites, output routines commonly initialize `status`, issue many
I/O calls, close the file, and print a success message without checking the
final status. Examples include:

- `src/gas_rates_mod.f90:258-366`
- `src/jtally_mod.f90:286-327`
- `src/ion_peel_mod.f90:402-514`
- `src/lines_mod.f90:68-73` and `:343-350`

`emis_write_state` also creates a new local `status = 0` rather than propagating
the caller's status.

Consequences include incomplete datasets, corrupt output, leaked handles, and
false success messages.

Recommended change:

- Preserve the first error in a dedicated variable.
- Do not call write/close operations on handles that were never created.
- Use separate cleanup error variables and never replace the primary error.
- Return a checked status from every output routine and stop collectively if a
  required scientific product cannot be written.
- Write to a temporary file and atomically rename only after successful close
  and validation.
- Add fault-injection tests for invalid paths, full filesystems, and failed
  dataset creation.

### 6.2 [P0, confirmed] The optional-HDF5 GNU build is broken

The non-HDF5 branch of `src/hdf5io_mod.f90:1639-1818` generates many complete
subroutines through multiline C-preprocessor macros. After preprocessing,
these become lines longer than the standard 132-character free-form limit.

A direct GNU diagnostic compile with the Makefile's `HDF5=0` flags fails at
the stub instantiations with line-truncation and malformed declaration errors.
This occurs in both release-style and debug-style GNU commands.

Recommended change:

- Replace the macro-generated subroutines with normal Fortran procedures or
  generate a `.f90` source file at build time with proper line wrapping.
- As a temporary workaround only, GNU could add `-ffree-line-length-none`, but
  removing the preprocessor dependence is safer.
- Compile both `HDF5=0` and `HDF5=1` in CI.

### 6.3 [P0, confirmed] GNU linking receives an Intel-only OpenMP flag

`Makefile:62` appends `-qopenmp` unconditionally. The detected `mpif90` in the
review environment invokes GNU Fortran, which rejects `-qopenmp`; GNU requires
`-fopenmp`.

No OpenMP directives were found in the main `src` tree, so the link flag may
not be needed for MoCHII itself. If SEDust requires OpenMP, its compile and
link interface should expose that requirement explicitly.

Recommended change:

- Select OpenMP through compiler-aware build logic.
- Use `-qopenmp` for Intel and `-fopenmp` for GNU only when required.
- Prefer a build-system OpenMP package check over literal flags.

### 6.4 [P1] Build configuration is not portable or dependency-safe

Issues in the current `Makefile` include:

- HDF5 defaults to a machine-specific `/data/opt/hdf5_intel`.
- Static HDF5 archive names and transitive libraries are hard-coded.
- CFITSIO and SEDust discovery are also path assumptions rather than checked
  dependencies.
- Fortran module dependencies are deliberately omitted at
  `Makefile:117-120`; `make -j` can compile a consumer before its `.mod` file.
- `default: clean` forces a full in-tree rebuild, but direct `make MoCHII.x`
  can still reuse stale module files.
- Intel, GNU, and linker user flags use inconsistent variable names:
  `Fextra`, `EXTRAFLAG`, and `extra`.
- The Intel release profile enables `-no-prec-div -fp-model fast=2`, while the
  GNU profile follows different floating-point semantics.
- Object and module files are written into `src`, mixing generated artifacts
  with sources.

Recommended change:

- Adopt CMake, Meson, or fpm, or generate accurate Make dependencies.
- Discover HDF5 through its compiler wrapper or package metadata and make
  FITS/HDF5 capabilities explicit.
- Use an out-of-tree build directory.
- Provide separate `debug`, `release-precise`, and opt-in `release-fast`
  profiles.
- Add a compiler/build matrix for GNU and Intel with HDF5 both on and off.

### 6.5 [P1] Floating-point reproducibility policy is implicit

The Intel default uses aggressive floating-point transformations. MPI
reductions also change summation order with task topology. Consequently,
compiler choice, optimization profile, and MPI process count may affect
iteration convergence and weak physical components.

Recommended change:

- Make precise floating-point semantics the scientific default.
- Keep fast-math as an explicit, validated performance option.
- Record compiler ID/version, compiler flags, MPI size, source revision,
  random seeds, QMC seed, and relevant data-table identifiers in every output.
- Define tolerance-based reproducibility requirements across compiler and MPI
  counts; use deterministic or compensated reduction where a stronger
  guarantee is needed.

### 6.6 [P1] Configuration and table parsing are incompletely validated

#### Namelist input

`src/setup.f90:19-36` uses fixed 128/256-character buffers and opens/reads the
namelist without `iostat` or `iomsg`. A missing file, malformed namelist, or
truncated path therefore produces compiler-dependent runtime diagnostics.
Every MPI rank also reads the file independently.

Several values lack explicit validation:

- positive `te_fixed`, `te_min`, `te_max`, `gas_tol`, and `gas_tol_te`;
- `xHI_init`, `xHeI_init`, `xHeII_init` ranges and He-stage sum;
- non-negative `He_abund`;
- `conv_crit` enumeration (`main.f90` treats anything except `vol` as `cell`);
- `refine_lbase <= refine_lmax`, valid refinement iteration, and enough
  post-refinement iterations;
- positive `nnu_fuv` when FUV is enabled.

`no_photons` is stored as a real and then converted to an `int64` packet
count. Non-integral values are silently truncated, and integers above the
exact range of `real64` cannot all be represented.

#### Atomic-data input

The table readers trust counts and indices from text files:

- `src/nlevel_mod.f90:80-105` checks `NLEV` and transition counts against
  maxima, but does not check `nc(k) <= 16`, level indices, negative counts, or
  each read status.
- `src/species_mod.f90:200-221` does not verify
  `nstage <= MAX_ST`, transition index bounds, or `DR n <= 9`.
- `src/cooling_mod.f90:104-113` trusts the term count and rows.

A malformed or newly generated table can therefore cause out-of-bounds access
or leave a partially initialized model.

Recommended change:

- Parse configuration and tables on rank 0 with `iostat`, `iomsg`, line
  numbers, schema checks, and finite/range validation, then broadcast typed
  validated data.
- Use integer namelist fields for integer quantities.
- Centralize validation in one routine and test every rejected input.

### 6.7 [P1] Regression checks are not uniformly enforceable

The repository contains useful scientific cases, but there is no top-level
test runner or dependency manifest. Several scripts print a gate result
without exiting nonzero on failure. Examples include:

- `tests/g1_stromgren/check_g1_stromgren.py`
- `tests/g2a_thermal/check_g2a.py`
- `tests/g4_refine/check_g4.py`
- `tests/pdr/check_pdr.py`

Some scripts skip missing result files and continue successfully. Several
reference calculations intentionally duplicate the same equations and
coefficients as the Fortran code. Those are useful integration checks, but a
shared equation error can appear in both the implementation and its oracle.

Recommended change:

- Provide one test command that builds, runs cases in fresh temporary
  directories, and returns nonzero on every failure or missing output.
- Convert checks to assertions with documented numerical tolerances.
- Separate unit tests, integration tests, conservation tests, and independent
  physics benchmarks.
- Run debug bounds/FPE builds, precise release builds, MPI sizes 1/2/many,
  GNU/Intel, and FITS/HDF5 backends in CI.
- Treat plots and printed diagnostics as artifacts, not pass/fail mechanisms.

### 6.8 [P1] MPI errors and fatal paths are inconsistent

Most MPI calls store `ierr` but do not inspect it. The code uses a mixture of:

- `call MPI_FINALIZE(ierr); stop`
- `error stop`
- plain `stop`
- `MPI_ABORT`

This makes rank-local parser, allocation, or I/O failures difficult to handle
collectively. The program also uses the legacy `use mpi` interface, which
provides less compile-time checking than `mpi_f08`.

`src/ion_peel_mod.f90:393-399` reduces whole image arrays using default-integer
`size(...)` counts. Very large images can exceed the traditional MPI count
range. Tally code already uses chunking in some paths, but this is not a
uniform policy.

Recommended change:

- Add one `mpi_check(ierr, context)` wrapper and check every MPI call.
- Use one collective fatal-error routine with a rank-qualified message.
- Migrate new or modified code to `mpi_f08`.
- Chunk all potentially large collectives or use large-count interfaces where
  available.
- Explicitly free derived communicators and shared-memory windows at shutdown.

### 6.9 [P2] Global state and broad module imports increase coupling

`src/define.f90` is a 900-plus-line public module containing constants, packet
types, global configuration, MPI state, grid-related types, and procedure
pointers. More than fifty source locations import `define`, often without an
`only` list.

Large modules such as `ion_band_mod`, `hdf5io_mod`, `fitsio_mod`, and
`memory_mod_mpi` combine several responsibilities and extensive rank/type
boilerplate. This does not directly make the numerical result wrong, but it
makes hidden dependencies and stale-state defects more likely.

Recommended change:

- Split constants, configuration types, packet types, MPI context, and runtime
  state into private-by-default modules.
- Use `use module, only: ...` consistently.
- Pass a simulation context or smaller explicit state objects to core
  algorithms instead of reading mutable globals.
- Keep pure rate functions and remapping kernels independent of MPI and I/O so
  they can be unit tested.
- Generate repetitive type/rank wrappers or use appropriate modern Fortran
  generic facilities.

### 6.10 [P2] Additional portability and operational issues

- `hdf5_state_type` stores HDF5 handles as `int64` and the source explicitly
  assumes `hid_t`/`hsize_t` are eight bytes. Preserve the HDF5-provided kinds
  inside the backend rather than relying on this ABI assumption.
- `io_detect_format` at `src/iofile_mod.f90:84-95` is case-sensitive, while
  other format parsing is case-insensitive. An uppercase `.H5` can be routed
  to the FITS backend.
- Several output and input filenames have fixed lengths of 128, 192, or 256
  characters, allowing silent truncation.
- Some runtime data files are resolved by trying `../../data` and then `data`,
  making behavior dependent on the current working directory. Use one
  executable-relative or explicitly configured data root.
- Large allocations rarely use `stat=`. A rank-specific allocation failure can
  terminate without a coordinated diagnostic.
- Long runs have no source-visible checkpoint/restart path. A checkpoint
  containing gas state, opacity state, iteration number, RNG state, and
  configuration fingerprint would reduce operational risk.

## 7. Recommended tests

### 7.1 P0 scientific correctness tests

1. **Final state/rate consistency**
   - Force a one-iteration or intentionally unconverged run.
   - Re-transport the written state with the same packet set.
   - Require written rates and dust heating to match the consistency pass.

2. **AMR conservation**
   - Refine and coarsen a grid with density and ionization discontinuities.
   - Compare H, every He stage, electron number, thermal energy, and
     `sum(rhokap V)` before and after.

3. **Live dust normalization**
   - Test `laursen09_live + taumax` and `laursen09_live + tauhomo`.
   - Verify that the first opacity refill preserves the normalization apart
     from the intended survival-factor change.

4. **Two-level limits**
   - Test C II and O I at densities far below and above the critical density.
   - Verify low-density scaling and the high-density Boltzmann population.

### 7.2 P1 scientific robustness tests

5. Verify that toggling `grain_pe` changes PE heating but not independent
   H-impact cooling.
6. For `L_diff/L_packet = 0.1, 0.6, 1.4`, require exact emitted energy and
   compare diffuse spectral moments with the normalized truncated model.
7. Sweep metal-stage ratios from `1e-300` to `1e300` and require finite,
   non-negative fractions summing to one.
8. Supply a synthetic multi-root cooling curve and verify stable-branch
   selection and root-status reporting.
9. In absorption-only boxes and slabs, verify:

   ```text
   source luminosity = escaped + gas absorbed + dust absorbed
   ```

10. Independently vary photon count and frequency-bin count to distinguish
    Monte Carlo error from group-discretization error.

### 7.3 P0/P1 software tests

11. Compile GNU and Intel configurations with `HDF5=0` and `HDF5=1`, in debug
    and precise release modes.
12. Run a real `make -j` or equivalent dependency-aware parallel build from a
    clean out-of-tree directory.
13. Inject HDF5/FITS open, create, write, and close failures and require a
    nonzero program exit without a success message or final-file rename.
14. Fuzz namelist and atomic-table counts, indices, truncated rows, NaNs, and
    overlong paths under bounds checking.
15. Run small deterministic tests at MPI sizes 1, 2, and multiple nodes and
    compare within declared reproducibility tolerances.
16. Make every existing `check_*.py` script fail the process when its gate
    fails or an expected product is missing.

## 8. Recommended implementation order

1. Fix HDF5 error preservation and enforce checked output status.
2. Restore GNU/HDF5-off builds and compiler-specific OpenMP handling.
3. Add the unconditional final radiation/rate consistency pass.
4. Implement conservative AMR remapping with conservation assertions.
5. Persist the live-dust normalization.
6. Correct the C II/O I denominator and decouple H-impact cooling from
   `grain_pe`.
7. Make the regression suite executable and failure-enforcing in CI.
8. Fix diffuse packet energy accounting and energy-weighted sampling.
9. Add log-sum-exp metal fractions and multi-bracket thermal solving.
10. Centralize validated parsing and MPI/I/O error handling.
11. Replace the build system or add correct module dependencies and
    dependency discovery.
12. Gradually split global state and large modules while keeping the new
    scientific regression suite green.

## 9. Conclusion

The central Monte Carlo path-length estimator, analytic direct attenuation,
photo-rate unit conversion, and H/He equilibrium structure are sound under
the reviewed assumptions. The most important scientific risks lie in final
nonlinear-state consistency, non-conservative AMR remapping, lost live-dust
normalization, and the dense-limit fine-structure cooling expression.

From a software-engineering perspective, output reliability and build
portability need immediate attention. The current HDF5 status handling can
hide failures, the GNU non-HDF5 configuration does not compile, and the GNU
link path receives an incompatible OpenMP flag. The existing tests provide
valuable cases, but not all of them enforce failure, so they do not yet form a
reliable safety net for scientific refactoring.

Fixing the P0 items and making the regression suite authoritative would create
a strong base for improving the thermal model, adding physics, and optimizing
parallel performance without sacrificing result integrity.
