# Quasi-Random Photon Launching in MoCHII

## Summary

MoCHII is a promising candidate for randomized quasi-Monte Carlo (RQMC)
sampling at photon emission. The recommended first implementation is an
Owen-scrambled Sobol sequence for the initial stellar packet variables,
while retaining the existing Mersenne Twister generator for dust scattering,
Russian roulette, and other variable-length histories.

For a single internal point source, the launch problem is only
three-dimensional:

1. source-frequency selection;
2. polar direction, represented by `cos(theta)`; and
3. azimuthal direction `phi`.

This is a favorable regime for RQMC. In MoCHII's absorption-only ionizing
transport, the attenuation along a launched ray is integrated analytically.
Consequently, much of the remaining stellar-field noise originates in the
sampling of launch frequency and direction. A low-discrepancy launch set can
therefore reduce noise in global photoionization rates, ionized volume, and
iteration-to-iteration changes without modifying the transport estimator.

The expected gain must nevertheless be measured. AMR cell boundaries,
ionization fronts, and shadow edges make the integrand discontinuous, so the
ideal smooth-integrand QMC convergence rate should not be assumed.

## Current sampling path

The stellar ionizing loop is in `src/ionizing_field_mod.f90` (it was in
`src/main.f90` when this was written). Each MPI rank handles a strided subset
of the global photon indices:

```fortran
do ip = mpar%p_rank+1, par%nphotons, mpar%nproc
   call gen_ion_photon(photon)
   call transport_ion_packet(photon)
end do
```

That stride is now one of the two schedules behind `photon_schedule_mod`
(`par%use_master_slave`); it is still the default and still the same index
sequence. Neither schedule changes the Sobol launch set, since a packet's
coordinates depend only on its global index.

`gen_ion_photon` in `src/ion_band_mod.f90` currently uses independent
Mersenne Twister draws (MT19937-64, `src/random_mt.f90`, one stream per MPI
rank) in the following order:

- a source component through the component-luminosity CDF, when several
  components are active (`multi_src`);
- the emission variables of that component's geometry: for a point source
  the isotropic direction (`mu`, `phi`); for the external rectangular
  surface the entry face (area-weighted CDF), the cosine-weighted incidence
  direction (`mu = sqrt(u)`, `phi`), and two surface coordinates; for the
  external sphere the entry point (two coordinates) and the cosine-weighted
  incidence direction; and
- a frequency **bin** through the discrete band CDF.

The emission geometry itself is fixed by the input (`par%source_geometry`
or the component list), not drawn. Two properties matter for the RQMC
design. First, the launch frequency is a *discrete bin index*
(`photon%inu`, one of `nnu_band` bins); there is no continuous within-bin
frequency, because opacity and rates are evaluated bin by bin. Second, the
present draw order (direction first, frequency second) differs from the
dimension assignment proposed below; this is irrelevant once the uniforms
are passed explicitly, but it means the RQMC path cannot simply reuse the
sequential draws.

For an isotropic internal point source, the direction mapping is

```text
mu  = 2 u_mu - 1
phi = 2 pi u_phi
```

followed by

```text
kx = sqrt(1 - mu^2) cos(phi)
ky = sqrt(1 - mu^2) sin(phi)
kz = mu.
```

This mapping already has the form required by a Sobol point. The principal
code change is therefore to pass a fixed vector of uniforms into the photon
generator instead of requesting them sequentially from the global PRNG.

After launch, `transport_ion_packet` in `src/raytrace_amr.f90` may draw a
variable number of random values for forced interaction, Henyey-Greenstein
scattering, subsequent optical depths, and Russian roulette. These draws
should remain pseudo-random in the first implementation.

## Why randomized QMC is appropriate

Ordinary Monte Carlo integration has an RMS error proportional to

```text
N^(-1/2),
```

independent of dimension under broad conditions. A low-discrepancy sequence
fills the unit hypercube more uniformly than independent random points. For a
smooth low-dimensional integrand, its deterministic integration error can
approach

```text
N^(-1) [log(N)]^d,
```

where `d` is the dimension. Photon transport through an AMR nebula is not
smooth enough to guarantee this rate. Nevertheless, randomized Sobol nets
often retain useful variance reduction for piecewise-smooth transport
problems.

Randomization is important. A deterministic Sobol sequence alone provides no
ordinary statistical error estimate and can form coherent patterns with a
Cartesian or octree grid. Owen scrambling preserves the low-discrepancy
structure while producing independent randomized replicates and an unbiased
RQMC estimator under the usual construction.

MoCHII's analytic absorption estimator makes launch-only RQMC especially
attractive. Given the launch frequency, position, and direction, the direct
absorption contribution is deterministic. RQMC is therefore applied to a
fixed, low-dimensional integral rather than to the entire branching packet
history.

## Recommended coordinate assignment

The meaning of every Sobol dimension must be fixed. Random-number calls must
not consume dimensions conditionally, because that would shift all later
coordinates and destroy the intended low-dimensional projections.

### Single internal point source

Use a three-dimensional scrambled Sobol point:

| Dimension | Variable | Mapping |
|---|---|---|
| 1 | Frequency | Inverse source-luminosity CDF |
| 2 | Polar direction | `mu = 2*u(2) - 1` |
| 3 | Azimuth | `phi = 2*pi*u(3)` |

This should be the first validation case.

### Multiple internal sources

Use separate dimensions for component and frequency:

| Dimension | Variable |
|---|---|
| 1 | Source component through the component-luminosity CDF |
| 2 | Frequency through the selected component's CDF |
| 3 | `mu` |
| 4 | `phi` |

Using one coordinate for both component and conditional frequency can create
undesirable correlations. Separate dimensions preserve better projections.

### External rectangular illumination

The external rectangular source requires additional dimensions:

- entry face;
- two coordinates on the selected face;
- cosine-weighted incidence angle;
- incidence azimuth; and
- frequency.

The dimension assignment should remain identical for every packet even when
different faces are selected. The face-dependent coordinate transformation
should use the same reserved coordinates.

### External spherical illumination

Reserve dimensions for:

- two coordinates defining the entry point on the sphere;
- cosine-weighted incidence angle;
- incidence azimuth; and
- frequency.

Again, the mapping should use a fixed dimension layout.

### Mixed component sets

MoCHII allows point sources and an external field in the same run
(`nsource >= 1` with `ext_intensity > 0`), so the general case needs one
superset layout that every component uses:

| Dimension | Variable |
|---|---|
| 1 | Source component |
| 2 | Frequency bin |
| 3 | Polar / incidence angle (`mu`) |
| 4 | Azimuth (`phi`) |
| 5 | Entry face (rectangular) or first surface coordinate (sphere) |
| 6 | Second surface coordinate |
| 7 | Third surface coordinate (rectangular only) |

A point component simply leaves dimensions 5-7 unused. Unused but fixed
dimensions are harmless: the projections onto the consumed subset remain
those of the scrambled net. What must never happen is a component-dependent
*reassignment* of meanings.

## Frequency sampling

The ionizing-band photon generator already selects the frequency bin with a
monotonic discrete CDF (`ion_cdf`, or `ion_cdf_src` per component), which is
exactly what RQMC needs. The source coordinate `u_nu` should be inverted
through that CDF.

Because the band frequency is discrete, inverting a stratified Sobol
coordinate through the step-function CDF automatically delivers
near-proportional packet counts in every bin (for `N = 2^m` the count in a
bin of probability `p` deviates from `pN` by at most order one). This is a
substantial part of the launch-noise budget: with plain Monte Carlo the bin
counts fluctuate binomially, and the rate integrals inherit that noise
bin by bin. The "specialized design" of explicitly allocating a
deterministic number of packets to each bin is therefore largely subsumed
by the inverse-CDF construction, without any packet-weight bookkeeping.

The alias-sampler caveat applies elsewhere: the multi-wavelength dust/SED
transport (`src/sed_mod.f90`) selects its wavelength bin with an alias
table (`rand_alias_choise`). Alias sampling is efficient for ordinary
random numbers, but its piecewise permutation of the unit interval destroys
the monotonicity that makes a stratified coordinate useful. If RQMC is ever
extended to the SED band, replace the alias lookup by one of:

1. a monotonic inverse CDF;
2. systematic or residual randomized allocation of packet counts to bins; or
3. exact per-bin stratification with a correctly adjusted packet luminosity.

The first option is the simplest and least intrusive. The ionizing-band
generator needs no such change.

## MPI and reproducibility

The Sobol point must be indexed by the global photon number `ip`, not by a
rank-local counter. A random-access interface should conceptually provide

```text
u(:) = sobol_scrambled(global_photon_index, scramble_key).
```

This has several benefits:

- the same physical packet set is generated for different MPI task counts;
- static, dynamic, or strided work assignment does not change the sample;
- restart and debugging are simpler; and
- the launch sequence is independent of each rank's Mersenne Twister state.

It is preferable to use a counter-based or random-access Sobol implementation
rather than advancing a separate stateful generator on each rank. Assigning
one contiguous or strided unsynchronized Sobol subsequence per rank can damage
low-dimensional uniformity or make results depend on MPI decomposition.

One caveat: the launch *set* becomes task-count independent, but any
post-launch draw (forced interaction, Henyey-Greenstein scattering, Russian
roulette) still comes from the rank-local Mersenne Twister, whose stream
depends on which packets a rank handles. Bitwise MPI-count independence of
the full run is therefore obtained only in absorption-only transport, where
the post-launch path is deterministic. With dust scattering enabled the
launch set is still reproducible, but the scattered histories are not.

Sobol nets have their cleanest balance properties when the total photon count
is a power of two. MoCHII should not require this, but documentation can
recommend `nphotons = 2^m` for RQMC production runs.

## Treatment across nonlinear iterations

Using the same scrambled launch set at every ionization/thermal iteration is
a form of common random numbers. It correlates the radiation estimates before
and after an opacity update, often reducing noise in their difference. This
can improve the stability of convergence measures and make front motion easier
to distinguish from sampling noise.

The recommended policy is:

1. select one scramble key at the beginning of a nonlinear solve;
2. keep that key fixed through all opacity/ionization/thermal iterations;
3. use a different scramble key for an independent replicate; and
4. validate the final solution across several independent scrambles.

A fixed scrambled set can also create a repeatable quadrature imprint that the
nonlinear solution adapts to. It may therefore make iteration differences look
quieter than the true sampling uncertainty. Convergence under one scramble is
not a substitute for replicate-to-replicate agreement.

For production uncertainty estimates, use at least four independent scrambles;
eight are preferable when cost permits. The scatter among replicate results
provides the RQMC error estimate.

## What should remain pseudo-random initially

### Dust scattering

Dust histories consume a variable number of uniforms:

- forced-first-interaction optical depth;
- Russian-roulette survival;
- HG polar angle;
- scattering azimuth; and
- the next free-flight optical depth.

Mapping these calls onto one sequential Sobol stream would cause the meaning
of a dimension to depend on the packet's previous history. Later dimensions
would also be high-dimensional and less effective. The first implementation
should therefore switch back to the existing Mersenne Twister immediately
after launch.

A later experiment could reserve a fixed block of dimensions for each
scattering order, but it would need a maximum-order policy and careful tests.
The expected benefit is smaller than for the initial direct field.

### Diffuse emission

A diffuse recombination packet requires more coordinates:

- emitting leaf;
- three coordinates within the leaf;
- direction;
- emission channel; and
- photon energy within the selected channel.

The emitting-leaf CDF is strongly discontinuous and changes as the gas state
iterates. RQMC may still help, but the problem is higher-dimensional and less
stable than stellar launch. Diffuse RQMC should be considered only after the
stellar implementation has a clear benchmarked benefit.

### Stochastic dust emission internals

SEDust and any rejection or variable-iteration samplers should remain
unchanged. The RQMC proposal concerns the MoCHII transport packet launch, not
a global replacement of every PRNG call.

## Expected benefits and limitations

The largest gains are expected for:

- direct stellar photoionization rates;
- integrated ionized volume and effective Strömgren radius;
- volume-integrated convergence measures;
- smooth radial profiles;
- global absorbed gas and dust luminosities; and
- early iterations run with a modest photon budget.

Smaller or inconsistent gains are expected for:

- the maximum error in a single I-front cell;
- narrow shadow boundaries and clump skins;
- highly localized line emissivities;
- high-resolution image pixels;
- scattering-dominated paths; and
- diffuse-field-dominated cells.

Low discrepancy does not eliminate physical discretization errors. Frequency
binning, AMR resolution, atomic-data fits, and nonlinear stopping tolerances
remain separate error sources.

No universal speedup factor should be assumed. A useful reduction in variance
can be reported as an effective photon-number gain: the number of ordinary
Monte Carlo packets required to reach the RMS error obtained by a given RQMC
run.

## Proposed implementation interface

A minimal user-facing design could add:

```fortran
par%launch_sequence = 'random'       ! 'random' or 'sobol'
par%qmc_scramble     = .true.
par%qmc_seed         = 12345
```

The default must remain `random` until the new path has passed all gates.

Internally, avoid replacing `rand_number()` globally. Instead, provide an
explicit launch sample:

```fortran
call launch_uniforms(ip, iteration_key, u)
call gen_ion_photon(photon, u)
```

or a dedicated QMC photon generator. Explicit data flow makes the dimension
assignment auditable and prevents transport routines from accidentally
consuming Sobol dimensions.

The scramble key should be reproducible from the input seed. If a fixed
scramble is used through the nonlinear solve, `iteration_key` should not
change between iterations.

### Implementation notes

No external library is needed; the whole generator is a few hundred lines
of self-contained Fortran.

- **Sobol point by index.** The unscrambled coordinate is
  `x_d(i) = XOR over set bits k of i of v_d(k)`, where `v_d(:)` are the
  direction numbers of dimension `d`. This is a random-access formula:
  cost O(popcount(i)) integer XORs, no stored state, no sequential
  dependence. Direction numbers for the at most eight dimensions needed
  here come from the Joe & Kuo (2008) table; the handful of required
  primitive polynomials and initial numbers can be embedded as constants.
- **Owen scrambling by hashing.** Full tree-based Owen scrambling is
  unnecessary; the practical standard is nested-uniform scrambling with a
  hash (Laine & Karras 2011; Burley 2020): reverse the bits of `x_d`,
  apply a keyed avalanche hash that lets each bit affect only
  less-significant bits, and reverse back. A distinct key per dimension is
  derived from `(qmc_seed, dimension, replicate)`. This is O(1) per
  coordinate and preserves the net's equidistribution while decorrelating
  replicates.
- **Integer-to-uniform mapping.** Map the scrambled `b`-bit integer `k` to
  `u = (k + 0.5) * 2^-b`, which keeps every coordinate strictly inside
  (0,1); the launch mappings themselves use no logarithms, but the open
  interval costs nothing and removes a failure class for later extensions.
- **Cost.** A few tens of integer operations per packet, negligible
  against ray traversal.

## Validation plan

### Gate Q0: unit-hypercube projections

For `N = 2^m`, inspect one- and two-dimensional projections of the generated
coordinates. Verify:

- every coordinate lies strictly inside the interval required by logarithms;
- source and frequency probabilities agree with their target CDFs;
- the angular distribution is isotropic;
- results are identical for different MPI task counts; and
- independent scramble keys produce different but individually balanced nets.

### Gate Q1: analytic absorption test

Use the existing fixed-state H/He attenuation gate (`tests/g0_gamma`) with a
single point source and no scattering. At each photon count, compare ordinary
Monte Carlo and scrambled Sobol replicates for:

- RMS error in `Gamma_HI` and `Gamma_HeI`;
- mean bias;
- radial-profile error;
- maximum and volume-weighted errors; and
- wall time.

This is the cleanest test because the exact reference is already available
and the post-launch path is deterministic.

### Gate Q2: Strömgren iteration

Run the same physical setup (`tests/g1_stromgren`; the `tests/uni_dda`
namelist-sphere variants are convenient for fast car-grid runs) with
identical packet counts and several seeds or scrambles. Compare:

- effective ionized radius;
- ionized volume;
- radial `x_HI`, `x_HeI`, and `x_HeII` profiles;
- cell-max and volume-integrated convergence histories; and
- the iteration number at which the stopping criterion is reached.

Keep one scramble fixed across iterations in the primary RQMC test. Add a
control in which the scramble changes every iteration to quantify the value
of common random numbers.

### Gate Q3: AMR and clumpy medium

Use an I-front-refined sphere (`tests/g4_refine`) and one clumpy/TNG test
(`tests/g4_tng`). Measure whether RQMC reduces integrated noise without
introducing visible angular spokes, grid-aligned structures, or persistent
front corrugation.

### Gate Q4: dust scattering control

Use QMC only at launch and the Mersenne Twister after the forced first event
(`tests/d_dusty`). Confirm that the direct field improves and that scattered
luminosities remain unbiased. Compare total gas absorption, dust absorption,
and escaped energy. Peel-off images need no special treatment: the direct
peel is deterministic given the launch variables, so it inherits the RQMC
launch set unchanged, and the scattered peel rides on the Mersenne Twister
history as before.

### Acceptance criteria

Adopt the QMC option if it satisfies all of the following:

1. no statistically significant bias relative to ordinary Monte Carlo;
2. MPI-count-independent launch samples and integrated results;
3. lower replicate RMS for at least the direct rates and ionized volume;
4. no coherent angular or grid-aligned artifacts; and
5. negligible sequence-generation overhead relative to ray traversal.

## Recommended development order

1. Implement random-access Owen-scrambled Sobol points indexed by global
   photon ID.
2. Apply them only to a single internal point source: frequency, `mu`, and
   `phi`.
3. Pass the analytic attenuation and Strömgren gates.
4. Add multiple internal sources with a fixed four-dimensional layout.
5. Add external rectangular and spherical launch mappings.
6. Evaluate inverse-CDF sampling in the SED path.
7. Consider diffuse packets only if the stellar results justify the added
   complexity.
8. Do not extend QMC into variable-length scattering histories without a
   separate design and validation campaign.

## Conclusion

Randomized quasi-Monte Carlo launching is likely to reduce MoCHII's direct
stellar-field noise, especially in absorption-dominated calculations where
the post-launch estimator is analytic. The lowest-risk design is deliberately
hybrid: scrambled Sobol coordinates for the fixed-dimensional initial launch,
followed by the existing Mersenne Twister for all conditional transport
events. Global photon indexing, fixed dimension semantics, a fixed scramble
during nonlinear iteration, and independent scrambled replicates are the key
requirements for correctness and useful uncertainty estimates.

The method should be introduced as an optional, gated variance-reduction
feature rather than as a global replacement for pseudo-random sampling.
