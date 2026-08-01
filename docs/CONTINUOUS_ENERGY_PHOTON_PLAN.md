# Plan: continuous-energy ionizing packets and exact-energy cross sections

- Date: 2026-07-28
- Target: the ionizing-band source, transport, and rate-estimator path
- Status: continuous slice released; both-band continuous unification + grouped removal planned (Section 19, aligned to MoCafe v2.00 1473f44)
- Primary modules: `src/ion_band_mod.f90`, `src/define.f90`,
  `src/gas_opacity_mod.f90`, `src/raytrace_amr.f90`,
  `src/gas_rates_mod.f90`, `src/species_mod.f90`, and
  `src/diffuse_mod.f90`

Implementation progress as of 2026-07-28:

- Safe sequence 16.1 complete: grouped HII20/HII40 baselines and the
  C II/He I threshold-window production-cross-section gate are frozen.
- Safe sequence 16.2 complete: grouped-default mode control, fail-fast
  continuous-mode guard, packet `energy_eV` plumbing, and rates-file metadata
  are implemented without changing grouped physics.
- Safe sequence 16.3 complete: AMR and DDA direct/scattered walks use the
  common `score_ion_path` interface. Opt-in grouped H/He and metal shadow
  estimators agree with the legacy `jt_ion` reconstruction at approximately
  `3e-15` or better; shadow-off and shadow-on solver arrays are bitwise
  identical.
- Safe sequence 16.4 complete: each transported packet receives an H/He,
  flattened-metal, and dust physics cache at its grouped bin-center energy.
  Metal stage fractions are cached per leaf and refreshed once per gas-state
  opacity update. H/He-only, metal, dust-scattering, AMR, DDA, and
  re-refinement dynamic-opacity shadows agree with `kap_ion` at approximately
  `3e-16` or better.
- Safe sequence 16.5 complete: an isolated continuous-energy sampler provides
  analytic inversion of piecewise-linear tabulated spectra, an adaptive
  band-limited Planck CDF shared by random and Sobol uniforms, and the
  five-uniform full-Planck Barnett--Canfield reference. Random/Sobol moments,
  quantiles, threshold fractions, the C II/He I window, multi-source,
  external, diffuse-stream, and bitwise 1/2/3/5/8-rank tests pass.
- Safe sequence 16.6 complete: both pseudorandom and Sobol diffuse launchers
  retain their sampled `eph` in continuous mode and explicitly select the
  historical bin center in grouped mode. A standalone policy gate and an
  actual grouped diffuse-packet smoke pass.
- Safe sequence 16.7 complete: the guarded single-point H/He-only slice uses
  continuous Planck or tabulated source energies, preserves diffuse energies,
  evaluates exact-energy H/He opacity, accumulates direct exact-energy
  ionization/heating estimators, and feeds those estimators to the solver.
  Planck, threshold-structured table, random/Sobol, diffuse, iterative-solver,
  1/3-rank, emitted/absorbed/escaped energy-closure, and 8/32
  diagnostic-bin tests pass.
- Safe sequence 16.8 step 17 complete: trace-metal ionization and heating
  estimators consume each packet's exact-energy cross-section cache while
  `ion_metal_abs=.false.` keeps metals out of transport opacity. Metal
  `Gamma`, heating, and stage fractions are bitwise identical with 8 and 32
  diagnostic bins.
- Safe sequence 16.8 step 18 complete: the flattened transition cache is
  threshold-sorted once and packet cross-section evaluation stops at the
  first inaccessible transition. The operation benchmark predicts a
  61--69% reduction for HII20/HII40 and a 56--64% reduction with every
  registry element enabled, without stochastic estimator noise.
- Safe sequence 16.8 step 19 complete: continuous transport combines the
  packet's exact-energy metal cross sections with the iteration-refreshed
  per-leaf stage cache. Iterated 8/32-bin states and rates are bitwise
  identical, 1/3-rank results agree, and an HII20 AMR/thermal smoke passes.
- Safe sequence 16.9 step 20 complete: `metal_ne` and `metal_heat` consume
  exact-energy metal ionization/heating arrays in the continuous thermal
  iteration. Full-coupling 8/32-bin states are bitwise identical and 1/3-rank
  outputs agree.
- Safe sequence 16.9 step 21 complete: 8-rank, `2^25`-packet Sobol
  production calculations of HII20 and HII40 converged in 27 and 38 gas
  iterations, respectively, and both completed their final consistency pass.
  The production gate passes energy closure, metadata, C III-tail recovery,
  Cloudy radial-profile non-regression, and temperature non-regression.
- Safe sequence 16.11 step 29 complete: full HII20/HII40 production runs on
  68 MPI ranks with `nnu_ion = 8, 16, 32, 64, 128` used the same `2^25`
  Sobol packet history.  All solver/state/rate arrays were bitwise identical
  across diagnostic resolutions; only `J_nu`, `E_bin`, and `dE_bin` changed.
- Safe sequence 16.11 RQMC ensemble portion complete: four independent
  scrambled-Sobol seeds passed for both production models.
- Safe sequence 16.11 node-local MPI/performance portion complete:
  HII20/HII40 100,000-packet, three-iteration continuous-metal smoke runs
  passed at 1, 2, 3, 5, 8, and 68 ranks.  The 68-rank tests used the 72-thread
  node while reserving four hardware threads.  Multi-node execution is also
  verified on lart2/lart3/lart4.
- Safe sequence 16.11 step 32 complete: the paper manuscript, its HII20/HII40
  tables, Figures 14--15, captions, and continuous-energy discussion now use
  the accepted production outputs.  The paper records `2^25` packets, Sobol
  seed 20260728, continuous Planck-CDF energy sampling, exact-energy cross
  sections, and the diagnostic-only role of `nnu_ion`.  Generated outputs are
  organized under `results/continuous_energy/`, with the paper reading the
  `production/` subset.
- Next: update `README.md` immediately before the single final push, then
  commit and push the completed release-validation sequence.

## 1. Required physical model

The input spectrum is necessarily represented numerically, either by an
analytic function evaluated on a quadrature grid or by an interpolated
tabulated spectrum. That representation must define a continuous probability
density from which each packet energy is sampled. The transport physics must
then use the sampled energy itself.

For MoCHII's equal-energy packets, the source probability density is

```text
p(E) = L_E(E) / integral L_E(E) dE .
```

After sampling `E_packet`, all relevant quantities must be evaluated at that
energy:

```text
sigma_s       = sigma_s(E_packet)
kappa_cell    = kappa(E_packet, cell)
Gamma tally   proportional to sigma_s(E_packet) / E_packet
heating tally proportional to sigma_s(E_packet)
                         * (1 - E_threshold,s / E_packet)
```

Energy bins may remain for spectrum representation, acceleration, diagnostic
histograms, and output. They must not determine the physical energy,
cross-section, opacity, ionization rate, or heating rate of a packet.

## 2. Why the current implementation must change

The source spectrum is currently integrated accurately enough within each
energy bin: `ion_band_mod` uses 32 midpoint samples per bin to construct the
bin luminosities. A packet then samples only a bin index, however, and
`photon_type` stores only `photon%inu`. It does not store the sampled physical
energy.

Downstream code consequently treats every packet in a bin as if its energy
were

```text
ion_e(inu) = sqrt(E_edge(inu) * E_edge(inu+1)).
```

The following calculations all use this representative energy:

- H I, He I, and He II transport opacity;
- optional metal transport opacity;
- dust absorption, scattering, and asymmetry;
- H/He photoionization and photoheating rates;
- metal photoionization and photoheating rates;
- secondary ionization;
- conversion from energy luminosity to photon number.

Diffuse recombination photons initially sample a continuous energy, but that
energy is immediately converted to `photon%inu` and discarded. They therefore
have the same grouped-energy limitation during transport and rate evaluation.

This design causes the C II/He I failure seen in the HII20 and HII40
benchmarks. The C II threshold is 24.383 eV and the He I threshold is
24.587 eV. The aligned grid merges the thresholds and represents the
23.330--24.587 eV group by 23.950 eV. Stellar luminosity between 24.383 and
24.587 eV is present in the group integral, but
`sigma_CII(23.950 eV) = 0`, so those photons cannot produce C III.

Separating those two edges would repair this particular threshold window, but
would remain a grouped-energy approximation. It is a useful short-term
diagnostic, not the fundamental solution.

## 3. Target architecture

The physical and diagnostic paths should be separated:

```text
source spectrum
    |
    +-- continuous inverse-CDF sampling --> photon%energy_eV
                                                |
                                                +--> exact sigma(E)
                                                +--> exact kappa(E, cell)
                                                +--> direct rate estimators
                                                |
                                                +--> output-bin lookup
                                                     photon%inu
                                                     diagnostic J_E only
```

The solver must consume the direct exact-energy estimators. A binned `J_E`
tally may still be written, but it must not be converted back into solver
rates with `sigma(E_bin)`.

## 4. Migration and compatibility

Introduce an explicit transition switch while the implementation is being
validated:

```fortran
par%ion_energy_mode = 'grouped'
par%ion_energy_mode = 'continuous'
```

The grouped path initially remains unchanged and supplies a regression
reference. The continuous path is developed beside it. After all acceptance
tests pass, make `continuous` the default. Retain `grouped` temporarily for
reproduction of historical runs, then decide separately whether to remove it.

Keep `grouped` as the default throughout development. A partially implemented
continuous mode must fail explicitly when it encounters an unsupported
feature:

```text
continuous mode + unimplemented metal opacity       -> fatal error
continuous mode + unimplemented secondary ionization -> fatal error
continuous mode + unimplemented dust scattering      -> fatal error
```

Never fall back silently to a bin-center calculation inside continuous mode.
Such a fallback would produce a hybrid result that appears to be
continuous-energy but is physically inconsistent.

The physical cutover must be atomic within each supported feature set. For
the initial H/He-only vertical slice, continuous source energy, exact-energy
opacity, direct ionization estimators, and direct heating estimators must
become active together. Do not expose a mode in which transport uses the
sampled energy while the gas solver still uses rates reconstructed from bin
centers.

Every rates file must record the selected energy mode and the source sampler
tolerance. Historical and new results must therefore be distinguishable
without consulting the input file.

## 5. Phase 0: freeze baselines and add focused tests

Before changing the packet representation, preserve the following baseline
outputs:

- HII20 and HII40 radial H, He, C, N, O, and Ne fractions;
- electron temperature and density;
- H/He photoionization and heating rates;
- metal ionization rates if exposed by a temporary diagnostic;
- total absorbed and escaped luminosity;
- major line luminosities;
- iteration history and convergence state;
- runtime, peak memory, and MPI rank count.

Add a threshold-window unit test with monochromatic packets:

| packet energy | C II ionization | He I ionization |
|---:|:---:|:---:|
| 24.38 eV | no | no |
| 24.40 eV | yes | no |
| 24.58 eV | yes | no |
| 24.60 eV | yes | yes |

Also require a 24.40 eV packet to traverse neutral-helium gas without He I
absorption while still contributing to the C II photoionization estimator.
This test directly guards the physical window that motivated the change.

## 6. Phase 1: store a physical packet energy

Add an ionizing-packet energy to `photon_type`:

```fortran
real(kind=wp) :: energy_eV = 0.0_wp
integer       :: inu       = 0
```

The fields have different roles:

- `energy_eV` is authoritative for all transport and rate physics;
- `inu` is derived from `energy_eV` and is used only for diagnostic spectral
  tallies and output.

Add assertions in continuous mode:

```text
energy_eV > 0
eion_min <= energy_eV <= eion_max
inu == ion_bin_of(energy_eV), when a diagnostic bin is requested
```

No continuous-mode physics routine may accept only `inu` where an energy is
required.

## 7. Phase 2: implement continuous source sampling

### 7.1 Equal-energy sampling

Preserve the equal-energy packet model:

```text
L_packet = L_band / N_packets .
```

Packet energies are sampled from the energy-luminosity density `L_E`, not
from `L_E/E`. Photon-number weighting enters the rate estimator through
division by the sampled photon energy.

### 7.2 Analytic spectra

For a Planck source, build a monotone continuous CDF with adaptive
quadrature. All active H, He, and metal ionization thresholds must be inserted
as integration knots so that discontinuities in downstream cross-sections do
not coincide with an unresolved sampler interval.

Recommended requirements:

- configurable relative integral tolerance, default `1e-8`;
- monotone CDF with exactly normalized final value;
- robust inversion by bracketed interpolation;
- no dependence on `nnu_ion` or the diagnostic output grid;
- the same uniform variate produces the same energy independent of MPI rank.

#### 7.2.1 Direct full-Planck sampler for pseudorandom mode

The pseudorandom path may also use the exact Barnett--Canfield series sampler
when the requested source is the full, untruncated Planck distribution. For
equal-energy packets, define

```text
x = E / (k_B T)
p(x) = (15 / pi^4) x^3 / [exp(x) - 1].
```

Using

```text
x^3 / [exp(x) - 1] = sum_(L=1..infinity) x^3 exp(-L x),
```

first draw one independent uniform `xi0` and choose

```text
L = min { l : sum_(j=1..l) j^-4 >= xi0 * pi^4 / 90 }.
```

Then draw four further independent uniforms and set

```text
x = -log(xi1 * xi2 * xi3 * xi4) / L
E_packet = k_B T * x.
```

This is sometimes described as a four-random-number Planck sampler because
four uniforms appear in the final product, but it consumes five independent
uniforms in total: one for the discrete mixture index `L` and four for the
shape-4 gamma variate. The `L` uniform must not be reused in the product.

This algorithm samples the energy Planck distribution `x^3/[exp(x)-1]`
required by equal-energy packets. A three-factor product would instead be
associated with the photon-number distribution and is not correct for
MoCHII's current packet weighting.

Implement this as an internal pseudorandom Planck-sampler algorithm, not as a
second top-level launch-sequence control. Validate it independently against
analytic Planck moments and the inverse-CDF sampler.

#### 7.2.2 Restricted ionizing band

The Barnett--Canfield formula above samples the full interval
`0 < E < infinity`, whereas MoCHII normally transports a conditional band
such as 13.598--100 eV. Generating the full Planck distribution and rejecting
out-of-band energies is exact but can be very inefficient:

- for the 20 kK HII20 source, only about 4.2% of the bolometric energy lies
  in the ionizing band, requiring about 24 full-Planck draws per accepted
  packet on average;
- for the 40 kK HII40 source, the corresponding fraction is about 41%.

The production sampler for a band-limited Planck source should therefore use
the conditional CDF

```text
F_band(E) =
    integral_(Emin..E) B_E(E',T) dE'
    / integral_(Emin..Emax) B_E(E',T) dE'.
```

This requires one uniform, has no rejection, and is shared naturally by the
pseudorandom and Sobol paths. The full-Planck series sampler remains useful
as:

- an exact unit-test reference;
- an efficient production path if a future configuration transports the
  complete Planck spectrum;
- an independent statistical cross-check of the conditional inverse CDF
  after applying a band cut.

A truncated Barnett--Canfield mixture could be constructed using
band-integrated mixture weights and truncated shape-4 gamma variates, but it
is not the initial implementation. It is more complex than the conditional
inverse CDF and provides no clear benefit for the standard MoCHII bands.

### 7.3 Tabulated spectra

First convert the input to ascending `(E, L_E)` points as the current code
does. Define the continuous spectrum by the documented interpolation rule.
For the present linear interpolation, integrate each interval analytically.
After the CDF selects an interval, invert the integral of the local linear
function analytically, using a numerically stable linear limit when its slope
is small.

Sampling is then exact for the chosen interpolation model. The resolution of
the input table remains a property of the supplied spectrum, but transport
bin boundaries do not alter the sampled distribution.

### 7.4 Multiple sources and external fields

Keep source-component selection proportional to each component's integrated
band luminosity. After selecting a component, sample energy from that
component's continuous CDF. Each component therefore retains its own spectral
shape without a shared grouped-energy approximation.

### 7.5 Random and QMC paths

Continuous photon-energy sampling must support both of MoCHII's existing
launch sequences:

```fortran
par%launch_sequence = 'random'  ! pseudorandom energy and launch coordinates
par%launch_sequence = 'sobol'   ! scrambled-Sobol energy and launch coordinates
par%qmc_seed        = 12345
```

Do not introduce a second, overlapping energy-sampler switch. Extend the
meaning of the existing `launch_sequence='sobol'` option so that frequency
means the continuously sampled physical energy, rather than only a discrete
energy-bin index.

Both paths must call exactly the same inverse-CDF routine:

```fortran
E_packet = source_energy_from_cdf(source_id, u_energy)
```

Only the origin of `u_energy` differs:

- `random`: a Mersenne-Twister uniform;
- `sobol`: the existing Owen-scrambled Sobol frequency coordinate.

This common inverse-CDF rule applies to all band-limited spectra. The optional
five-uniform Barnett--Canfield implementation in Section 7.2.1 is confined to
the pseudorandom full-Planck case and its unit tests; it must not replace the
one-coordinate conditional inverse CDF for the normal band-limited or Sobol
paths.

Keep the current QMC coordinate assignments:

| launch configuration | energy coordinate |
|---|---|
| single internal point source | `u(1)` |
| slab source | `u(1)` |
| single external source | `u(2)` |
| multiple or mixed source components | `u(2)`, after `u(1)` selects the component |
| diffuse recombination stream | the existing channel/energy coordinates `u(7:9)` |

The source-component CDF and continuous energy CDF must remain separate.
For a multi-source launch, `u(1)` first selects a source in proportion to its
integrated luminosity; `u(2)` is then inverted through that source's own
continuous spectral CDF.

No additional QMC dimension is required. Reusing the current frequency
coordinate preserves the launch dimension, the global-packet-index mapping,
and MPI-task-count independence. It also prevents a source-sampling change
from shifting the Sobol dimensions assigned to direction and position.

The QMC transform must be monotone: a larger frequency coordinate must never
produce a smaller sampled energy for one source. Exact CDF inversion therefore
preserves the one-dimensional low-discrepancy property in energy. Do not use
rejection sampling in the Sobol path, because a variable number of uniforms
would destroy the fixed coordinate semantics and MPI-independent mapping.

The diffuse QMC path already uses a separate scrambled stream. Preserve that
separation. Store the continuously sampled diffuse energy in
`photon%energy_eV`; do not reduce it to `photon%inu` before transport.

Write the effective energy-sampling sequence to the output metadata, so a
result records whether its continuous energies came from pseudorandom or
scrambled-Sobol uniforms.

## 8. Phase 3: preserve diffuse-photon energies

`diffuse_mod` already samples physical continuum and line energies. Store the
result instead of discarding it:

```fortran
photon%energy_eV = eph
photon%inu       = ion_bin_of(eph)  ! diagnostic only
```

Apply this to:

- H I ground-state recombination continuum;
- He I ground-state recombination continuum;
- He II ground-state recombination continuum;
- fixed-energy He I cascade photons;
- the He I two-photon continuum.

Dust scattering is elastic in the current model, so it changes direction but
does not change `energy_eV`.

## 9. Phase 4: evaluate packet physics at the sampled energy

Cross-sections depend on packet energy but not on the cell. Compute them once
per packet and reuse them during every cell crossing.

A transport-local cache should contain at least:

```fortran
type ion_packet_physics_type
   real(wp) :: energy_eV
   real(wp) :: sigma_HI
   real(wp) :: sigma_HeI
   real(wp) :: sigma_HeII
   real(wp) :: sigma_metal(MAX_METAL_TRANSITIONS)
   real(wp) :: dust_abs
   real(wp) :: dust_sca
   real(wp) :: dust_g
end type ion_packet_physics_type
```

The exact layout may be flattened using the species registry, but it must
avoid re-evaluating every analytic cross-section in every crossed cell.

The cache is filled with:

```text
sigma_HI(E_packet)
sigma_HeI(E_packet)
sigma_HeII(E_packet)
species_sigma(element, transition, E_packet)
dust properties at lambda = 1.23984 / E_packet [micron]
```

Threshold behavior must come only from the cross-section functions evaluated
at `E_packet`, never from the diagnostic bin index.

## 10. Phase 5: replace grouped opacity in continuous mode

For each packet and leaf, calculate

```text
kappa(E, leaf) =
    n_H * [
        x_HI * sigma_HI(E)
      + A_He * (x_HeI * sigma_HeI(E) + x_HeII * sigma_HeII(E))
      + sum_elements A_Z sum_stages x_Zi * sigma_Zi(E)
    ] * distance2cm
  + kappa_dust_abs(E, leaf)
  + kappa_dust_sca(E, leaf), when scattering is enabled.
```

The packet cache supplies the cross-sections; the leaf supplies densities and
ion fractions. This is inexpensive for H/He. For metals, flatten the active
transitions once during species setup and use a compact dot product.

Before exact metal opacity is enabled, add a per-leaf metal-stage population
cache that is refreshed once after each gas-state update:

```text
n_CI(leaf), n_CII(leaf), n_CIII(leaf), ...
n_NI(leaf), n_NII(leaf), ...
```

Do not call `species_fractions` during every packet-cell crossing. Re-solving
the metal cascade inside the transport walk would multiply the most expensive
part of the ionization calculation by the number of path segments. Exact
metal opacity must combine packet-cached cross-sections with these cached
per-leaf stage populations.

The existing analytic direct edge walk remains valid. It requires a constant
opacity inside one leaf, not a grouped energy. Replace

```fortran
kap = kap_ion(inu, leaf)
```

with a continuous-mode call that combines the cached packet cross-sections
with the leaf state.

Keep `kap_ion(nbin,nleaf)` only for grouped compatibility during migration.
Once continuous mode is established, removing it can recover substantial
memory.

During migration, implement the dynamic opacity function first in shadow
mode. Evaluate it at `E = ion_e(inu)` and compare it with the existing
`kap_ion(inu,leaf)` array before allowing it to transport continuous-energy
packets.

## 11. Phase 6: accumulate exact-energy rate estimators

Exact-energy opacity alone is insufficient. If the solver later derives rates
from `J_bin * sigma(E_bin)`, bin dependence returns. The solver rates must be
accumulated directly from packet path estimators.

The path-scoring infrastructure must be landed before continuous opacity is
activated. Initially run the new estimator in grouped shadow mode with
`E_packet = ion_e(inu)` and compare its rates with the existing
`J_bin * sigma(E_bin)` post-processing. This separates scoring and MPI
reduction errors from the later physical change in packet energy.

For a path contribution `P = L_packet * weight * effective_path_length` in a
leaf of volume `V`, accumulate for absorber `s`:

```text
dGamma_s = P / (V * distance2cm^2)
           * sigma_s(E_packet) / (E_packet * eV_to_erg)

dHeat_s  = P / (V * distance2cm^2)
           * sigma_s(E_packet)
           * max(1 - E_threshold,s / E_packet, 0)
```

Use the same effective path contribution already generated by the analytic
zero-variance direct walk:

```text
P = L_packet * weight
    * [exp(-tau_in) - exp(-tau_out)] / kappa
```

For explicitly followed scattered paths, use the existing Lucy path-length
contribution. Route both through one rate-tally helper so direct and scattered
transport cannot diverge.

Accumulate:

- H I, He I, and He II `Gamma`;
- H I, He I, and He II heating;
- each active metal transition's `Gamma` and heating;
- hard-photoelectron terms required by secondary ionization;
- dust absorbed-energy estimators;
- FUV estimators such as `G0`, when the packet is in the relevant range.

Reduce these arrays over MPI in the same deterministic ownership pattern as
the current radiation tally. The gas and thermal solvers must consume these
direct estimators in continuous mode.

## 12. Phase 7: make spectral bins diagnostic-only

Retain an optional `J_E` histogram:

```text
diagnostic_bin = ion_bin_of(E_packet)
J_tally(diagnostic_bin, leaf) += path contribution
```

This output remains useful for plotting and comparison, but it must not feed
the following continuous-mode calculations:

- gas opacity;
- H/He ionization;
- metal ionization;
- photoheating;
- secondary ionization;
- dust heating;
- thermal balance.

Rename future input controls to make the distinction explicit, for example:

```fortran
par%nnu_output
par%source_cdf_tol
```

The legacy `nnu_ion` name may be retained temporarily for input
compatibility, but its continuous-mode meaning must be documented as output
resolution only.

## 13. Phase 8: output and restart metadata

Write at least the following metadata:

```text
IONEMODE = grouped | continuous
CDFRTOL  = source-sampler relative tolerance
NUBINUSE = diagnostic_only | solver_and_diagnostic
ENRGSAMP = random | sobol
QMCSEED  = scramble seed when ENRGSAMP = sobol
```

If packet histories or restart state are ever serialized, store
`energy_eV`. Reconstructing a packet energy from `inu` is not valid in
continuous mode.

The rates output should distinguish direct exact-energy solver rates from
diagnostic binned `J_E`.

## 14. Validation

### 14.1 Source sampler

For Planck spectra at 20 kK and 40 kK, compare sampled and quadrature values
for:

- total band luminosity;
- mean packet energy;
- median and selected energy quantiles;
- luminosity fractions above every active threshold;
- luminosity in 24.383--24.587 eV.

Run both pseudorandom and QMC sampling. Statistical errors must scale as
expected, and QMC results must be independent of MPI task count at launch.

Test the full-Planck Barnett--Canfield implementation separately:

- confirm that it consumes five independent uniforms per variate;
- compare its sampled mean, variance, and selected quantiles of
  `x = E/(k_B T)` with analytic Planck values;
- compare a large sample with the full-range inverse-CDF result;
- apply the HII20 and HII40 energy cuts and confirm that accepted samples
  follow the same conditional distribution as the band-limited inverse CDF;
- measure the rejection acceptance fractions and keep rejection out of the
  HII20 production path.

For the QMC path specifically:

- compare the sampled continuous-energy quantiles with the analytic CDF;
- verify luminosity fractions above every ionization threshold;
- verify the 24.383--24.587 eV luminosity fraction for both 20 kK and 40 kK
  Planck spectra;
- assert that packet index `p` and `qmc_seed` produce the same energy at 1,
  2, 3, 5, and 8 MPI ranks;
- assert that changing only the diagnostic energy bins does not change the
  sampled Sobol energies;
- compare several independent scramble seeds, rather than treating one Sobol
  realization as an uncertainty estimate;
- compare random and scrambled-Sobol ensemble means to detect an inversion or
  coordinate-mapping bias.

For tabulated spectra, test constant, linear, sharply broken, and
threshold-adjacent spectra against analytic interval integrals.

### 14.2 Cross-section and opacity tests

Test cross-sections immediately below, at, and above every H/He threshold and
the C II/He I pair. Compare one-cell optical depths against analytic values:

```text
tau(E) = n_s * sigma_s(E) * path_length.
```

Test mixtures of H I, He I, and C II to prove that a 24.40 eV packet sees C II
but not He I opacity.

### 14.3 Rate-estimator tests

Use a uniform optically thin cell and monochromatic source. Compare
photoionization and heating rates with analytic expressions. Repeat with two
energies inside one diagnostic bin; the continuous result must equal the sum
of the two exact contributions and must not equal a bin-center approximation
unless by coincidence.

Verify direct and scattered path estimators separately.

### 14.4 Binning independence

Run one physical setup with diagnostic grids of 8, 16, 32, 64, and 128 bins.
Using the same continuous packet energies, require solver quantities to be
identical to roundoff for a fixed packet history. Only the diagnostic `J_E`
histogram may change.

At minimum compare:

- H/He front positions;
- volume-averaged ion fractions;
- radial C II and C III profiles;
- electron temperature;
- heating and cooling totals;
- major line luminosities.

### 14.5 HII20 acceptance test

Cloudy c25 gives approximately:

| radius | Cloudy C III fraction |
|---:|---:|
| 1.5 pc | 0.110 |
| 2.0 pc | 0.051 |
| 2.5 pc | 0.023 |

The continuous-energy MoCHII solution must retain a nonzero C III tail
outside the He II/He I transition. Agreement with Cloudy need not be exact
because recombination and charge-exchange data differ, but the result must no
longer collapse to approximately `1e-15` merely because helium is neutral.

The C III profile must also be invariant under changes to the diagnostic
energy grid.

### 14.6 HII40 acceptance test

Require:

- preservation of the already similar inner C III profile;
- recovery of C III beyond the He II region;
- focused comparison over 4.2--4.5 pc;
- no degradation of H, He, O, N, Ne, or temperature agreement.

### 14.7 Energy and MPI tests

Require:

- absorbed plus escaped luminosity to close at least as well as the grouped
  implementation;
- identical continuous source energies for the same QMC sequence at 1, 2, 3,
  5, and 8 ranks;
- no omitted or double-counted rate contributions after MPI reduction;
- successful multi-node execution;
- finite results at exact thresholds and spectrum endpoints.

## 15. Performance plan

The main risk is repeated cross-section evaluation. Prevent it by computing
all energy-only quantities once per packet. Cell traversal should then perform
only density/fraction combinations and compact dot products.

Benchmark separately:

- source CDF inversion;
- packet cross-section cache construction;
- H/He-only opacity per crossed leaf;
- metal opacity per crossed leaf;
- direct rate-tally accumulation;
- MPI reduction time and memory.

Compare continuous and grouped mode using the same HII20 and HII40 photon
counts. A target overhead of no more than 20--30% is reasonable initially,
but physical correctness is the release criterion. Optimize only after the
exact-energy tests pass.

Potential optimizations, in order:

1. flatten active metal transitions;
2. skip cross-sections below threshold without calling the fit;
3. cache packet energy-only data once;
4. vectorize the leaf composition dot products;
5. partition rate arrays consistently with the existing node-local work;
6. make the diagnostic `J_E` tally optional when it is not requested.

Do not reintroduce group-center cross-sections as a performance shortcut in
continuous mode.

## 16. Safe implementation and commit sequence

The component descriptions above explain the target design; they are not the
safe landing order. In particular, exact-energy opacity must not be activated
before direct rate scoring exists. Use the following order.

### 16.1 Freeze the grouped baseline

1. Record HII20/HII40 ion fractions, temperatures, rates, lines, energy
   conservation, iteration history, runtime, memory, and MPI configuration.
2. Add the 24.38, 24.40, 24.58, and 24.60 eV C II/He I threshold-window
   tests.

Gate: no production code has changed.

### 16.2 Add mode control and packet-energy plumbing

3. Add `ion_energy_mode`, metadata, and fatal guards for continuous features
   not yet implemented. Keep `grouped` as the default.
4. Add `photon%energy_eV`, but initialize it to `ion_e(photon%inu)` in
   grouped mode. Thread the field through stellar, diffuse, transport, and
   peel call paths without using it to change physics.

Gate: grouped results must reproduce the baseline. Any differences here are
implementation errors, not intended physics changes.

### 16.3 Refactor path scoring before changing energy

5. Route analytic direct-walk and explicit scattered-walk contributions
   through one `score_ion_path` interface.
6. Add direct H/He rate and heating arrays in grouped shadow mode, using
   `E_packet = ion_e(inu)`. Compare them with the current rates reconstructed
   from `jt_ion`.
7. Add the corresponding metal-rate shadow arrays and MPI reductions, still
   at bin-center energies.

Gate: shadow and legacy rates agree to the tolerance expected from their
different floating-point accumulation order. The legacy arrays still drive
the solver.

### 16.4 Prepare exact opacity without activating it

8. Add the per-packet H/He, metal, and dust cross-section cache, initially
   evaluated at `ion_e(inu)`.
9. Add and validate the per-leaf metal-stage population cache. Refresh it
   once per gas-state update.
10. Implement dynamic packet-energy opacity in shadow mode. At every tested
    bin center, compare it with `kap_ion(inu,leaf)`.

Gate: dynamic and precomputed grouped opacities agree. No packet is yet
transported with a continuous energy.

### 16.5 Implement and test source samplers in isolation

11. Add the piecewise-linear tabulated-spectrum conditional inverse CDF.
12. Add the band-limited Planck conditional inverse CDF for both
    `launch_sequence='random'` and `launch_sequence='sobol'`.
13. Add the full-Planck Barnett--Canfield pseudorandom reference sampler.
14. Add multi-source, external-field, and separate diffuse-stream sampler
    tests.

Gate: sampled moments, quantiles, threshold fractions, the C II/He I window,
QMC scramble behavior, and MPI rank independence all pass. The samplers are
not yet connected to production transport.

Implementation result (2026-07-28): complete. `energy_sampler_mod` contains
no diagnostic-bin dependency. Its tabulated sampler analytically integrates
and inverts each linear interval, while its Planck sampler adaptively refines
the energy-luminosity density with all supplied thresholds retained as knots.
The same monotone inverse accepts either pseudorandom or Sobol uniforms. The
standalone `tests/energy_sampler` gate also checks the independent
five-uniform Barnett--Canfield reference and reproduces the identical
global-index energy array at 1, 2, 3, 5, and 8 MPI ranks. Production launch
code remains disconnected from these samplers at this gate.

### 16.6 Preserve diffuse energy

15. Store the already sampled diffuse `eph` in `photon%energy_eV`. Keep
    grouped mode using `ion_e(inu)` until the continuous cutover.

Gate: grouped diffuse runs remain unchanged; standalone tests confirm that
continuous diffuse energies survive launch.

Implementation result (2026-07-28): complete. Both diffuse launch paths call
one fail-closed energy-assignment policy after deriving `inu`: grouped mode
stores `ion_e(inu)`, while continuous mode stores the sampled continuum,
fixed-line, or two-photon `eph`. The setup guard still prevents continuous
production transport before 16.7. The standalone survival gate covers all
representative diffuse energy ranges, and the grouped plumbing gate launches
real diffuse packets whose packet-cache assertion confirms bin-center
compatibility.

### 16.7 Activate one complete H/He vertical slice

16. For a guarded H/He-only, dust-free configuration, activate all of the
    following in one commit:

    - continuous stellar and diffuse energy;
    - cross-sections at `photon%energy_eV`;
    - dynamic exact-energy H/He opacity;
    - direct exact-energy H/He ionization estimators;
    - direct exact-energy H/He heating estimators.

Initially disallow metals, dust, secondary ionization, and unsupported
external features in continuous mode.

Gate: monochromatic-cell, threshold, Strömgren, energy-conservation,
diagnostic-bin-independence, random/QMC, and MPI tests pass.

Implementation result (2026-07-28): complete for the deliberately guarded
single internal point-source configuration. Analytic Planck and linearly
interpolated tabulated spectra use the same monotone inverse-CDF interface;
tabulated inputs are clipped to the transported band and integrated
analytically. Pseudorandom and Sobol launch paths store the sampled energy,
and diffuse launchers retain their existing continuum/line energy. The packet
cache evaluates H/He cross-sections at that energy, dynamic H/He opacity
drives both AMR and DDA walks, and direct path estimators supply the solver
Gamma/heating arrays. `jt_ion` is written only as a diagnostic histogram and
`NUBINUSE=diagnostic_only` records this separation.

The production gate verifies bitwise solver/state identity between 8 and 32
diagnostic bins for Planck, threshold-structured tabulated, diffuse, and
iterated H/He cases; pseudorandom and Sobol paths; and 1/3-rank agreement.
It also reduces emitted, absorbed, and escaped ionizing luminosities across
MPI, requires closure to `5e-13`, and records `L_EMIT`, `L_ABS`, and `L_ESC`
in the rates metadata. The existing threshold and isolated source-sampler/MPI
gates also pass.
Continuous mode remains fail-fast for multiple/external/slab sources, FUV,
secondary ionization, EUV dust, He I metastable photoionization, and peel-off
imaging until their later atomic cutovers.

### 16.8 Add metals in two steps

17. Enable exact-energy metal ionization and heating estimators with
    `ion_metal_abs=.false.`. This trace-metal configuration is sufficient to
    test recovery of the HII20 C III tail without coupling metal opacity back
    into transport.

Implementation result (2026-07-28): complete. The packet cache evaluates all
active metal-transition cross-sections once at `photon%energy_eV`. The common
path scorer accumulates exact `sigma/E` ionization and excess-energy heating
coefficients, MPI reduction preserves the existing ownership model, and
`species_gamma_compute` feeds those reduced arrays to the metal cascade in
continuous mode. Rates files expose per-element `Gamma_*_stages` and
`Heat_*_stages` arrays in addition to stage fractions. Active metal thresholds
are also mandatory Planck-CDF quadrature knots. The production gate verifies
that all three outputs are finite and bitwise identical for 8 and 32
diagnostic bins. `ion_metal_abs=.true.`, `metal_ne=.true.`, and
`metal_heat=.true.` remain fail-fast so this step does not prematurely
activate opacity or the 16.9 thermal coupling.

18. Benchmark metal scoring strategies before enabling full metal opacity:

    - all active transitions;
    - threshold-sorted early exit;
    - compact element/stage vector loops;
    - only if necessary, an unbiased stratified metal-transition estimator.

Implementation result (2026-07-28): complete. The existing flattened packet
cache is already the compact vector loop. Sorting that registry by threshold
permits an exact early exit while retaining the same `(element, stage)`
mapping. The reproducible operation-count benchmark
`tests/continuous_energy/benchmark_metal_scoring.py` integrates the
band-limited Planck packet-energy PDF. For the 12 HII C/N/O/Ne/S transitions,
the mean number of cross-section evaluations falls from 12 to 3.742 at 20 kK
and 4.668 at 40 kK, reductions of 68.8% and 61.1%. With all 32 registry
transitions, it falls to 11.420 and 13.990, reductions of 64.3% and 56.3%.
Because at most 32 deterministic evaluations remain and the exact compact
loop is already inexpensive, stochastic transition subsampling is rejected:
it would add variance and convergence diagnostics without a justified cost
benefit. Grouped opacity shadows and continuous metal `Gamma`, heating, and
stage-fraction bin-independence remain unchanged after the optimization.

19. Enable exact metal opacity using the per-leaf stage cache and the
    packet-cross-section cache.

Implementation result (2026-07-28): complete. The existing dynamic opacity
path now serves continuous production with `ion_metal_abs=.true.`. It combines
packet-cached exact-energy cross sections with the stage fractions refreshed
by `gas_opacity_fill` after each gas update; no cascade solve occurs during a
packet-cell crossing. The final consistency pass follows the same refresh and
transport sequence. Three-iteration 8/32-bin tests give bitwise-identical
H/He states, metal `Gamma`/heating, and metal stage fractions, while 1/3-rank
outputs agree to `5e-14`. A 100,000-packet HII20 AMR smoke including the
thermal solver and final consistency pass completes successfully.

That larger smoke exposed ordinary double-precision summation drift in the
independent emitted/absorbed/escaped diagnostic totals. Rank-local energy
accounting now uses at least 18 decimal digits before conversion for the MPI
reduction, retaining the `5e-13` closure gate; the HII20 smoke closure is
`1.46e-16`.

Gate: the selected estimator meets accuracy requirements at acceptable cost,
HII20/HII40 metals pass, and no metal cascade is solved inside a packet-cell
walk.

### 16.9 Restore full thermal coupling

20. Connect exact-energy metal photoheating, metal electron contribution, and
    the thermal iteration.

Implementation result (2026-07-28): complete. Both existing thermal consumers
already read `egam`/`eheat`, which are authoritative direct exact-energy
estimators in continuous mode. The continuous setup guard now permits
`metal_ne` and `metal_heat`. A three-iteration full-coupling gate enables
exact metal opacity, metal electrons, metal photoheating, and the thermal
solver together; H/He state, electron density, temperature, metal
`Gamma`/heating, and stage fractions are bitwise identical with 8 and 32
diagnostic bins. One- and three-rank outputs agree to `5e-14`.

21. Run fully converged HII20 and HII40 calculations, including the final
    consistency pass.

Gate: the C III tail is recovered, other species and temperatures do not
regress, and diagnostic energy-bin changes do not alter the physical solution
for a fixed packet history.

Implementation result (2026-07-28): complete. Both runs used
`launch_sequence='sobol'`, `qmc_seed=20260728`, `no_photons=2^25`,
`ion_energy_mode='continuous'`, 32 diagnostic bins, explicit case-A diffuse
transport including He I channels, exact metal opacity, metal electrons, and
metal photoheating.  Both ran on eight MPI ranks and required convergence
before output.  HII20 converged in 27 iterations with final energy closure
`1.996e-16`; HII40 converged in 38 with closure below printed precision.

The C III shell medians, continuous / grouped / Cloudy c25, were:

| model | radius [pc] | continuous | grouped | Cloudy c25 |
|---|---:|---:|---:|---:|
| HII20 | 1.50 | 0.0890 | `1.01e-15` | 0.1096 |
| HII20 | 2.00 | 0.0409 | `5.00e-16` | 0.0509 |
| HII20 | 2.50 | 0.0185 | `2.37e-15` | 0.0235 |
| HII40 | 4.20 | 0.257 | `1.76e-3` | 0.544 |
| HII40 | 4.35 | 0.210 | `6.88e-9` | 0.261 |
| HII40 | 4.50 | 0.156 | `1.22e-10` | 0.201 |

`tests/continuous_energy/check_hii_production.py` passes for both models.
The mean ionized-gas temperatures move only -0.54% (HII20) and -0.09%
(HII40) from the frozen grouped calculations.  All checked H/He, C, N, O,
and Ne radial-profile errors against Cloudy either improve or remain within
the predefined non-regression allowance.

### 16.10 Add secondary and dust physics one feature at a time

22. Add exact-energy secondary ionization.
23. Add FUV and `G0` estimators.
24. Add grain photoelectric heating.
25. Add exact-energy dust absorption.
26. Add elastic dust scattering while preserving `photon%energy_eV`.
27. Add exact-energy peel-off imaging support.

Each commit must remove the corresponding continuous-mode fatal guard only
after its analytic, energy, and MPI tests pass.

### 16.11 Separate diagnostics and perform the final cutover

28. Disconnect continuous-mode solver rates from binned `J_E`. Keep `J_E`
    diagnostic-only.
29. Run diagnostic grids of 8, 16, 32, 64, and 128 bins with a fixed packet
    history and require the physical state to be invariant.

Implementation result (2026-07-29): complete for both HII20 and HII40.  Each
resolution used the same `2^25` Sobol history and ran on 68 MPI ranks.  The
metadata gate and bitwise comparison passed for all physical datasets,
including H/He/metal rates, fractions, electron density, temperature, and
stage populations.  The only intentionally resolution-dependent datasets
were the diagnostic `J_nu` histogram and its `E_bin`/`dE_bin` grid.
30. Complete random/QMC ensemble, multiple-scramble, MPI 1/2/3/5/8 rank,
    multi-node, energy-conservation, performance, and memory tests.

Implementation result (2026-07-29, RQMC portion): complete.  Seeds 101, 202,
303, and 404 were run at 68 MPI ranks with the same `2^25` packet count.
All eight outputs converged, completed their final consistency pass, and
passed energy closure.  The sample relative scatter was 0.002% in HII20
ionized-gas temperature and below 0.001% in HII40; C III shell scatter was
0.261%, 0.619%, 0.630% at 1.5, 2.0, 2.5 pc in HII20 and 0.088%, 0.074%,
0.100% at 4.2, 4.35, 4.5 pc in HII40.

Implementation result (2026-07-29, node-local MPI/performance portion):
complete.  HII20/HII40 state and rate datasets at 1, 2, 3, 5, 8, and 68
ranks agree with the 68-rank reference to `rtol=5e-12`, `atol=5e-14`; all
closures pass.  Wall times [s] for 1/2/3/5/8/68 ranks were
171.45/110.46/90.17/74.86/63.36/58.09 (HII20) and
162.09/104.79/87.31/72.22/64.17/58.73 (HII40).  The measured maximum resident
set size was 0.55--0.76 GiB.  This matrix used one 72-thread node; a
multi-node allocation remains required to complete that specific release
sub-gate.

Implementation result (2026-07-29, multi-node MPI portion): complete.  The
same HII20/HII40 smoke ran across `lart2,lart3,lart4` with 68 ranks per node
(204 ranks total), retaining four threads per node.  Both outputs pass energy
closure and agree with the node-local 68-rank reference to `rtol=5e-12`,
`atol=5e-14`.  Wall times were 177.59 s (HII20) and 185.05 s (HII40), with
maximum resident sets 0.70 GiB and 0.69 GiB.  The small 100,000-packet smoke
is launch/communication dominated, so these timings validate correctness and
resource behavior rather than claiming strong scaling at production size.
31. Make continuous mode the default.
32. Update the paper manuscript itself (`MoCHII_paper/mochii.tex`) from the accepted
    continuous HII20/HII40 production outputs, then update the user
    documentation.  Regenerate every HII20/HII40 numerical value in every
    Table included in the paper, and revise every affected sentence,
    paragraph, stated deviation, conclusion, and cross-reference in the
    paper body.  Also regenerate Figures 14--15 and their captions, and
    revise the surrounding interpretation of the C III tail.  Record the
    photon count, Sobol scramble/ensemble policy, continuous-energy sampling,
    exact-energy cross sections, and diagnostic-only role of `nnu_ion`; do
    not leave grouped-output values or bin-center language in any updated
    table, caption, or discussion.  Preserve the frozen grouped results only
    where they are explicitly labelled as historical comparison baselines.
    Immediately before the single push after this stage, update the
    `Last updated:` date and time in `README.md` and include that edit in the
    final commit.

Implementation result (2026-07-29): complete.  The generated outputs reside
under `results/continuous_energy/production/` (paper inputs),
`bin_independence/`, `rqmc_ensemble/`, `mpi_scaling/`, and `mpi_multinode/`.
The production checker and both paper generators read the production paths.

Implementation result (2026-07-29, steps 31 and 35): complete.  The default
`ion_energy_mode` is `continuous`; every pre-existing MoCHII regression and
example input declares `ion_energy_mode = 'grouped'` explicitly to preserve
its established configuration.  In continuous mode `nnu_ion` is documented
and recorded as diagnostic spectral-tally resolution, and
`ion_align_edges` affects only that diagnostic binning.  Both retain their
transport/solver-bin meanings in explicit grouped mode.

### 16.12 Cleanup only after release validation

33. Decide whether to retain grouped mode for historical reproduction.
34. Only then consider removing `kap_ion`, the grouped solver-rate path, or
    other compatibility structures.
35. Reinterpret or rename `nnu_ion` as diagnostic output resolution and
    reduce `ion_align_edges` to a diagnostic-binning concern.

Do not delete the legacy path, remove `kap_ion`, change the default, or
regenerate the paper before the complete continuous mode passes its release
gates.

## 17. Completion criteria

The change is complete only when all of the following hold:

1. every stellar and diffuse ionizing packet carries a physical sampled
   energy;
2. all gas and dust cross-sections are evaluated at that energy;
3. transport opacity is evaluated at that energy;
4. solver ionization and heating rates come from exact-energy path
   estimators;
5. diagnostic spectral binning cannot change the physical solution for a
   fixed packet history;
6. the C II/He I window tests pass;
7. the HII20 and HII40 C III discrepancies caused by grouped thresholds are
   removed;
8. energy conservation, MPI behavior, and other benchmark species do not
   regress;
9. output metadata unambiguously distinguishes grouped historical results
   from continuous-energy results.

## 18. Out of scope

This plan does not require:

- changing the atomic cross-section fits themselves;
- changing carbon recombination or charge-exchange data;
- changing the ionization-equilibrium product-chain solver;
- changing spatial-grid resolution;
- increasing the number of diagnostic energy bins as a substitute for the
  continuous treatment.

Atomic-data comparisons with Cloudy should follow the continuous-energy
change. Before this correction, grouped-energy errors and atomic-data
differences cannot be cleanly separated.

## 19. Finalized plan: unify both bands under continuous spectral sampling (2026-08-01)

**Principle.** Monte Carlo RT -- dust and photoionization alike -- samples each
packet's energy/wavelength continuously from the source spectrum and evaluates
every cross section, opacity, albedo, phase function, and rate at that exact
value. Spectral bins exist only to **accumulate** the radiation field and the
absorbed-energy spectrum; the `J_nu`/`jt_ion` diagnostics and SEDust correctly
consume that binned accumulation (it is not a defect -- `sedust_compute_write`
operates on a binned absorbed-energy spectrum, not on photons). The defect is
bin-center *transport*: quantizing the transported packet's energy to a bin
center (`ion_e(inu)`, `kap_ion(inu)`) and reading cross sections by bin index
(`ion_dust_*(inu)`, `sed_sext(il)`). "grouped" is that defect, not a mode to
keep. The goal is one continuous sampler spanning both bands, bins left
diagnostic-only.

**Reference: MoCafe v2.00 `1473f44` "continuous photon wavelength draw"** did
exactly this for the dust/wavelength side and is the pattern MoCHII mirrors
(MoCHII already implements it for the ionizing/energy side via
`energy_sampler_mod`). Its five contracts, verified against the MoCafe source:

1. **Continuous sampler** (`src/spectrum_sampler_mod.f90`):
   `spectrum_sampler_type{lam, dens, cdf, total}`, built by
   `tabulated_spectrum_sampler`/`planck_spectrum_sampler` (both forced to carry
   the transfer bin edges as CDF nodes), sampled by `sample_wavelength(sampler,
   u)` -- the **exact analytic inverse** of the piecewise-linear-density CDF (a
   stable quadratic root per interval, not a power-law approximation).
   `band_luminosity_fraction` is exact when both ends are nodes. Source and
   external field get separate sampler instances; multiple sources get
   `src_spectrum(:)`.
2. **Exact-value evaluators** (`*_at`): `sed_cext_at` (log-log),
   `sed_sext_at = cext_at/cext_ref`, `sed_albedo_at`/`sed_hgg_at` (linear in
   `ln lambda`), all interpolating the retained cross-section table at the
   photon's own wavelength. The bin-center arrays survive only for output writers.
3. **Diagnostic bin address**: `photon%lambda` is continuous; `photon%il =
   sed_bin_of(lambda)` (closed-form on the log grid) addresses only the output
   image planes and the `jt_sum` tally.
4. **Continuous transport, discrete accumulation**: the packet transports at its
   exact wavelength (`s_ext/albedo/hgg = *_at(lambda)`); `jt_sum(il,cell)` is the
   binned radiation-field histogram, while a separate **`jt_abs(cell)`**
   accumulates the wavelength-integrated absorbed power with each photon's own
   `s_ext*(1-albedo)` -- never a bin-average cross section, which would bias the
   dust luminosity and temperature. Dust emission then normalizes to
   `rhokap*jt_abs` and takes its shape from `jt_sum`, drawing the re-emission
   wavelength inside a bin with `sample_power_law_bin` (reusing the
   bin-selecting uniform's remainder, no new QMC dimension).
5. **Fixed-dimension QMC contract**: the wavelength keeps its uniform slot
   (stellar `uq(3)`, dust-emission `ud(5)`); the `_u(uf)` sampler variants invert
   the monotone CDF on that one slot, so Sobol stratification and rank-invariance
   hold with no added dimension.

MoCafe ships two self-checking gates (`tests/test_spectrum_sampler.f90`,
`tests/test_power_law_bin.f90`) MoCHII will mirror at the foundation stage (S0
in Section 19.2): stratified quantile count within one draw, `cdf(sample(u)) == u`,
blackbody band fraction vs precise integration, bin-index bracketing, the size
of the bin-center bias, and the power-law-bin integral against `logarithmic_mean`.

This supersedes the earlier "parity then delete" framing: grouped is a defect to
remove, not a mode to match. Section 16.12 items 33--34 are hereby closed as
*remove*. The stages below give each feature its exact-energy path first only
because removing grouped before that would drop working capability, not because
grouped is a baseline worth preserving.

**Execution progress (2026-08-01).** Done and verified (all changes additive /
minimal, grouped physics bit-unchanged, uncommitted):
- **S0** dust cross sections at exact energy + `dust_abs_shadow` accumulator +
  continuous `heat_dust` wiring (`gas_opacity_mod`, `ion_score_mod`, `main`,
  `ion_packet_mod`, `gas_rates_mod`). Gate: grouped dust-heat shadow matches the
  binned form to 2--6e-16; a continuous dusty run closes energy to 0.
- **P3 (dust absorption)** and **P4 (dust scattering)**: the `ion_add_dust`
  guard is removed (`setup.f90`); continuous dust scattering reproduces grouped
  (`d_dust_scat`: 6.23e6 vs 6.19e6 scattered segments, closure 0).
- **P6 (He I metastable)**: guard removed; the 2^3S/10830 diagnostic runs in
  continuous on the binned `jt_ion` spectrum (soft-UV band floor same as
  grouped, i.e. needs P5's `add_fuv`).
- **P7 (peel-off)**: guard removed; direct and dust-scattered peel both produce
  images in continuous (`sum(scatt) > 0` with dust present).
- **P2 (secondary ionization)**: exact-energy nine-term path estimator
  (`sec_shadow` in `ion_score_mod`) + `ion_score_apply_sec`; `secion_apply`
  reused unchanged; guard removed. Gate: the estimator reproduces the grouped
  binned nine terms to 2--8e-16; continuous ON vs OFF shifts `<T[NpNe]>` by
  -1.25 K (cooler with secondary ionization, the correct direction).
- **P5 (FUV band + G0 + grain photoelectric heating)**: the continuous source
  sampler now extends down to `par%efuv_min` when `add_fuv` (Planck knots gain
  `eion_min`; tabulated lower bound lowered), and `ion_Ltot` is renormalized to
  `par%luminosity / f_ion` with `f_ion = 1 - CDF(eion_min)` so the ionizing
  segment carries `par%luminosity` and the FUV part rides on the spectrum shape,
  exactly matching the grouped `sum(ion_lum)`. Source FUV photons take their true
  FUV diagnostic bin (`ion_bin_of(.., include_fuv=.true.)`; the diffuse-recomb
  clamp is preserved via the optional argument). G0 is a transport-path estimator
  (`fuv_field_shadow = Sum(path_lum)` over segments with `energy_eV < eion_min`,
  `ion_score_apply_g0` scaling it by `fac/1.6e-3`), and grain PE consumes the
  resulting `g0_fuv`; `grain_pe` gains a physical-consistency guard (needs
  `add_fuv` + `ion_add_dust`) and the `add_fuv` continuous guard is removed. A
  stale pre-P5 override (`ion_Ltot = par%luminosity` inside the single-component
  fast path) was found by the A/B and removed: `setup_continuous_source_sampler`
  is now the single authority for the continuous `ion_Ltot`, so `log_fast_path`
  no longer misreports the band total in continuous. Gates: grouped G0 path
  shadow PASS at ~1e-15 (alongside dust/sec/H-He/metal, all ~1e-15); a continuous
  `pdr_full` run launches the FUV band (`32 ionizing + 16 FUV bins`), closes
  energy to 0/1e-16, reports `band total with FUV = 6.664e38, L_FUV/L_ion = 1.097`
  identical to grouped, and tracks the grouped iteration-by-iteration convergence
  (the ~7.6 cell-max `|dTe|/Te` at this low photon count is a property of the PDR
  front shared by both modes, not a P5 regression; the `vol` criterion that
  actually gates reads ~6-10%).

- **P1 (external / multiple / slab source generality)**: the scalar
  `source_energy_sampler` is promoted to a per-component array indexed exactly
  like `comp_kind`/`src_cdf`/`ion_cdf_src` (internal points 1..nsource, then the
  external field; a single-component run is `ncomp = 1`). Each component builds
  its own sampler from its own spectrum: `build_point_source_sampler` follows the
  `point_shape`/`comp_shape` chain (`src_spectrum_file` column `is` >
  `src_tstar(is)` Planck, multi only > `ion_spectrum` > `par%tstar`) and
  `build_external_field_sampler` follows the `build_ext_bins` chain (ISRF preset >
  `ext_spectrum` file > `ext_tstar` Planck > `ion_spectrum` > `par%tstar`). ISRF
  presets are made continuous by sampling the analytic `preset_je` onto a dense
  4096-point grid and inverting it like any tabulated spectrum. The six launch
  sites (`gen_ion_photon` single/multi, `gen_ion_photon_qmc` point/slab/external/
  multi) draw the frequency from the component sampler on the **same QMC slot**
  as grouped, so `ion_qmc_ndim` and rank are unchanged. The continuous `ion_Ltot`
  override in `setup_continuous_source_sampler` is dropped so both modes use
  `sum(ion_lum)` (the pre-P5 single-source override was already gone; the
  external-only/slab/multi cases it mis-set are now correct), and a mode-agnostic
  `if (par%luminosity <= 0) par%luminosity = ion_Ltot` fills the derived-
  luminosity report. The `setup.f90` continuous guard is reduced to "no source
  spectrum at all". Gates: the external decoy test closes the silent-wrong-
  spectrum trap (`ext_tstar=4e4` and the `ext_spectrum` 4e4 file both give
  absorbed 7.06e31, agreeing to 1.5e-5, while stripping `ext_tstar` to expose the
  2e4 decoy gives 1.03e32 -- 46% off); component-spectrum separation holds (two
  sources at 3e4/5e4 K sum to the mixed run to 0.03%); the grouped shadow gates
  (dynamic-opacity, H/He, metal, dust-heat, G0) all still PASS at ~1e-15 after the
  shared-code change; and the single-point continuous run is unchanged (band total
  6.664e38, emitted moved +9e-7 from the override removal only).
- **P1b (sampler-based `ion_Ltot`, the former open item §6a)**: the continuous
  band total and the component-selection CDF no longer come from the binned
  `sum(ion_lum)`/`comp_L` (NSUB=32 midpoint quadrature) but from the samplers
  themselves, so the bin grid is finally diagnostic-only on the normalization too.
  `setup_continuous_source_sampler` collects a `band_L` per component and sets
  `ion_Ltot = sum(band_L)`, rebuilding `src_cdf` from the same numbers; the binned
  `ion_lum`/`ion_cdf`/`ion_cdf_src` are untouched (grouped keeps them). Each
  component's rule reproduces what grouped holds fixed, with the FUV/ionizing
  split taken from `ionizing_fraction = 1 - CDF(eion_min)` instead of the
  quadrature: a shape spectrum gives `scale/f_ion` (`par%luminosity`, or
  `src_lum(is)` when multi), an absolute file with no rescale gives the sampler's
  `total_integral`, an external absolute field gives `pi*A*ext_intensity` when
  rescaled (analytically exact) or `pi*A*total_integral` otherwise, a legacy
  external shape gives `pi*ext_intensity*A/f_ion`, and slab gives
  `sum(slab_Fz)*A_zface/f_ion`; a FUV-only spectrum asked to carry an ionizing
  luminosity aborts exactly as grouped's `finalize_source` does. Gates: the
  continuous `pdr_full` emitted returns to **6.663991e38** (the sampler value; it
  read the binned 6.663997e38 after P1) with closure 0, and the run now prints
  both the binned diagnostic and the authoritative `continuous band total
  (sampler)`; grouped is bit-unchanged (the five shadow gates report identical
  values, band total 6.6640e38 / `L_FUV/L_ion` 1.097); and the external analytic
  rule checks out -- `ext_intensity = 2.0e-5` gives band total 5.743177e34 against
  2.871588e34 at 1.0e-5, a ratio of exactly 2.000000. Every component rule was
  then exercised against the luminosities grouped prints for the same input, with
  the continuous `emitted` equal to the sampler band total and the closure at
  roundoff in all five: `two_point` 4.000000e37 against `L_ion` 3.0e37 + 1.0e37
  (exact, `f_ion = 1` without `add_fuv`), `mixed` 1.002872e37 against 1.0e37 +
  `L_enter` 2.8716e34, `st_draine` 6.124264e35 against the preset's `pi*A*J`
  (the 4096-point sampler integral replacing the NSUB=32 binned one), `slab`
  7.302827e28 against `L_slab` 7.3028e28, and `st_extje2` as above.

Remaining: **D** (re-baseline + grouped removal + the `.md`/`.tex` sweep) --
**not to be started yet** (user hold).

### 19.1 Current MoCHII state

The ionizing band already follows the reference pattern: `energy_sampler_mod`
draws a continuous energy, `ion_bin_of` labels the diagnostic bin, and the H/He
and metal rates come from exact-energy path estimators. What remains is (a) the
dust/SED band, still on the old bin pattern MoCHII inherited from MoCafe
(`sample_sed_lambda -> il`, `sed_sext(il)`, and the bin-indexed
`ion_dust_*(inu)` used during ionization) with no `jt_abs`-style exact-wavelength
absorbed accumulator; and (b) the grouped bin-center transport mode plus the
`continuous_reason` guard that pins six feature families to it. `setup.f90`'s
guard refuses continuous, with `error stop 1`, whenever the input uses any of

| guarded feature (`par%`)            | grouped-pinned inputs that use it |
|-------------------------------------|-----------------------------------|
| external / `nsource>1` / `slab`     | 26 (external+slab), 14 (`nsource>1`) |
| `ion_add_dust`                      | 23 |
| `add_fuv`                           | 21 |
| `ion_peel`                          | 11 |
| `hei_metastable`                    | 6  |
| `use_sec_ion`                       | 1  |

Of the 176 inputs that pin `ion_energy_mode = 'grouped'` (159 under `tests/`,
17 under `examples/`), 77 use at least one feature continuous refuses; the
other 99 use none but carry grouped bin-center baselines that continuous would
change. So deleting grouped today would delete the only working transport path
for external/multiple sources, slab geometry, dust-coupled ionization, FUV/G0,
grain photoelectric heating, secondary ionization, He I metastable, and
peel-off imaging. Each of these must be given its exact-energy path before the
grouped bin-center transport can be removed -- not to preserve grouped, but to
keep the capability while the transport is corrected.

### 19.2 Work stages (each gives one feature its exact-energy path, then drops its guard)

A four-cluster read-only source audit (2026-08-01) established that **every
guarded feature can be made continuous with no fundamental blocker.** Continuous
sampling is not tied to the ionizing (LyC) band: `energy_sampler_mod` inverts a
CDF over whatever energy range the sampler is built on, so an FUV-only source is
sampled continuously exactly like an ionizing one -- its CDF simply lives in the
FUV band. The ISRF presets (`draine`/`habing`/`mathis`, which
`ion_band_mod.f90`'s `draine_je`/`mathis_jl_si` place below 0.0912 um and which
force `add_fuv`) are therefore not a special case: they are FUV sources that
fold into P5 (FUV-band sampling + transport), like any other `add_fuv` input.
The only reason they cannot run continuously today is the same one that blocks
every FUV source -- the sampler currently floors at `par%eion_min` -- which P5
removes. The audit's structural finding is why the rest is tractable: continuous rate
estimators read only `photon%ionphys%energy_eV` and the cached cross sections
(`ion_score_mod.f90:131-157`, `species_mod.f90:343-357`) with no geometry or
source dependence; all transport opacity routes through `ion_packet_opacity`,
which returns the energy-based dynamic opacity in continuous (AMR and car/DDA
alike); and the Tier-2 diagnostics (`nebcont_mod`, `lines_mod`, `sh95_mod`,
`nlevel_mod`, `recomb_mod`) reference neither `jt_ion` nor `ion_e` -- they read
the converged gas state only. So most stages are launch-side or one added path
estimator, not transport surgery.

Stages are grouped by the work each needs. Each stage's gate is the existing
shadow harness (`par%ion_shadow_rates`, continuous vs grouped to rtol) plus
energy closure and MPI/diagnostic-bin invariance; where continuous is genuinely
more accurate the grouped number is retired rather than matched (below). No
stage removes its guard until its analytic, energy, and MPI tests pass, and P2
and P5 must land their estimator **before** the guard lifts (their `apply`
routines are silent no-ops in continuous until then).

**Foundation (S0) -- exact-energy dust cross sections and the split accumulator.**
MoCHII does not need MoCafe's new sampler module: `energy_sampler_mod`'s
`sample_energy_cdf` is already the same exact analytic inverse-CDF as MoCafe's
`sample_wavelength` (binary search + the stable quadratic root
`2*area/(y0+disc)` with a small-slope linear branch), and the continuous
ionizing photon already carries `energy_eV`. What the dust band lacks is the
`*_at`/`jt_abs` half of the MoCafe pattern, in energy space:
(i) retain the dust `kext` table at module scope (`gas_opacity_mod` currently
builds `ion_dust_*` into bins from a local table and discards it) and add
energy-space evaluators `ion_dust_{sabs,ssca,g}_at(E)` (energy -> wavelength
`1.23984/E` -> interpolate, exactly as the bin arrays are built) replacing the
bin-center `ion_dust_*(inu)` lookups;
(ii) add a `jt_abs`-style exact-energy absorbed-power accumulator so `heat_dust`
and the dust temperature/IR come from each photon's own cross section, not a bin
average;
(iii) wire `ion_packet_opacity` and the dust-heating estimator to the
`*_at(energy_eV)` values.
MoCHII does **not** transport dust re-emission -- SEDust is an output-time SED
built from `heat_dust` plus the binned `jt_ion` (`main.f90:334-337`) -- so
MoCafe's `sample_power_law_bin` contract does not apply here, and there is no new
dust-wavelength sampler or `sed_bin_of` to add (the ionizing `ion_bin_of`
already labels the diagnostic bin). Gate S0: `ion_dust_*_at(ion_e(b))` reproduces
the existing `ion_dust_*(b)` at every bin center to round-off (same table), and
the `jt_abs` accumulator closes against the absorbed luminosity. Smaller than
MoCafe's because MoCHII already has the sampler and does not re-emit dust
photons; the MoCafe self-check gates (`test_spectrum_sampler`) still port over as
a standing check on `energy_sampler_mod` itself.

**Nearly free -- a guard flip plus verification (no new physics):**

- **P6 -- He I metastable.** The guard is conservative, not a missing feature:
  `jt_ion` is populated in continuous (`ion_score_mod.f90:123`, before the
  continuous early-return), and `hei_metastable_run` is a feedback-free output
  diagnostic that recomputes `Phi_3` from `jt_ion` and the converged gas state
  (`hei_metastable_mod.f90:110-144`). Minimal path: remove the `hei_metastable`
  guard and run on the diagnostic spectrum (bin-center `Phi_3`, acceptable for a
  no-feedback diagnostic). Optional exact-energy parity: cache `sigma_3` and add
  a path estimator like H/He. Soft dependency: the soft-UV part below the band
  floor (`4.78 eV`) is only captured with `add_fuv` (P5) -- the same limitation
  grouped already has, so not a blocker. Gate: `tests/hei_meta`,
  `examples/hei_metastable`.
- **P7 (direct) -- peel-off imaging.** Direct peel is physically ready:
  brightness is `Lpacket*wgt*exp(-tau)` and the LOS `tau` uses the energy-based
  `ion_packet_opacity` after `build_ion_packet_physics(pobs)`
  (`ion_peel_mod.f90:261-271`, `raytrace_amr.f90:760`); `ion_e(inu)` enters only
  the diagnostic channel/cube binning. Only the `setup.f90:504` guard blocks it;
  no new routines. Dust-scattered peel is not part of P7 -- it rides on P4
  (`main.f90:255-256` gates the scatter-peel hook on `ion_add_dust`). Gate:
  `tests/peel`, `tests/peel_ext`.

**Self-contained -- one added exact-energy path estimator each:**

- **P2 -- Secondary ionization** (item 22). The `si_*`/`heat_hard_*` band
  integrals live only in the grouped binned loop (`gas_rates_mod.f90:159-184`);
  continuous returns early (`:117-128`) and `secion_apply` then reads zeros --
  a silent no-op. Add an exact-energy estimator cloning the H/He skeleton
  (`ion_score_mod.f90:131-142`): accumulate the `E0 = h nu - E_th > 40 eV`
  hard-heating and secondary-ionization terms from `photon%ionphys%energy_eV`,
  reduce, apply; reuse `secion_apply`'s x-dependent part unchanged
  (`gas_rates_mod.f90:214-265`). Exact energy is more accurate (per-photon
  threshold, not bin-center). No dependence on other P-stages. Remove the
  `use_sec_ion` guard. Gate: `tests/g2_hii/sec_on`, plus a physical check where
  the two legitimately differ near threshold.
- **P3 -- Dust absorption heating** (items 24--25, absorption part). Extinction
  is already correct (dust is in the transport opacity), but the absorbed energy
  is scored nowhere in continuous: `heat_dust` has its only home in the binned
  loop (`gas_rates_mod.f90:186-192`) and stays zero, so dust temperature, IR,
  and SED are all zero in continuous. Add a leaf estimator
  `dust_absorbed_path(il) += path_lum*photon%ionphys%dust_abs` in
  `score_ion_path` + reduce + an apply that reproduces the grouped expression.
  Also interpolate the dust cross sections at `energy_eV`: the `ionphys` dust
  cache is still bin-indexed (`ion_packet_mod.f90:53-56`, `ion_dust_*(inu)`).
  `dust_temp`/IR then work unchanged (they read only `heat_dust`). Remove the
  `ion_add_dust` guard. Gate: `tests/d_dusty` absorption cases; dust heating
  closes.

**Rides on P3 -- opens with the same `ion_add_dust` guard:**

- **P4 -- Elastic dust scattering** (item 26). Already implemented and
  energy-mode-independent: the scattering loop (`raytrace_amr.f90:248-322`)
  changes direction and weight only, and `photon%ionphys`/`energy_eV` are built
  once before the loop (`:258`) and preserved (elastic). No recalculation
  needed; it becomes reachable when P3 opens `ion_add_dust`. Scatter-peel
  imaging then also becomes reachable (needs P7's `ion_peel` too). Gate:
  `tests/d_dusty/d_dust_scat`, `tests/peel/peel_scat`.

**Largest -- new launch-side infrastructure:**

- **P1 -- Source generality.** Estimators and transport are already
  geometry-independent, so this is purely launch-side. Promote the single scalar
  `source_energy_sampler` (`ion_band_mod.f90:51`) to a per-component array and
  build each from its own spectrum (external `ext_spectrum`/`ext_tstar`,
  multiple-source columns via `read_*_table(..., is, nsource, ...)`, Planck with
  the metal-threshold knots). Wire the continuous branches that are missing: the
  non-QMC multiple-source path (`ion_band_mod.f90:1518-1544`, no continuous
  branch today), the non-QMC external sampler selection (`:1490-1517` currently
  feeds every geometry the point-only sampler -- a silent-wrong-spectrum trap),
  and the three QMC branches (slab `:1603`, external `:1613`, multi `:1632`),
  keeping the frequency uniform on the **same slot** (`ion_qmc_ndim` does not
  depend on the energy mode, so rank stays invariant). **Fix the `ion_Ltot`
  override** (`ion_band_mod.f90:206-208` and `:268`): it sets `ion_Ltot =
  par%luminosity` inside the single-component path, which is correct for a point
  source but overwrites the true `sum(ion_lum)` for slab and external-only.
  Relax the guards (point-only, no-external, internal-spectrum-required) into a
  general internal/external sampler build (FUV-only ISRF presets need no special
  handling here -- they arrive with P5's FUV-band sampler). Gate:
  `tests/ext_field`, `tests/multi_src`, `tests/slab`, `tests/peel_ext` in
  continuous, matched to grouped to the shadow rtol as `nnu_ion` refines, with
  rank-invariance at 1/2/3/5/8.
- **P5 -- FUV band + G0 + grain photoelectric heating** (item 23). The real
  prerequisite is coverage, not just an estimator: the continuous source sampler
  floors at `par%eion_min` (`ion_band_mod.f90:264-273`), so no FUV photon is
  ever launched. Extend the sampler down to `par%efuv_min` when `add_fuv`, and
  add the FUV luminosity to `ion_Ltot` as the grouped path does. Then add a
  `G0` path estimator (`fuv_field_path(il) += path_lum` for
  `energy_eV < eion_min`) reproducing the grouped `g0_fuv`
  (`gas_rates_mod.f90:193-195`). **Grain PE belongs here, not with P3**: it
  consumes `g0_fuv`, not `heat_dust` (`thermal_mod.f90:93,102-104`), so it needs
  the FUV field of P5 with dust present (P3's `ion_add_dust` open for `rhokap`).
  `grain_pe` currently has no setup guard at all and is a silent no-op in
  continuous -- add a guard consistent with P5. Remove the `add_fuv` guard.
  Because the sampler now covers the FUV band, the FUV-only ISRF presets
  (`draine`/`habing`/`mathis`) also become continuous here -- they are FUV
  sources with an FUV-band CDF, no longer a special case. Gate:
  `tests/d_dusty/*fuv*`, `tests/pdr`, `tests/g0_gamma`, and an ISRF-preset case
  (`examples/isrf_cloud`).

Dependency summary: P6 and P7-direct are guard flips; P1 and P2 are independent;
P3 -> P4 (and P3 + P7 -> scatter-peel); P3 + P5 -> grain PE. Where a stage
should converge to grouped (bin-center -> exact as `nnu_ion` grows), the gate is
the shadow rtol; where continuous is genuinely more accurate (per-photon
secondary-ionization threshold, the C III window, exact-energy dust/FUV), the
grouped number is retired rather than matched.

**Before D, add regression gates for three paths that are continuous-safe by
inspection but never exercised in continuous:** solution-driven re-refinement
(`refine_front`, off in every continuous input), the car/DDA backend (continuous
tests use AMR only, but car raytrace routes through `ion_packet_opacity`), and
`xy_periodic`/plane-parallel slab. No defect was found in any; they lack
coverage, not correctness.

### 19.3 Deletion stage D (only after S0 and P1--P7 all pass)

**D1 -- re-baseline the 176 inputs.** Remove the `ion_energy_mode = 'grouped'`
line from every pinned input (they then take the continuous default), rerun,
and regenerate each recorded baseline (`.log`, and the JSON/`_new.log`
snapshots). This is the largest surface and must come last, because until every
guard is gone some of these inputs will not run in continuous. Numbers change
from bin-center to exact energy; that is the intended physics, so each gate's
anchor is re-recorded, not matched to the old grouped value.

**D2 -- delete grouped-only code.** With grouped gone the packet energy is
always the sampled energy, so:

- `src/ion_energy_policy_mod.f90` -- whole module; its callers in
  `ion_band_mod` and `diffuse_mod` set `photon%energy_eV = sampled_energy`
  directly.
- the shadow-comparison harness (its sole purpose was proving continuous ==
  grouped): `par%ion_shadow_rates`; `hhe_shadow`/`hhe_coeff` and the
  `species_shadow_*`/`metal_shadow%coeff` estimators; `packet_opacity_shadow_*`;
  the `jt_ion` reconstruction path in `ion_score_mod`.
- `ion_cdf`/`ion_cdf_src` (grouped source bin-CDF) once P1's per-component
  samplers replace them; the grouped branches of `gen_ion_photon` and
  `gen_ion_photon_qmc`; the `kap_ion(inu,leaf)` per-bin opacity table path in
  `ion_packet_mod`.
- the grouped branches in `gas_rates_mod` and `diffuse_mod`, and the
  `if (continuous)` guards in `species_mod`/`ion_score_mod` become unconditional.
- `src/setup.f90` mode-string validation + the `continuous_reason` guard block;
  `par%ion_energy_mode` and `par%ion_shadow_rates` in `src/define.f90`.

**D3 -- keep (now diagnostic-only, drop "grouped" from comments).**
`ion_e`/`ion_nu`/`eedge`/`ion_bin_of` and `photon%inu` feed the `jt_ion`/`J_nu`
diagnostic spectrum; `energy_sampler_mod` is the continuous core; `nnu_ion` and
`ion_align_edges` keep only their diagnostic-resolution meaning (Section 16.12
item 35).

**D4 -- retire grouped-specific tests.**
`tests/continuous_energy/{grouped_plumbing.in,check_grouped_plumbing.py,
check_grouped_baseline.py,grouped_baseline.json}` and the
`diffuse_energy_policy_test.f90` policy gate lose their subject.

**D5 -- final step of the plan: update every `.md` and `.tex` file.** A
repository-wide documentation sweep is the last stage. In *every* tracked
`.md` and `.tex` -- `docs/*.md`, `docs/*.tex`, `MoCHII_paper/*.tex`, `CLAUDE.md`,
`README.md`, and any other -- remove the "two modes"/`grouped`/bin-center/
`ion_energy_mode`/`ion_shadow_rates` language and state the single continuous
spectral model: both bands sampled continuously, all cross sections at the exact
energy/wavelength, bins diagnostic-only, and the `jt_abs`-style exact-wavelength
absorbed accumulator feeding dust/SEDust. Known specifics: drop the
`ion_energy_mode`/`ion_shadow_rates` rows from the `par%` reference in
`docs/MoCHII_UserGuide.tex`; revise `docs/MoCHII_physics.tex` (charge-exchange,
transport, and cooling sections that describe bin-center energies),
`docs/MoCHII_cooling_analysis.tex`,
`docs/PHOTON_ENERGY_SAMPLING_CODE_COMPARISON.tex`, `MoCHII_paper/mochii.tex`
(already continuous -- confirm no bin-center wording survives), `CLAUDE.md`,
`docs/PLAN.md`, and the plan documents that reference the old two-mode design;
mark this document complete. Preserve frozen grouped results only where they are
explicitly labelled as historical baselines. Verify with a grep that no
`.md`/`.tex` mentions `ion_energy_mode`/`grouped`/bin-center transport except as
labelled history, then rebuild every `.tex` (`docs/make_pdf.sh` and the paper
build) and inspect the layout. Refresh `README.md`'s `Last updated:` in the same
push.

### 19.4 Risks and latent traps

- **Ordering is load-bearing:** D1 must follow every P-stage, or pinned inputs
  hit a still-present guard. The parity stages are independent except
  P4/scatter-peel and grain PE (need P3), and grain PE also needs P5.
- **A guard must not be lifted before its estimator lands.** `secion_apply`
  (P2) and the grain-PE/`g0_fuv` path (P5) are silent no-ops in continuous until
  their exact-energy estimators exist; removing the guard first would produce
  quietly wrong (not crashing) results.
- **Latent luminosity bug in P1.** `ion_Ltot = par%luminosity`
  (`ion_band_mod.f90:206-208`, `:268`) is inside the single-component path and
  overwrites the true `sum(ion_lum)` for slab and external-only sources. It is
  masked by the guard today; fix it as part of P1, not after.
- **Silent-wrong-spectrum trap in P1.** The non-QMC single-component continuous
  branch (`ion_band_mod.f90:1490-1517`) already calls `sample_energy_cdf` for
  every geometry but the sampler only knows the internal point spectrum, so
  lifting the guard without wiring per-component samplers samples the wrong
  spectrum silently.
- **`grain_pe` has no setup guard** (absent from the `continuous_reason` list);
  it runs in continuous but is a no-op because `g0_fuv = 0`. Add a guard in P5.
- **QMC determinism:** P1 changes the busiest launch path; keep the frequency
  uniform on the same slot (`ion_qmc_ndim` is mode-independent) or every QMC
  baseline shifts for the wrong reason. Verify rank- and scramble-invariance
  before touching baselines.
- **Baseline churn is physical, not a bug:** ~99 inputs that already run in
  continuous will still change (bin-center -> exact). Re-anchor, do not chase
  bit-identity with the retired grouped numbers.
- **`sedust` SED shape draws from binned `jt_ion`** (`sedust_mod.f90:112-180`)
  while its normalization is exact via `heat_dust`; in continuous the shape is
  the diagnostic-bin spectrum, a minor inconsistency to document (or lift to a
  `jt_ion`-refinement) when P3 opens dust.
- **Reversibility:** grouped is deleted only at D2; keep P1--P7 as separate
  reviewed commits so the removal can be staged and, if a parity gap surfaces
  late, deferred without unwinding the parity work.

## 20. Sampling references

- C. Barnett and E. Canfield, *Sampling a Random Variable Distributed
  According to Planck's Law*, Lawrence Radiation Laboratory report (1970),
  [NTIS record](https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/DE97050749.xhtml).
- J. E. Bjorkman and K. Wood, *Radiative Equilibrium and Temperature
  Correction in Monte Carlo Radiation Transfer*, ApJ 554, 615 (2001),
  [arXiv:astro-ph/0103249](https://arxiv.org/abs/astro-ph/0103249).
