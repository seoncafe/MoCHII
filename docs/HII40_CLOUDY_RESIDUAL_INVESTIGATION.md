# HII40 C2+ / He I Front Residual Investigation

**Status:** active investigation (2026-07-29)  
**Reference:** Figure 14 in `paper/mochii.tex`; Cloudy c25 HII40 output in
`tests/cloudy_c25/hii40_c25.*`.

## Question

After continuous packet-energy transport was introduced, the HII40 C2+
profile agreed much more closely with Cloudy than the earlier grouped result.
Nevertheless, the C2+ decline and the associated He I front remained inward
of the Cloudy solution.  This note records the numerical and physical tests
performed to identify the origin of that residual.

## Common diagnostic

All reported MoCHII front positions use volume-weighted spherical shells of
width 0.01 pc, followed by linear interpolation between shell centres.  The
diagnostic fields are `x_HeI = 0.5`, `x_c_stages(:,3) = 0.5` (C2+), and
`x_HI = 0.1`.  The analysis is implemented in
`tests/continuous_energy/analyze_hii40_residuals.py`.

| Case | r(He0 = 0.5) [pc] | r(C2+ = 0.5) [pc] | r(H0 = 0.1) [pc] | C2+(4.20 pc) | C2+(4.35 pc) | C2+(4.50 pc) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Cloudy c25 | 4.30531 | 4.22215 | 4.71311 | 0.54343 | 0.26062 | 0.20087 |
| MoCHII production (L6) | 4.12276 | 4.05918 | 4.62707 | 0.25609 | 0.20980 | 0.15623 |
| Diffuse field disabled (OTS) | 4.18159 | 4.10080 | 4.47390 | 0.25048 | 0.16583 | 0.01359 |
| He I excited diffuse channel disabled | 4.02014 | 3.96893 | 4.51111 | 0.22819 | 0.17517 | 0.04835 |
| Alternate MOCREC atomic data | 4.11660 | 4.03970 | 4.62831 | 0.21084 | 0.17152 | 0.12684 |
| L7 spatial grid | 4.08411 | 4.05686 | 4.66714 | 0.25191 | 0.20884 | 0.15817 |
| Ground continua at threshold | 4.09112 | 4.02432 | 4.57811 | 0.24231 | 0.19657 | 0.13390 |
| Ground continua from Milne | 4.12339 | 4.05970 | 4.62628 | 0.25621 | 0.20976 | 0.15602 |

Thus the production He I 50% front is 0.18255 pc inward of Cloudy, and the
C2+ 50% front is 0.16297 pc inward.  The co-location of these residuals is
the main reason to prioritize the He-ionizing diffuse field over C-specific
photon sampling or cross sections.

## Effects excluded or bounded

### Diagnostic energy bins

The production state is independent of the number of diagnostic bins
(`nnu_ion = 8, 16, 32, 64, 128`): the physical HDF5 arrays are bitwise
identical.  In continuous mode, a stellar or diffuse packet carries its
sampled energy and evaluates its cross sections at that energy; `inu` is only
the radiation-field output bin.  The bin-independent gate is
`tests/continuous_energy/check_hii_bin_independence.py`.

### Source-energy sampling and Monte Carlo noise

The source spectrum is sampled from its continuous CDF, and the continuous
mode uses the sampled energy for H, He, and metal opacities and rates.  Four
RQMC seeds give standard deviations of only 3.6e-5 pc for the He I 50% front,
1.4e-5 pc for the C2+ 50% front, and 6.6e-6 pc for the H I 10% front.  Those
are negligible compared with the 0.16--0.18 pc residual.

For the 40,000 K blackbody, numerical integration gives
`Q(He I)/Q(H I) = 0.10803`.  The narrow C II--He I interval from 24.383 to
24.587 eV contains 0.00495 of Q(H I), or 4.58% of Q(He I).  Photons above
100 eV are only 4.57e-10 of Q(H I).  The high-energy cutoff and this narrow
threshold interval cannot account for the observed front shift.

### Spatial resolution

The L6 leaf full width at the front is about 0.296 pc; the original He I
residual is about 0.62 L6 cell widths.  A fully converged L7 calculation
(2,097,152 leaves, 50 iterations) does not move the C2+ front outward:
`r(C2+=0.5)` changes from 4.05918 to 4.05686 pc.  The He I front moves inward
to 4.08411 pc.  Resolution is therefore not the explanation.

### H/He recombination-rate data

Replacing the relevant atomic data with the alternate MOCREC set converged
in 38 iterations.  It moves the He I and C2+ fronts inward by 0.00616 and
0.01948 pc, respectively, and makes the C2+ tail smaller.  This rate choice
does not explain the residual in the required direction.

## Diffuse-field tests

The explicit diffuse field is important, but all-off versus on is not a
single-process comparison: turning it off changes the hydrogen front and
removes the outer C2+ tail.  In particular, OTS moves the He I front outward
by 0.05883 pc but moves the H I front inward by 0.15317 pc and reduces C2+ at
4.50 pc from 0.15623 to 0.01359.

The optional He I excited-recombination channel is a cleaner paired test.
Enabling it moves the He I front outward by 0.10262 pc and the C2+ front
outward by 0.09025 pc.  It has the required sign and is material, but it does
not close the remaining gap: the production He I front is still 0.18255 pc
inside Cloudy.

## Current ground-continuum model and active test

`src/diffuse_mod.f90` currently emits ground-state H II, He II, and He III
recombination photons with the approximation

    E = E_th + kT_e [-ln(u)],

and weights their emitted luminosity by the matching mean energy
`E_th + kT_e`.  This is a near-threshold exponential approximation, not an
exact Milne free-bound spectrum including the energy-dependent inverse
photoionization cross section and free-bound Gaunt factor.

To measure the importance of this approximation before introducing a more
complex physical sampler, `par%diffuse_energy_model` now accepts:

- `exponential` (default): the existing approximation;
- `threshold`: all three ground continua emitted at their thresholds.

The threshold mode is a sensitivity bracket, not a replacement physical
model.  Its converged HII40 calculation moves the He I front inward by
0.03164 pc and the C2+ front inward by 0.03486 pc relative to the exponential
default; both changes worsen agreement with Cloudy.  Thus the finite energy
width of the present diffuse-continuum approximation has the physically
important sign, but it supplies only a small fraction of the remaining
0.16--0.21 pc discrepancy.  The runner is
`tests/continuous_energy/run_hii40_residual_experiments.sh threshold` and its
generated files are kept under `results/continuous_energy/residual_investigation/`.

## The Milne sampler: implemented, and it does not explain the residual

`par%diffuse_energy_model = 'milne'` now samples each ground continuum from
the detailed-balance photon spectrum `E^2 sigma_pi(E) exp[-(E-E_th)/kT]`,
tabulated per temperature and inverted with the transport cross sections
(`src/milne_recomb_spectrum_mod.f90`, gated by `tests/milne`).  The converged
68-rank HII40 calculation reached the same 38 gas iterations as the production
run.  Against the exponential default:

| front | exponential | milne | shift | Cloudy gap closed |
|---|---:|---:|---:|---:|
| r(He0 = 0.5) | 4.12276 | 4.12339 | **+0.00063** | 0.35% |
| r(C2+ = 0.5) | 4.05918 | 4.05970 | **+0.00052** | 0.32% |
| r(H0 = 0.1) | 4.62707 | 4.62628 | −0.00079 | — |

The residual to Cloudy moves from 0.18255 to 0.18192 pc for He I.  **The exact
free-bound spectrum is worth about 0.3% of the discrepancy.**  This settles the
question the threshold bracket raised: spectral width in the ground continua is
not where the missing 0.16-0.18 pc comes from.

The signs are instructive and are not a single "harder or softer" statement.
Each front follows the mean energy of the continuum that drives it:

| continuum | exponential `<E>` | Milne `<E>` | effect |
|---|---:|---:|---|
| H I (13.598 eV) | 14.4601 | 14.4274 | softer, H I front moves inward |
| He I (24.587 eV) | 25.4491 | 25.4603 | harder, He I front moves outward |

The H I continuum softens because the hydrogenic cross section falls as
`E^-3`, so the `E^2 sigma` weight decreases with energy.  The He I continuum
hardens because the Verner He I fit falls more slowly near threshold, leaving
`E^2 sigma` still rising there.  An estimate made from the H I channel alone
predicts the right magnitude but the wrong sign for the He I front.

## A lead the Milne work exposed

`milne_setup` integrates the same spectrum with the Milne prefactor restored
and reports the implied `alpha_1(T)` against the `alpha_A - alpha_B` that sets
the packet count.  These are independent data — Verner cross sections through
detailed balance on one side, fitted recombination coefficients on the other —
and are reported rather than reconciled:

| continuum | 5000 K | 10000 K | 20000 K |
|---|---:|---:|---:|
| H I | 0.990 | 0.995 | 1.001 |
| He I | **1.034** | **1.028** | **1.020** |
| He II | 0.998 | 0.996 | 0.996 |

H I and He II agree to about 0.5%, which is the expected level for independent
fits.  He I disagrees by 2-3.4%.  The unexplained residual is the He I front,
so an internal inconsistency of that size in exactly the He I ground-continuum
data is worth pursuing, though it is a lead rather than a demonstrated cause:
a 3% error in the He I ground recombination rate is not obviously enough to
move a front by 0.18 pc.

## Next decision

Ground-continuum spectral shape is now excluded, along with diagnostic
binning, Monte Carlo noise, spatial resolution, and the H/He recombination
rate set.  What remains, in the order the evidence supports:

1. **The He I excited cascade.** Enabling it already moves the He I front
   0.103 pc outward, the largest single effect measured, and it is still the
   crudest piece of physics in the diffuse field: fixed AGN3 branching, and a
   flat two-photon spectrum in place of the real 2^1S shape.  A level-resolved
   cascade with a physical two-photon distribution is the obvious next step.
2. **He I ground-continuum atomic data**, given the `alpha_1` disagreement
   above.
3. **Diffuse transfer** itself, as distinct from the emitted spectrum.
