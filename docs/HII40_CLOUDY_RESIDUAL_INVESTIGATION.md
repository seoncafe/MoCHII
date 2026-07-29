# HII40 C2+ / He I Front Residual Investigation

**Status:** active investigation (2026-07-30)  
**Reference:** Figure 14 in `paper/mochii.tex`; Cloudy c25 HII40 output in
`tests/cloudy_c25/hii40_c25.*`.

## Question

After continuous packet-energy transport was introduced, the HII40 C2+
profile agreed much more closely with Cloudy than the earlier grouped result.
Nevertheless, the C2+ decline and the associated He I front remained inward
of the Cloudy solution.  This note records the numerical and physical tests
performed to identify the origin of that residual.

## Central conclusion

**The residual is not in the diffuse field --- neither in what it emits nor in
how it is transported.**

*Emission.*  Four independent corrections have been made to the emitted diffuse
spectrum: the Milne free-bound shape of the three ground continua, the Drake
two-photon shape, the temperature-dependent He I cascade branch fractions, and
the local fate of the He I 584 A photon.  Between them they span the spectral
shape of every continuum, the weights of every branch, and the destination of
the channel's strongest line, and the last of them changes that channel's
content by 9-11% rather than by a fraction of a percent.  Their combined effect
on `r(He0 = 0.5)` is **-0.00306 pc**: 1.7% of the 0.18 pc gap, and in the wrong
direction.  No remaining refinement of the emitted spectrum can plausibly
supply 0.18 pc, because the corrections already made include one large enough to
have done so and it did not.

*Transport.*  This was the leading candidate when the note was last revised, on
the grounds that it was the one axis never varied.  It is now bounded from two
directions.  The front's response to the diffuse field, calibrated two ways from
the 584 A measurement, is 0.233-0.585 pc per unit fractional change in the total
diffuse band luminosity, so 0.18 pc requires a **31-77%** change in that
luminosity.  And the *entire* dynamic range of the diffuse treatment --- fully
explicit transport at 4.11970 pc against no explicit transport at all
(on-the-spot) at 4.18159 pc --- is **0.062 pc**, with Cloudy a further 0.124 pc
beyond the more favorable end.  Every conceivable transport variant lies inside
that bracket, so none of them reaches the reference.

Two further candidates from the previous revision have now been measured rather
than argued about, and both are bounded far below the gap: the He I ground-state
photoionization cross section (a ~70% change would be needed) and the He I
recombination rate (~23%), against a spread of at most 5% among the published
determinations of the recombination coefficient and 0.5% between the cross
section and its laboratory measurement.

What survives is not an explanation but a signature: the discrepancy is
He-specific and compounds with depth
([The radial signature](#the-radial-signature)).  **The residual is still
unexplained.**  What has changed is that the space of simple explanations has
been closed systematically, and what remains demands an He-specific effect at
the 15-20% level that is present in none of the rates, cross sections or fields
checked here.

## What is excluded, and how tightly

Each row is a measurement or a construction, with the bound it places on
`r(He0 = 0.5)`.  The residual to be closed is **0.18561 pc**.

| candidate | how it was bounded | bound |
| --- | --- | --- |
| Diffuse-field **emission** | four independent corrections to the emitted spectrum, one of them removing 9-11% of the strongest channel's energy | sum **-0.00306 pc**, wrong sign: 1.7% of the gap |
| Diffuse-field **transport** | front response calibrated two ways from the 584 A run: 0.233 pc (differential) and 0.585 pc (channel switched off) per unit fractional change in the diffuse band luminosity | closing the gap needs a **31-77%** change in that luminosity |
| The diffuse treatment **as a whole** | explicit transport (4.11970 pc) against on-the-spot (4.18159 pc), the two limits of the treatment | full range **0.062 pc**; Cloudy 0.124 pc outside its favorable end |
| `sigma_pi(He I)` | paired run at x0.97 | **+0.00799 pc**; **~70%** change needed |
| He I recombination rate | paired run at x0.97 | **+0.02410 pc**; **~23%** change needed |
| Dust | off by construction in the benchmark (`par%dust_model='none'`, `par%ion_add_dust=.false.`) | identically zero |
| Source photon budget | `Q(>13.598 eV)` and `Q(>24.587 eV)` against Cloudy c25 | agree to **1.00005** and **1.00003** |
| H ionization structure | `x(H0)` profile ratio to Cloudy across the ionized zone | **1-2%** at every radius |
| Cloudy's own resolution at the front | 40 zones between 4.0 and 4.6 pc, median width 0.0106 pc, maximum 0.039 pc | reference front resolved **17x** finer than the residual |

Earlier in the investigation the diagnostic energy binning, Monte Carlo noise,
spatial resolution and the H/He recombination rate set were excluded as well;
those measurements are under
[Effects excluded or bounded](#effects-excluded-or-bounded).

None of these bounds is an argument that the residual is small.  It is 0.18561
pc, about 4.4% of the front radius, and it is real.  They are the statement that
no single quantity in the list can produce it.

## The radial signature

The one positive result of the investigation is the *shape* of the discrepancy.
`x(H0)` agrees with Cloudy to 1-2% throughout the ionized zone, while the
`x(He0)` ratio grows monotonically outward:

| r [pc] | 1.6 | 2.0 | 2.5 | 3.0 | 3.5 | 3.8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `x(He0)` MoCHII/Cloudy | 1.007 | 1.015 | 1.031 | 1.059 | 1.109 | **1.196** |
| `x(H0)` MoCHII/Cloudy | 0.981 | 0.984 | 0.989 | 1.000 | 1.005 | 1.009 |

This is not a normalization offset.  A wrong rate, cross section or source
luminosity appears as a radius-independent factor; what is observed instead
starts at 0.7% and accumulates with depth to 20%, in helium only.

The accumulation is expected once something starts it.  He^0 carries most of the
opacity above 24.6 eV in this zone --- at 30 eV,
`kappa(He0)/[kappa(He0) + kappa(H0)]` is 0.843 at 1.6 pc and rises to 0.956 at
4.0 pc --- so a slight excess of He^0 attenuates the He-ionizing field slightly
more, which leaves more He^0 further out.  The He ionization structure amplifies
its own error, and the consequence is that **the seed is small**: 0.7% at
1.6 pc, well inside the accuracy of any individual rate.

Two things follow.  First, hunting for a large error in a single He I rate or
cross section is misdirected, and the sensitivity runs below confirm that
quantitatively.  Second, a budget statement: moving the front outward by 4.4% in
radius is 13.7% in volume, so a *uniform* explanation would require `Q(He0)`
higher by ~14% or `alpha_B(He I)` lower by ~20%.  `Q(He0)` is verified against
Cloudy to 0.003%, and no published `alpha_B(He I)` is 20% away from MoCHII's.
No uniform explanation exists, and what is needed instead is something acting
progressively on the He^0 population --- that is, inside the He ionization
network.

## Common diagnostic

All reported MoCHII front positions use volume-weighted spherical shells of
width 0.01 pc, followed by linear interpolation between shell centers.  The
diagnostic fields are `x_HeI = 0.5`, `x_c_stages(:,3) = 0.5` (C2+), and
`x_HI = 0.1`.  The analysis is implemented in
`tests/continuous_energy/analyze_hii40_residuals.py`.

| Case | r(He0 = 0.5) [pc] | r(C2+ = 0.5) [pc] | r(H0 = 0.1) [pc] | C2+(4.20 pc) | C2+(4.35 pc) | C2+(4.50 pc) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Cloudy c25 | 4.30531 | 4.22215 | 4.71311 | 0.54343 | 0.26062 | 0.20087 |
| MoCHII baseline: exponential continua, statistical cascade (L6) | 4.12276 | 4.05918 | 4.62707 | 0.25609 | 0.20980 | 0.15623 |
| Diffuse field disabled (OTS) | 4.18159 | 4.10080 | 4.47390 | 0.25048 | 0.16583 | 0.01359 |
| He I excited diffuse channel disabled | 4.02014 | 3.96893 | 4.51111 | 0.22819 | 0.17517 | 0.04835 |
| Alternate MOCREC atomic data | 4.11660 | 4.03970 | 4.62831 | 0.21084 | 0.17152 | 0.12684 |
| L7 spatial grid | 4.08411 | 4.05686 | 4.66714 | 0.25191 | 0.20884 | 0.15817 |
| Ground continua at threshold | 4.09112 | 4.02432 | 4.57811 | 0.24231 | 0.19657 | 0.13390 |
| Ground continua from Milne | 4.12339 | 4.05970 | 4.62628 | 0.25621 | 0.20976 | 0.15602 |
| + Drake two-photon shape | 4.12309 | 4.05954 | 4.62599 | 0.25618 | 0.20968 | 0.15589 |
| + Wood branch fractions | 4.12333 | 4.05969 | 4.62635 | 0.25621 | 0.20979 | 0.15608 |
| + He I 584 A conversion (current production) | 4.11970 | 4.05637 | 4.61633 | 0.25380 | 0.20723 | 0.15171 |

Thus the baseline He I 50% front is 0.18255 pc inward of Cloudy, and the
C2+ 50% front is 0.16297 pc inward; the current production state, with all
four diffuse-field corrections in place, sits at 0.18561 and 0.16578 pc --- that
is, marginally *further* from Cloudy than the baseline, because the fourth
correction moves inward by more than the first three move outward.
The co-location of these residuals is the main reason to prioritize the
He-ionizing diffuse field over C-specific photon sampling or cross sections.

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
not close the remaining gap: the He I front is still 0.18255 pc inside
Cloudy.

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
| r(H0 = 0.1) | 4.62707 | 4.62628 | −0.00079 | --- |

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
the packet count.  These are independent data --- Verner cross sections through
detailed balance on one side, fitted recombination coefficients on the other ---
and are reported rather than reconciled:

| continuum | 5000 K | 10000 K | 20000 K |
|---|---:|---:|---:|
| H I | 0.990 | 0.995 | 1.001 |
| He I | **1.034** | **1.028** | **1.020** |
| He II | 0.998 | 0.996 | 0.996 |

H I and He II agree to about 0.5%, which is the expected level for independent
fits.  He I disagrees by 2-3.4%.  The unexplained residual is the He I front,
so an internal inconsistency of that size in exactly the He I ground-continuum
data was worth pursuing.

It has since been pursued, and it is closed as a candidate: the inconsistency is
real and is between two independent atomic datasets rather than being an
implementation error, and the front is far too insensitive to it to matter.  See
[He I atomic data](#he-i-atomic-data-sensitivity-measured-inconsistency-located).

## The He I cascade: implemented, and it does not explain the residual either

The cascade was the first item on the earlier candidate list, and it has now
been done in three separately measured steps.  The third came out of doing the
first two, is an order of magnitude larger than either, and has its own section
below.

**Step 1, the two-photon shape.**  The `2^1S` continuum was sampled flat over
[13.598, 20.62] eV.  It now follows the tabulated Drake, Victor & Dalgarno
(1969) distribution in `data/HeI2phot.dat` --- the same table `nebcont_mod`
had been using all along, so the nebular continuum MoCHII wrote and the
packets it transported had disagreed about one process.  The fixed constant
0.56 H-ionizing photons per decay turns out to be right to 0.6% (0.5564), but
the flat surrogate overestimated the mean in-band photon energy by 6%
(17.109 against 16.110 eV) and the radiated in-band energy by 6.5%.

**Step 2, the branch fractions.**  The cascade was split by the AGN3
low-density statistical weights, which carry no temperature dependence.
`src/hei_cascade_mod.f90` now takes all three branches from one internally
consistent set of effective recombination coefficients --- Benjamin, Skillman
& Smits (1999) Table 1 with J. S. Mathis, in the form tabulated in K. Wood's
`ionize` (`probset.f`), and used with the same sum normalization by
CMacIonize, which cites Wood, Mathis & Ercolano (2004) section 3.3:

    alpha(2^3S) = 2.10e-13 (T/1e4)^-0.381
    alpha(2^1S) = 2.06e-14 (T/1e4)^-0.451
    alpha(2^1P) = 4.17e-14 (T/1e4)^-0.695     [cm^3 s^-1]

| branch | statistical | this set (8000 K) | 5000 K | 10^4 K | 1.5e4 K |
|---|---:|---:|---:|---:|---:|
| 2^3S | 0.750 | 0.762 | 0.741 | 0.771 | 0.787 |
| 2^1S | 0.083 | 0.076 | 0.076 | 0.076 | 0.075 |
| 2^1P | 0.167 | 0.162 | 0.183 | 0.153 | 0.138 |

This is a modest temperature-dependent correction to the statistical weights,
not a reordering of them.  The fractions are normalized by the sum of the
three coefficients rather than forced onto `alpha_B(He II)`, which is Wood's
own choice and made for the reason his source states --- the three fits have
unrelated temperature dependences.  Their sum is 0.767, 0.896, 0.970 and
1.134 of MoCHII's `alpha_B(He II)` at 5000, 8000, 10^4 and 1.5e4 K; that
shortfall is a statement about two independent atomic data sets and is
reported at setup rather than absorbed.  The fits are nebular, so they are
clamped outside 3e3-2e4 K.  The 2^3S coefficient is deliberately *not* shared
with `hei_metastable_mod`, whose Oklopcic & Hirata (2018) fit
`2.10e-13 (T/1e4)^-0.778` has the same prefactor but a steeper exponent (9%
higher at 8000 K): that one is an absolute production rate for the metastable
population, this one is a member of a set whose *ratios* are the branch
fractions, and mixing sources would break the consistency the fractions
depend on.

**Audit note: what this replaced, and why it was wrong.**  A first
implementation read `nebcont_mod`'s `6.23e-14 (T/1e4)^-0.827` as
`alpha_eff(2^1S)` and closed 2^1P against `alpha_B(He II)`, obtaining
0.746 : 0.224 : 0.031 at 8000 K, and argued from that inversion that case-B
trapping of the 584 A resonance line sends nearly all singlet recombinations
through the two-photon continuum.  The coefficient is not the 2^1S branch:
`2.06e-14 + 4.17e-14 = 6.23e-14` exactly, so it is the singlet *total*, as
MOCASSIN's `emission_mod.f90` states at the line it comes from ("assume all
HeI singlets finally end up in the 2^1S / use total recombination cofficient
to all singlets").  The closure derivation, the 0.224/0.031 fractions and the
inversion argument are withdrawn.  The observation behind the argument
survives and is now step 3 below: `6.23e-14` is the limit in which every
584 A decay converts, and MoCHII does not implement that conversion at all.

**What the channel radiates.**  For each case-B He II recombination
reaching an H-ionizing decay, channel 4 radiates

| configuration | in-band energy | in-band photons |
|---|---:|---:|
| statistical weights + flat two-photon | 19.200 eV | 0.9634 |
| after step 1 (Drake two-photon) | 19.149 eV (-0.27%) | 0.9630 |
| after step 2 (8000 K branch fractions) | 19.223 eV (+0.39%) | 0.9663 (+0.34%) |

**Measured effect.**  Step 1 moves the fronts inward, step 2 outward, and
both are of order 1e-4 pc against a 3.6e-5 pc seed-to-seed noise floor:

| step | dr(He0=0.5) | dr(C2+=0.5) | dr(H0=0.1) | line ratios |
|---|---:|---:|---:|---|
| Drake two-photon shape | -0.00030 | -0.00016 | -0.00029 | 0.10% median |
| Wood branch fractions | +0.00024 | +0.00015 | +0.00036 | 0.061% median |

Step 2 raises `L(Hbeta)` by 0.0121% and changes the line ratios by 0.061% in
the median, the largest being C I 609 um at -1.01% and C I 370 um at -0.85%
--- the expected sign for 0.34% more H-ionizing photons pushing the neutral
zone back.  The paired HII20 run responds in the same direction and more
weakly: `r(He0 = 0.5)` moves by -0.00001 pc after step 1 and +0.00006 pc
after step 2, `r(C2+ = 0.5)` is unchanged to 1e-5 pc, `L(Hbeta)` rises by
0.0027%, and the line ratios move by 0.016% in the median.  Its He II lines
are excluded from the HII20 statistics as usual,
being at 1e-8 of `L(Hbeta)` where the seed noise is 14%.
`tests/continuous_energy/check_hii_production.py` passes on both benchmarks.

## Step 3: the fate of the He I 584 A photon

This is the largest of the four corrections and the one that settles the
diagnosis.  It is planned in
[HEI_584_CONVERSION_PLAN.md](HEI_584_CONVERSION_PLAN.md), which also records
what each reference code does.

**The defect.**  About one sixth of every case-B He II recombination reaches
the ground state through `2^1P`, radiating He I Ly-alpha at 584 A = 21.22 eV.
MoCHII emitted that photon into transport as an H-ionizing packet in every
case.  But 584 A is a *resonance* line: the surrounding He^0 is optically
thick to it, so the photon scatters many times before anything irreversible
happens, and two irreversible outcomes compete during that walk --- absorption
by a hydrogen atom, which photoionizes it and deposits 7.62 eV, or collisional
transfer `2^1P -> 2^1S` in the last helium atom that absorbed it, after which
the decay is the two-photon continuum instead.  Emitting a 21.22 eV packet
unconditionally is one of the two limits of that competition, and on these
grids it is the wrong one.

**The implementation.**  `src/hei_cascade_mod.f90` gains the
escape-probability closure of Wood, Mathis & Ercolano (2004) section 3.3,

    hei_584_ionizes_H(T, xHI, xHeI) = [1 + 77 x(He0)/(sqrt(T) x(H0))]^-1

and `hei_branch_energies` in `src/diffuse_mod.f90` folds it into the branch
weights:

    e3s = f3s*19.82
    e1p = f1p*p584*21.22
    e1s = (f1s + f1p*(1 - p584))*e_2ph_per_decay

Because `p584` is a function of the leaf's own state, the split enters at build
time and **no new Sobol dimension is needed** --- the assignment (7 = channel,
8 = branch, 9 = two-photon energy) is unchanged, so every stored RQMC
comparison stays valid.  A converted decay radiates 20.62 eV as the two-photon
pair, not 21.22 eV; the 0.60 eV difference leaves as `2^1P -> 2^1S` at
2.06 micron, which is non-ionizing and correctly absent from the band.

`gamma_2q_HeI` in `src/nebcont_mod.f90` now uses the same split,
`alpha_2q = alpha(2^1S) + [1 - p584] alpha(2^1P)`.  **This closes a
long-standing internal inconsistency**: `nebcont_mod` had carried MOCASSIN's
`6.23e-14 (T/1e4)^-0.827`, which is the `p584 -> 0` limit, while `diffuse_mod`
sat at `p584 -> 1`, so the continuum MoCHII wrote and the packets it
transported disagreed about the fate of the same decay.  Both now describe one
nebula.

**The gate**, `tests/hei_cascade/check_hei_584.f90`, checks the closure against
the algebraically distinct forms in K. Wood's `ionize` (`mcionize.f`,
`h0find.f`) and in CMacIonize (`PhysicalDiffuseReemissionHandler`) and finds
agreement to **2.2e-16**; the limits `x(H0) -> 0` and `x(He0) -> 0` are exact;
`p584` is monotone in both fractions; forcing `p584 = 1` reproduces the step-2
output; and conversion strictly reduces the channel's in-band output.  All
three cascade gates pass, and
`tests/continuous_energy/check_hii_production.py` passes on both benchmarks.

**Measured effect.**  Paired 68-rank runs against step 2 at the same Sobol
seed, with the standard 0.01-pc volume-weighted shell diagnostic:

| | HII40 | HII20 |
|---|---:|---:|
| dr(He0 = 0.5) | **-0.00363** | -0.00022 |
| dr(C2+ = 0.5) | -0.00332 | -0.00001 |
| dr(H0 = 0.1) | **-0.01002** | (density bounded) |
| dL(Hbeta) | -0.1997% | -0.0704% |
| mean Te(x_HI < 0.95) | 9080.2 -> 9053.4 K | 8523.0 -> 8522.3 K |
| diffuse band luminosity | -1.56% | -0.22% |

In absolute terms `r(He0 = 0.5)` moves from 4.12333 to **4.11970** pc in
HII40 and from 1.23889 to **1.23867** pc in HII20.  This is the largest of the
four corrections by a wide margin: six times the Milne ground continua and
twelve times the step-2 branch fractions.

**The predictions, graded.**  Three predictions were recorded before the runs,
as the plan requires; two held and the one that failed is the informative one.

- **Sign --- correct.**  Both fronts moved inward, as the loss of channel-4
  photons requires.
- **Magnitude of order 0.01 pc --- correct.**  The HII40 H0 front moved
  -0.01002 pc and the He0 front -0.00363 pc.
- **"More in HII20 than in HII40" --- WRONG.**  HII40 moved sixteen times
  further.  The prediction was read off `p584` alone, which measures only what
  fraction of the `2^1P` branch converts: 99.6% in HII20 against 88% in HII40.
  Front displacement instead scales with the branch's **absolute** contribution
  to the ionizing budget.  Inverting the measured luminosity drops against the
  change in emission per recombination puts channel 4 at **17.6% of the HII40
  diffuse field and 2.0% of the HII20 one** --- a 20 kK spectrum makes very few
  photons above 24.6 eV, so He^+ is scarce and the channel feeding on its
  recombinations is small.  A larger converted fraction of a much smaller
  channel moves the front less.

  The lesson generalizes: a local branching probability predicts what happens
  per event, never how far a front moves.  The second factor is the channel's
  share of the ionizing budget, and that is set by the hardness of the source.

## What the four corrections together say

| correction | dr(He0 = 0.5) | Cloudy gap closed |
|---|---:|---:|
| Milne ground continua | +0.00063 | +0.35% |
| Drake two-photon shape | -0.00030 | -0.16% |
| Wood branch fractions | +0.00024 | +0.13% |
| He I 584 A conversion | **-0.00363** | **-1.99%** |
| **all four** | **-0.00306** | **-1.68%** |

`r(He0 = 0.5)` against the Cloudy c25 reference at 4.30531 pc:

| state | r(He0 = 0.5) | residual |
|---|---:|---:|
| baseline (exponential continua, statistical cascade) | 4.12276 | -0.18255 |
| after steps 1-2 and Milne | 4.12333 | -0.18199 |
| after the 584 A conversion (production) | 4.11970 | **-0.18561** |

**All four corrections are physically better, and together they move MoCHII
slightly away from Cloudy.**  That is not an argument against any of them.
Each reproduces an external reference implementation or an exact relation, and
the acceptance criterion in `CLAUDE.md` is physical correctness, not agreement
with a particular reference calculation.  Impact size is a tool for planning
and validating a change, never a vote on whether it belongs.

What the result does settle is the diagnosis, and it is the central conclusion
of this note: **the 0.18 pc residual is not in what the diffuse field emits.**
The first three corrections could be dismissed as sub-percent detail --- the
two-photon shape changed the channel's energy by -0.27%, the branch fractions
by +0.39%.  The fourth cannot: it removed 9-11% of the channel's energy and
7-8% of its photons, replaced a 21.22 eV line with a 13.6-20.6 eV continuum of
mean 16.1 eV, and still bought only 2% of the gap.  For comparison, switching
the whole He I cascade off moves `r(He0 = 0.5)` by -0.103 pc, so the channel is
large; it is the *refinements* to what it emits that are not.  The emission
axis is exhausted.

## How sensitive the front is to the diffuse field

The 584 A run also calibrates the front against the diffuse field as a whole,
which is what converts every later measurement into a bound.  Two calibrations
are available from it, and they bracket the answer:

| calibration | change in diffuse band luminosity | shift in `r(He0 = 0.5)` | implied slope |
| --- | ---: | ---: | ---: |
| differential: channel 4 lost 8.87% of its energy, and it is 17.6% of the diffuse field | -1.56% | -0.00363 pc | 0.233 pc |
| channel 4 switched off entirely | -17.6% | -0.103 pc | 0.585 pc |

The slope is in pc per unit fractional change in the total diffuse band
luminosity.  The differential figure is the local derivative and the switched-off
figure includes the nonlinear response, so the true slope for a large change is
between them.  Either way, **0.18 pc requires the diffuse band luminosity to
change by 31-77%.**  That is the yardstick used below: any candidate whose effect
on the diffuse field is a few percent cannot be the explanation, whatever its
sign.

## He I atomic data: sensitivity measured, inconsistency located

The `alpha_1` lead from the Milne work has been followed to its end.  It splits
into two questions --- is MoCHII's He I data wrong, and would it matter --- and
the answers are *no implementation error* and *no*.

**Would it matter: two paired sensitivity runs.**  HII40, 68 ranks, the same
Sobol seed as production, 38 gas iterations, outputs under
`results/continuous_energy/hei_sensitivity/`.  Each run scales one He I quantity
by 0.97 --- a 3% change, the size of the entire observed disagreement --- and the
shifts are measured with the standard 0.01-pc volume-weighted shell diagnostic:

| variant | dr(He0 = 0.5) | dr(C2+ = 0.5) | dr(H0 = 0.1) | pc per unit change | needed to close 0.18561 pc |
| --- | ---: | ---: | ---: | ---: | ---: |
| `sigma_pi(He I)` x 0.97 | **+0.00799** | +0.00518 | +0.00006 | 0.27 | **~70%** |
| He I recombination x 0.97 | **+0.02410** | +0.02677 | -0.00285 | 0.80 | **~23%** |

Both have the sign that helps, and the recombination rate is the stronger of the
two by a factor of three.  That ordering is expected: lowering
`sigma_pi(He I)` acts in two opposing ways --- less He^0 is ionized locally, but
the He-ionizing photons penetrate further --- and they partly cancel, which the
near-zero `dr(H0 = 0.1)` of that run also shows, whereas lowering the
recombination rate reduces He^0 everywhere at once.  But 3% buys
4.3% and 13.0% of the gap respectively, and the changes required to close it,
70% and 23%, are an order of magnitude beyond the spread of the published data:
the three determinations of `alpha_1(He I)` below differ by at most 5%, and
`alpha_B(He I)` differs between MoCHII and BSS99 by 2.9%.  **He I atomic data
cannot close the gap.**

**Is MoCHII's cross section wrong: no.**  `sigma_HeI` reproduces Verner's own
`phfit2.f` to a ratio of **1.00000** at every energy between 24.6 and 105 eV,
with a maximum absolute difference of 0.000%, so there is no implementation
error to find.  Its threshold value is 7.436 Mb against the measured
7.40 +/- 0.07 Mb of Samson et al. (1994) --- 0.5% high and inside the
experimental uncertainty.  The normalization is right, so any disagreement can
only be in the energy dependence.

**Which side of the inconsistency is the outlier: the Milne integral.**  Because
`recomb_mod` forms `alphaB_HeII = alphaA_HeII - alpha1_HeII`, the "rates" side of
the `milne_setup` cross-check *is* `alpha1_HeII` itself, an `rr_mao` fit, so the
comparison is a direct one between two determinations of the same coefficient.
A third, independent determination breaks the tie: BSS99 with J. S. Mathis, in
the form carried by K. Wood's `ionize`,
`alpha_1(He I) = 1.54e-13 (T/1e4)^-0.486`.

| T | BSS99 / Wood | `rr_mao` (MoCHII) | Milne integral (Verner sigma) |
| --- | ---: | ---: | ---: |
| 5000 K | 2.157 | 2.184 (+1.3%) | 2.259 (**+4.8%**) |
| 10^4 K | 1.540 | 1.565 (+1.6%) | 1.609 (**+4.5%**) |
| 2x10^4 K | 1.100 | 1.128 (+2.6%) | 1.151 (**+4.7%**) |

Units are 1e-13 cm^3 s^-1 and the percentages are against the BSS99 column.  The
two recombination determinations cluster within 2.6% of each other while the
Milne integral sits 4.5-4.8% above both, so it is the Milne integral that stands
apart --- and He I is the only non-hydrogenic case among the three continua,
which is consistent with the residual being in the Verner He I fit's energy
dependence rather than in either recombination coefficient.

The code reports the discrepancy at setup and rescales neither side.  That is the
correct handling: this is an inconsistency between two independent atomic
datasets, not a defect, and absorbing it into one of them would hide a real
uncertainty in the data.

## The source photon budget, and dust

The ionizing budget entering the calculation is verified directly against the
reference.  For the benchmark's 40 kK blackbody (`par%tstar=4.0e4`,
`par%luminosity=1.2734e39` over 13.598-100 eV):

| | MoCHII | Cloudy c25.00 | ratio |
| --- | ---: | ---: | ---: |
| `Q(>13.598 eV)` | 4.2656e49 s^-1 | 4.2653e49 (`Q > 1 Ryd`) | 1.00005 |
| `Q(>24.587 eV)` | 4.6082e48 s^-1 | 4.6081e48 (`Q > 1.8 Ryd`) | 1.00003 |

The two codes are illuminating the same nebula with the same number of photons
to 0.005%, and in particular the He-ionizing budget --- the quantity that would
have to rise by ~14% for a uniform explanation of the front --- agrees to
0.003%.  The source normalization, the spectral shape at the two thresholds, and
the band definition are all excluded.

Dust is excluded by construction rather than by measurement: the Lexington
benchmarks run with `par%dust_model='none'` and `par%ion_add_dust=.false.`, so
there is no grain opacity competing for ionizing photons in either calculation
being compared.

## Next decision

Excluded: diagnostic binning, Monte Carlo noise, spatial resolution, the H/He
recombination rate set, the emitted diffuse spectrum as a whole, the transport of
the diffuse photons, the He I ground-state photoionization cross section, the
He I recombination rate, dust, the source photon budget, the H ionization
structure, and the resolution of the reference calculation at the front.  The
bounds are tabulated in
[What is excluded, and how tightly](#what-is-excluded-and-how-tightly).

What remains, in the order the evidence supports:

1. **The He ionization network itself** --- specifically, what Cloudy's resolved
   He I model atom contains that MoCHII's effective treatment does not.  This is
   where the radial signature points: the discrepancy is He-specific, absent
   from hydrogen, seeded at the 0.7% level and amplified by He^0's own dominance
   of the opacity above 24.6 eV, so it is a property of how the He^0 population
   is computed rather than of any single rate feeding that computation.
   Confirming or refuting it requires reading Cloudy's He I implementation and
   identifying which of its channels MoCHII represents by an effective
   coefficient.  That is a **separate task** and is the open item to record.
2. **Collisional redistribution among the n = 2 terms**, which is the most
   concrete instance of item 1 that MoCHII could already act on.  The branch
   fractions are effective coefficients in the low-density limit and `p584` is an
   escape-probability closure; `hei_metastable_mod` already carries the
   `q_31a`/`q_31b` rates that would couple the terms at finite density.  It is
   density dependent, so it is unlikely to dominate at the benchmark densities,
   but it is testable without leaving the code.
3. **A difference outside the He treatment entirely** --- the metal opacity at
   the front, or a Cloudy-side modeling choice not yet identified.

The quantitative statement to keep in view while pursuing any of these: closing
0.18561 pc means moving the front out by 4.4% in radius, hence 13.7% in volume,
which a uniform explanation would have to supply as `Q(He0)` higher by ~14% or
`alpha_B(He I)` lower by ~20%.  `Q(He0)` agrees with Cloudy to 0.003% and no
published `alpha_B(He I)` is 20% from MoCHII's, so **neither is supported by the
data**.  The explanation, if there is one on the MoCHII side, has to be
progressive in depth rather than a change of scale.
