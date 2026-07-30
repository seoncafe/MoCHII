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
how it is transported --- and it is not in the metal opacity, nor in MoCHII's
H/He recombination coefficients.**

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
recombination rate (~23%).  The cross section reproduces Verner's own `phfit2` to
a ratio of 1.00000 and the laboratory threshold value to 0.5%.  The recombination
rate is the one quantity in the network on which independent codes visibly
disagree --- four of them span 11-13% in `alpha_A(He I)`
([The He ionization network, code by code](#the-he-ionization-network-code-by-code))
--- but MoCHII's own value sits within 0.3% of Cloudy's, so MoCHII's share of
that spread is nothing like the 23% required.

*Three results added in this revision.*

1. **A real defect was found in `alpha_1`, the ground-level direct-capture
   coefficient, and fixed.**  It is not the front's explanation, but it was
   physics: the code subtracted a level-resolved fit from an unrelated total.
   Correcting it moves the front outward by **+0.00811 pc**, so the residual is
   now **-0.17750 pc** ([The `alpha_1` defect](#the-alpha_1-defect-found-and-fixed)).
   Finding it also **withdraws** a conclusion recorded earlier the same day,
   that the Milne integral was the outlier among three determinations of
   `alpha_1(He I)`.
2. **The metal opacity at the front has been measured** for the first time, and
   it is not a trace effect: removing it moves the front inward by 0.615 pc.
   But the response is strongly asymmetric and saturating --- adding opacity
   gains only 0.103 pc per unit --- so closing the residual would take **+180%**,
   and an independent Opacity Project comparison says MoCHII's cross sections
   are if anything 2-14% too *high*
   ([Metal opacity](#metal-opacity-measured-and-not-the-explanation-either)).
   The earlier statement that a ~30% error in metal opacity would be worth
   ~0.18 pc is **withdrawn**: 30% is worth 0.031-0.072 pc.
3. **The He recombination network has been compared across four codes.**  No
   code truncates the recombination cascade --- all four sum over every level,
   and MoCHII's total agrees with Cloudy's explicit level-by-level sum to 0.3%.
   What the comparison isolates instead is two-electron atomic structure: the
   one-electron species agree to 0.7% (H I) and 0.3% (He II), where an analytic
   solution leaves nothing to disagree about, while He I spreads by 11-13% ---
   three determinations clustered within 2.5% and Verner & Ferland (1996), the
   set MOCASSIN uses, alone 11% high
   ([The He ionization network](#the-he-ionization-network-code-by-code)).

What survives is not an explanation but a signature: the discrepancy is
He-specific and compounds with depth
([The radial signature](#the-radial-signature)).  **The residual is still
unexplained.**  What has changed is that the space of *single-parameter*
explanations is now closed by measurement rather than by argument.  What remains
untested is a combination of several few-percent differences acting together
through that same compounding --- which is exactly what no single-parameter
experiment can reproduce.

## What is excluded, and how tightly

Each row is a measurement or a construction, with the bound it places on
`r(He0 = 0.5)`.  The residual to be closed is **0.17750 pc** --- it was
0.18561 pc until the `alpha_1` correction below moved the front outward by
0.00811 pc.

| candidate | how it was bounded | bound |
| --- | --- | --- |
| Diffuse-field **emission** | four independent corrections to the emitted spectrum, one of them removing 9-11% of the strongest channel's energy | sum **-0.00306 pc**, wrong sign: 1.7% of the gap |
| Diffuse-field **transport** | front response calibrated two ways from the 584 A run: 0.233 pc (differential) and 0.585 pc (channel switched off) per unit fractional change in the diffuse band luminosity | closing the gap needs a **31-77%** change in that luminosity |
| The diffuse treatment **as a whole** | explicit transport (4.11970 pc) against on-the-spot (4.18159 pc), the two limits of the treatment | full range **0.062 pc**; Cloudy 0.124 pc outside its favorable end |
| `sigma_pi(He I)` | paired run at x0.97 | **+0.00799 pc**, i.e. 0.27 pc per unit; **~70%** change needed |
| He I recombination rate | paired run at x0.97 | **+0.02410 pc**, i.e. 0.80 pc per unit; **~23%** change needed |
| `alpha_A`, `alpha_1`, `alpha_B` (He I) | each against Cloudy c25.00, which derives them level by level from its own cross sections through the Milne relation | **0.3%, 0.34%, 1.9%**; at 0.80 pc per unit that is at most 0.015 pc |
| **Metal opacity** | paired runs at x0.70 and x1.30, plus off and metals-off bounds | **0.103 pc per unit upward** (0.615 pc downward); **~180%** increase needed, and saturating |
| Metal recombination continua as a diffuse source | summed ground-continuum emission above 24.587 eV from the nine metal stages that produce it | **0.29%** of the He I ground-continuum channel, ~0.001 pc |
| He I **excited-level** photoionization | 2^3S population and `Q(>4.767 eV)` for the benchmark spectrum, with Cloudy's own 2^3S threshold cross section | **0.2-0.6%** extra ionization; structurally absent from MoCHII (band floor 13.598 eV) |
| Dust | off by construction in the benchmark (`par%dust_model='none'`, `par%ion_add_dust=.false.`) | identically zero |
| Source photon budget | `Q(>13.598 eV)` and `Q(>24.587 eV)` against Cloudy c25 | agree to **1.00005** and **1.00003** |
| H ionization structure | `x(H0)` profile ratio to Cloudy across the ionized zone | **1-2%** at every radius |
| Cloudy's own resolution at the front | 40 zones between 4.0 and 4.6 pc, median width 0.0106 pc, maximum 0.039 pc | reference front resolved **17x** finer than the residual |

Earlier in the investigation the diagnostic energy binning, Monte Carlo noise,
spatial resolution and the H/He recombination rate set were excluded as well;
those measurements are under
[Effects excluded or bounded](#effects-excluded-or-bounded).

None of these bounds is an argument that the residual is small.  It is 0.17750
pc, about 4.1% of the front radius, and it is real.  They are the statement that
no single quantity in the list can produce it --- not that the list is complete.
The one hypothesis left standing is a *combination* of several of these
few-percent differences acting together through the depth compounding of
[the radial signature](#the-radial-signature), and by construction no
single-parameter run can either produce it or exclude it.

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
quantitatively.  Second, a budget statement: moving the front outward by 4.1% in
radius is 12.9% in volume, so a *uniform* explanation would require `Q(He0)`
higher by ~13% or `alpha_B(He I)` lower by a comparable amount.  `Q(He0)` is
verified against Cloudy to 0.003%.  The `alpha_B` side is the more interesting
of the two now that four codes have been compared: the code-to-code spread in
`alpha_A(He I)` is 11-13%, which at the measured 0.80 pc per unit is worth
0.09 pc --- half the residual.  But MoCHII sits within 0.3% of Cloudy on that
coefficient, so the spread is real and MoCHII's share of it is not.  No uniform
explanation exists on MoCHII's side, and what is needed instead is something
acting progressively on the He^0 population --- that is, inside the He
ionization network.

## Common diagnostic

All reported MoCHII front positions use volume-weighted spherical shells of
width 0.01 pc, followed by linear interpolation between shell centers.  The
diagnostic fields are `x_HeI = 0.5`, `x_c_stages(:,3) = 0.5` (C2+), and
`x_HI = 0.1`.  The analysis is implemented in
`tests/continuous_energy/analyze_hii40_residuals.py`.

**A volume-convention error in that script, found and corrected.**  `LeafSize`
in the gas output is the leaf **full** width --- `src/gas_rates_mod.f90:408`
writes exactly that in the EXTNAME comment --- and the script formed the cell
volume as `(2*size)**3`, eight times too large.  Because the volume enters only
as a *relative* weight among the leaves inside one 0.01-pc shell, a constant
factor of eight cancels identically, so **every front position ever reported
from this script is unaffected**; re-running after the correction reproduces all
of them to the last digit.  The expression is now `size**3` with a comment.
`paper/tables/make_wood_tables.py:259` already used `(size * dist)**3`
correctly, `tools/python/mochii_output.py` documents the convention correctly,
and the g1/g4 gates derive their own half-widths from the grid geometry, so
nothing else carried the error.  It mattered only once the same volumes were
used for an *absolute* quantity, the photon budget below.

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
| + He I 584 A conversion | 4.11970 | 4.05637 | 4.61633 | 0.25380 | 0.20723 | 0.15171 |

Thus the baseline He I 50% front is 0.18255 pc inward of Cloudy, and the
C2+ 50% front is 0.16297 pc inward; with all four diffuse-field corrections in
place they sit at 0.18561 and 0.16578 pc --- that is, marginally *further* from
Cloudy than the baseline, because the fourth correction moves inward by more than
the first three move outward.
The co-location of these residuals is the main reason to prioritize the
He-ionizing diffuse field over C-specific photon sampling or cross sections.

The **current production state** adds the Milne-derived `alpha_1` of
[the section below](#the-alpha_1-defect-found-and-fixed), which moves
`r(He0 = 0.5)` from 4.11970 to **4.12781 pc** and leaves the residual at
**0.17750 pc**.  That change was measured as a paired difference against the row
above, so the remaining columns are not re-tabulated here.

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

It has since been pursued to the end, and it turned out to be **a defect in
MoCHII, not an inconsistency between datasets**: `alpha_1` was being subtracted
out of a total that came from a different determination.  The lead did not
explain the front, but it was worth following for its own sake.  The diagnosis,
the fix, and the withdrawal of the first (wrong) reading of this table are in
[The `alpha_1` defect](#the-alpha_1-defect-found-and-fixed).  After the fix this
cross-check reports 0.01% instead of 3.44%.

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
| after the 584 A conversion | 4.11970 | -0.18561 |
| after the Milne-derived `alpha_1` (production) | 4.12781 | **-0.17750** |

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

## He I atomic data: sensitivity measured, defect located

The `alpha_1` lead from the Milne work has been followed to its end.  It splits
into two questions --- is MoCHII's He I data wrong, and would it matter --- and
the answers are **yes, in `alpha_1`** and *no, not enough to matter for the
front*.  The sensitivity below is what settles the second question, and it was
measured first; the defect is [the next section](#the-alpha_1-defect-found-and-fixed).

**Would it matter: two paired sensitivity runs.**  HII40, 68 ranks, the same
Sobol seed as production, 38 gas iterations, outputs under
`results/continuous_energy/hei_sensitivity/`.  Each run scales one He I quantity
by 0.97 --- a 3% change, the size of the entire observed disagreement --- and the
shifts are measured with the standard 0.01-pc volume-weighted shell diagnostic:

| variant | dr(He0 = 0.5) | dr(C2+ = 0.5) | dr(H0 = 0.1) | pc per unit change | needed to close 0.17750 pc |
| --- | ---: | ---: | ---: | ---: | ---: |
| `sigma_pi(He I)` x 0.97 | **+0.00799** | +0.00518 | +0.00006 | 0.27 | **~70%** |
| He I recombination x 0.97 | **+0.02410** | +0.02677 | -0.00285 | 0.80 | **~23%** |

Both have the sign that helps, and the recombination rate is the stronger of the
two by a factor of three.  That ordering is expected: lowering
`sigma_pi(He I)` acts in two opposing ways --- less He^0 is ionized locally, but
the He-ionizing photons penetrate further --- and they partly cancel, which the
near-zero `dr(H0 = 0.1)` of that run also shows, whereas lowering the
recombination rate reduces He^0 everywhere at once.  But 3% buys 4.5% and 13.6%
of the gap respectively, and the changes required to close it are 70% and 23%.
For the cross section that is far outside anything defensible --- `sigma_HeI`
reproduces its own source exactly and the laboratory threshold value to 0.5%.
For the recombination rate the margin is narrower than it looked while only
fitted coefficients were being compared: the code-to-code spread in
`alpha_A(He I)` is **11-13%**, about half of the 23% required, and MoCHII lies
within 0.3% of Cloudy inside that spread
([The He ionization network, code by code](#the-he-ionization-network-code-by-code)).
**He I atomic data cannot close the gap** --- with the caveat that the
recombination rate is the one entry on the excluded list whose *external*
uncertainty reaches the same order as the requirement, even though MoCHII's own
value does not.

**Is MoCHII's cross section wrong: no.**  `sigma_HeI` reproduces Verner's own
`phfit2.f` to a ratio of **1.00000** at every energy between 24.6 and 105 eV,
with a maximum absolute difference of 0.000%, so there is no implementation
error to find.  Its threshold value is 7.436 Mb against the measured
7.40 +/- 0.07 Mb of Samson et al. (1994) --- 0.5% high and inside the
experimental uncertainty.  The next section confirms from the other side that the
cross section was never the problem: that same Verner cross section, put through
the Milne relation, reproduces Cloudy's independent level-resolved
`alpha_1(He I)` to 0.34%.

**WITHDRAWN: "which side of the inconsistency is the outlier: the Milne
integral."**  Earlier on 2026-07-30 this section concluded that the Milne
integral was the outlier, on the following evidence.  A third determination of
`alpha_1(He I)` was brought in --- BSS99 with J. S. Mathis, in the form carried
by K. Wood's `ionize`, `alpha_1(He I) = 1.54e-13 (T/1e4)^-0.486` --- and compared
with the `rr_mao` fit in `recomb_mod` and with the Milne integral over the Verner
cross section:

| T | BSS99 / Wood | `rr_mao` (MoCHII) | Milne integral (Verner sigma) |
| --- | ---: | ---: | ---: |
| 5000 K | 2.157 | 2.184 (+1.3%) | 2.259 (**+4.8%**) |
| 10^4 K | 1.540 | 1.565 (+1.6%) | 1.609 (**+4.5%**) |
| 2x10^4 K | 1.100 | 1.128 (+2.6%) | 1.151 (**+4.7%**) |

Units are 1e-13 cm^3 s^-1 and the percentages are against the BSS99 column.  The
arithmetic is correct and the table stands.  The **conclusion drawn from it does
not**: because the two fitted coefficients clustered within 2.6% of each other
and the Milne integral sat 4.5-4.8% above both, the Milne integral was called
the outlier and the disagreement was attributed to the energy dependence of the
Verner He I fit.

That was a majority vote among three numbers, and it was wrong.  Cloudy c25.00
computes the same coefficient a fourth way --- level by level from *its own*
photoionization cross sections through the exact Milne relation --- and lands
within **0.34%** of MoCHII's Milne integral, 1.6032e-13 against 1.60861e-13 at
10^4 K.  Two fits agreeing with each other is not independent confirmation when
both can be low for the same reason, and here they were: **the two fits are the
low side, not the Milne integral.**  What made the misreading possible was
treating "two out of three" as evidence, with no determination in the comparison
that was derived from first principles rather than fitted.

The former closing paragraph of this section --- "the code reports the
discrepancy at setup and rescales neither side; this is an inconsistency between
two independent atomic datasets, not a defect" --- is withdrawn with it.  It was
a defect, and it is fixed next.

## The `alpha_1` defect, found and fixed

### What was wrong

`src/recomb_mod.f90` took the ground-level direct-capture rate `alpha_1` from the
Mao & Kaastra (2016) **level-resolved** fit (`rr_mao`), took `alpha_A` from the
**Badnell/Strathclyde TAMOC total** (`rr_badnell` + DR), and then defined

    alphaB_* = alphaA_* - alpha1_*

That subtraction crosses two unrelated determinations, so any offset between them
appears in `alpha_B` --- and `alpha_A - alpha_B` is what sets the number of
diffuse ground-continuum packets.  A source-mixing structure of this kind is
wrong independently of how large the offset happens to be.

**Cloudy c25.00 settles which side is right.**  Cloudy derives He-like and H-like
recombination level by level from its own photoionization cross sections through
the exact Milne relation (`~/CLOUDY/c25.00/source/iso_radiative_recomb.cpp`; all
Cloudy paths below are relative to that tree) and caches the result in
`data/he_iso_recomb.dat` and `data/h_iso_recomb.dat`.  In those files
row 0 of each element is the ground level and the **last row is the total** ---
verified on hydrogen, where row 0 = 1.5840e-13 and the last row = 4.1642e-13
against the known `alpha_1s` = 1.58e-13 and `alpha_A` = 4.18e-13.

At 10^4 K, in cm^3 s^-1:

| | MoCHII `alpha_A` | Cloudy `alpha_A` | MoCHII `alpha_1` (Mao) | Cloudy `alpha_1` | MoCHII Milne integral |
|---|---:|---:|---:|---:|---:|
| H I | 4.19330e-13 | 4.1642e-13 | 1.58733e-13 | 1.5840e-13 | 1.57902e-13 |
| He I | 4.37188e-13 | 4.3583e-13 | **1.56500e-13** | **1.6032e-13** | **1.60861e-13** |
| He II | 2.18953e-12 | 2.1817e-12 | 6.53656e-13 | 6.5145e-13 | 6.51162e-13 |

Badnell's totals are accurate to 0.3-0.7% on all three, and the Mao `n = 1`
values to 0.2-0.3% for H I and He II --- but **He I is 2.38% low**, and the
subtraction propagates that into `alpha_B(He I)` as **+1.88%**.  MoCHII's own
Milne integral matches Cloudy to 0.3% / 0.34% / 0.04%.  The Milne relation is
exact and the two codes' He^0 ground cross sections agree to the 0.4% implied by
that agreement, so the Milne route is the correct one.

### The fix

`tools/fitting/fit_alpha1_milne.py` (new) evaluates the Milne relation on **the
cross sections the transport actually absorbs with** --- exact hydrogenic for
H I and He II (`sigma_hydrogenic`, Z = 1 and 2; note `sigma_HI` in
`src/photo_xsec.f90` is exact hydrogenic, not VFKY96) and VFKY96 for He I --- at
241 temperatures over 10^3-10^5 K, and refits the existing seven-coefficient
`rr_mao` functional form, so no Fortran structure changes.  New coefficients
`(a0, b0, c0, a1, b1, a2, b2)` in `src/recomb_mod.f90`:

| function | coefficients |
|---|---|
| `alpha1_HII` | 1.87007e-2, 1.29077, -2.17465e-3, 1.29189e1, 8.05060e-1, 8.44399e-2, 5.93940e-4 |
| `alpha1_HeII` | 3.26862e-3, 6.97361e-1, 3.14907e-4, 5.76341e1, 1.18534, 2.58417e1, 9.90134e-1 |
| `alpha1_HeIII` | 2.61525e-1, 1.39527, -3.22079e-4, 6.40708e1, 8.97582e-1, 5.01551e-1, 1.21533e-3 |

Fit residual, max `|fit/Milne - 1|`: 0.019% / 0.001% / 0.006% over 10^3-10^5 K,
and 0.0075% / 0.0004% / 0.0021% over the production range 3e3-3e4 K.  At 10^4 K
the fitted `alpha_1` is 1.57910e-13 / 1.60861e-13 / 6.51172e-13, i.e. **within
0.35% of Cloudy on all three channels**, with He I corrected from -2.38% to
+0.34%.

**The code's own cross-check confirms it.**  `milne_setup` reported "largest
`alpha_1` disagreement 3.44%" before the fix and **0.01%** after: the Milne
integral and the `alpha_A - alpha_B` that sets the packet count now come from one
set of cross sections, which is the whole point of the change.

**The gate.**  `tests/milne/check_alpha1.f90` (new, wired into
`tests/milne/run_check.sh`) re-implements the Milne quadrature independently from
`photo_xsec`, checks the three channels at 81 temperatures over 10^3-10^5 K to
1e-3, checks `alpha_1 < alpha_A` everywhere (max `alpha_1/alpha_A` = 0.560 /
0.496 / 0.445), and checks all three against the hardcoded Cloudy c25.00 values
to 1%.  Reinstating the old Mao He I coefficients makes it fail on both counts
(3.886e-2 from the quadrature, 2.383% from Cloudy) and it passes again on
reverting, so the defect cannot be reintroduced silently.

The namelist value was `'badnell_mao'` and is now **`'badnell_milne'`**, because
`alpha_1` is no longer Mao data.  The first reason recorded here for keeping the
old name --- that renaming it would break existing input files and examples ---
was wrong and is withdrawn: no input file in the repository sets `recomb_model`
at all, the value is reached as a default, so the rename touched twenty mentions
and needed no namelist migration.  Only the fit FORM of `alpha_1` remains Mao &
Kaastra's, which is why `recomb_mod` still calls it `rr_mao`.

### Measured effect, and the prediction graded

Recorded before the runs: *"`alpha_1(He I)` rises 2.4%, so channel 2 emits more;
`alpha_B(He I)` falls 1.9%, so there is less case-B recombination.  Both push the
front outward, toward Cloudy, by of order 0.015-0.02 pc."*

Measured on paired HII40/HII20 runs at 68 ranks with the same Sobol seed,
converged in the same 38 and 27 gas iterations as before:

| | HII40 | HII20 |
|---|---:|---:|
| `dr(He0 = 0.5)` | **+0.00811 pc** | **+0.00210 pc** |
| Cloudy residual | -0.18561 -> **-0.17750 pc** | --- |
| `dL(Hbeta)` | -0.3713% | -0.4159% |
| mean Te(`x_HI` < 0.95) | 9053.4 -> 9051.7 K | 8522.3 -> 8528.7 K |

`tests/continuous_energy/check_hii_production.py` passes on both, with energy
closure 1.5e-16 and 0.0.

- **Sign: correct.**  Outward, toward Cloudy, on both benchmarks.
- **Magnitude: over-predicted by about a factor of two.**  +0.0081 pc against a
  predicted 0.015-0.02 pc, so the measurement is about half the low end of the
  range.  Better than the metal prediction below, but still high.

**`L(Hbeta)` moves the wrong way, and that has to be said plainly.**  It falls
0.37-0.42%, taking HII40 from about 3% below the published 2.01-2.10e37 to
slightly further below.  The likely route is channel 1: `alpha_1(H I)` drops
0.52% from the Mao value to the Milne one, so the H I ground continuum recycles
marginally fewer H-ionizing photons.  This is **not** an argument against the
change.  An exact relation was being violated by a source-mixing subtraction; the
acceptance criterion is physical correctness, and the direction in which one
observable happens to move is not a vote on it.  It is recorded here rather than
buried because the next person to look at `L(Hbeta)` needs to know where 0.4% of
it went.

Production, `paper/tables/wood_hii{40,20}.tex` and
`paper/figures/wood_hii{40,20}.pdf` are regenerated from these runs; the
superseded set is preserved under
`results/continuous_energy/superseded_step3/`.  Benchmark values:

| | before | now | published Lexington (9 codes) | Cloudy c25.00 |
|---|---:|---:|---|---|
| HII40 `<T[NpNe]>` | 8162 K | **8169 K** | 7720-8199, **inside** | 8210 (-0.50%) |
| HII20 `<T[NpNe]>` | 6917 K | **6925 K** | 6402-6749, above | 6910 (**+0.22%**) |
| HII40 T_inner | 7535 K | 7536 K | --- | --- |
| HII20 T_inner | 6936 K | 6938 K | --- | --- |
| HII40 `L(Hbeta)` | 1.94e37 | 1.93e37 | 2.01-2.10e37 | --- |
| HII20 `L(Hbeta)` | 4.70e36 | 4.68e36 | 4.83-5.04e36 | --- |
| HII40 R_out | 1.44e19 cm | 1.44e19 cm | --- | --- |
| HII20 R_out | 8.89e18 cm | 8.88e18 cm | --- | --- |

The HII20 verdict is unchanged in kind: above the 1995-era published range, but
so is Cloudy c25.00, and MoCHII now agrees with that modern reference to 0.22%.

### Scope: only He I was defective

The same source-mixing structure could in principle affect any species for which
MoCHII carries both a total and a ground-level coefficient.  It does not: the
metal element files carry **totals only**, no ground-level entry, which is
correct because metals do not feed the diffuse field.  The defect was therefore
confined to H I, He I and He II, and of those only He I was outside its own
uncertainty.  The sweep that established this covered every registry entry as
well; see
[the registry-wide audit](#the-he-ionization-network-code-by-code) below.

## Metal opacity: measured, and not the explanation either

### What "trace" actually means here

`docs/PLAN.md:315` lists "metals feeding back on `n_e` / opacity (trace
approximation)" as a non-goal, but the benchmark inputs run with
`par%metal_ne = .true.`, `par%ion_metal_abs = .true.` and
`par%metal_heat = .true.`, so both feedbacks are on.  `species_opacity_add` sums
over **all** ionization stages with VFKY96 cross sections and self-consistent
stage fractions, iterated with the opacity feedback --- not a trace treatment in
any operational sense.  Those cross sections were checked against Verner's own
`phfit2.f` for the eight ions that carry the He-ionizing window (C I, C II,
C III, N II, O II, Ne I, S II, S III) and agree to **1e-6 relative**, which is
the single-precision round-off of `phfit2` itself.

What is genuinely absent is **metal recombination continua as a diffuse source**:
the diffuse field has exactly four channels (`dif_ch(4, nleaf)`), all H/He.
Summing the metal ground-state recombination continua whose photons lie above
24.587 eV (C IV, N III, N IV, O III, O IV, Ne III, Ne IV, S III, S IV) gives
**0.29%** of the He I ground-continuum channel, about 0.001 pc.  Metals as
secondary-ionization absorbers are also absent, but `par%use_sec_ion` is off in
the benchmark.

Two code comments --- at the `ion_metal_abs` declaration in `src/define.f90` and
in `src/species_mod.f90` --- used to justify the metal opacity as *"Negligible
next to H/He above 13.6 eV."*  **Measurement contradicts that**: removing it
moves the front 0.615 pc.  That justification has been withdrawn in the code;
both sites now carry the measured statement instead --- the 0.615 pc, the
hardening mechanism below, and the asymmetry (0.615 pc per unit removed against
0.103 pc per unit added).

### The measured lever

HII40, 68 ranks, same Sobol seed, standard 0.01-pc shell diagnostic:

| variant | r(He0 = 0.5) [pc] | dr vs production | pc per unit fractional change |
|---|---:|---:|---:|
| production | 4.11970 | --- | --- |
| metal `sigma` x 0.70 | 4.04764 | **-0.07206** | 0.240 |
| metal `sigma` x 1.30 | 4.15063 | **+0.03093** | 0.103 |
| metal absorption off | 3.50470 | -0.61500 | 0.615 |
| metals off entirely | 1.93335 | -2.18635 | (also removes all metal cooling; Te 8162 -> 18681 K) |

**The response is strongly asymmetric and saturating.**  Removing metal opacity
costs 0.615 pc per unit; adding it gains only 0.103 pc per unit.  Extrapolating
the upward branch, closing 0.17750 pc needs **+180% metal opacity (x2.8)**, and
the saturation means even more than that.  Not available.

**The mechanism is spectral hardening.**  Metal stages with thresholds inside
21-40 eV (Ne I 21.56, S II 23.33, C II 24.38, N II 29.60, S III 34.83, O II
35.12 eV) preferentially remove the *softest* He-ionizing photons.
`sigma(He I)/sigma(H I)` rises steeply with energy --- 6.1 at 25 eV, 9.0 at 35,
13.9 at 55, 19.9 at 95 eV --- so the hardened surviving field is taken up by
He^0 rather than H^0 and the He^+ zone extends.  Remove the metal opacity and the
field stays soft, H^0 takes a larger share, and the zone collapses.  A thermal or
hydrogen-structure explanation is ruled out by the same runs: between production
and `ion_metal_abs = .false.`, `r(H0 = 0.1)` moves from 4.61633 to 4.61551 pc and
the mean Te from 8162 to 8183 K, both negligible.

### Cloudy's metal photoionization, read from its source

Cloudy's base is the same Verner `phfit`, but `source/opacity_createall.cpp` adds
channels MoCHII does not have:

| Cloudy | MoCHII |
|---|---|
| N^+ **excited-level** photoionization (Henry 1970; 9e-18, `nu^-1.75`) | absent |
| O I valence replaced wholesale by `ofit` when PHFIT96 is in use | Verner as-is |
| O I 1S excited level (4.64e-18; "OP data very sparse") | absent |
| O^2+ 1D (3.8e-18, "fit to TopBase Opacity Project cs") and 1S (5.5e-18) | absent |
| photoionization into excited states of O^+ (explicit analytic form) | absent |
| Mg^+ excited level ("fit to opacity project data") | absent |
| K rescaled x5, the source noting that Verner and the Opacity Project differ hugely | not applicable |

So **Cloudy photoionizes metals out of excited levels and MoCHII only out of the
ground level.**  As opacity these are small: the excited-level populations are
1e-5 to 1e-4 of the ground, so the O^2+ 1D opacity is of order 1e-4 of the O^+
ground opacity.

An independent Opacity Project comparison was possible for C, N and O through
sirocco's TOPbase tables
(`~/RT_Codes/sirocco/xdata/atomic/topbase_cno_phot_extrap.dat`; the second index
there is the ion stage, verified because Z = 6 stage 2 has E_th = 24.323 eV
against the observed C II 24.383 eV).  Photon-weighted with a 40 kK blackbody
over 24.587-54.4 eV:

| ion | OP / VFKY96 |
|---|---:|
| C II | 0.946 |
| N II | 0.978 |
| O II | 0.863 |

OP is 2-14% **lower**.  With the measured upward derivative that is -0.005 to
-0.03 pc: the wrong direction.  Ne and S are absent from that table and were not
checked.

### The prediction, graded

Recorded before the runs: *"metal opacity x 1.30 moves `r(He0 = 0.5)` outward,
x 0.70 inward, of order 0.1-0.2 pc."*

- **Sign: correct, both ways.**
- **Magnitude: wrong** --- measured +0.03093 and -0.07206 pc, so 1.5 to 6 times
  smaller than predicted.
- **Why it failed:** the prediction extrapolated the on/off bound of 0.615 pc as
  if the response were symmetric, and it is not.  An on/off bound measures what
  is lost by *removal* and says almost nothing about what could be gained by
  *addition*.  **The earlier statement that "a ~30% error in metal opacity would
  be worth ~0.18 pc" is withdrawn**: 30% is worth 0.031-0.072 pc.

## The He ionization network, code by code

`alpha_A(He I)` at 10^4 K, in cm^3 s^-1:

| code | `alpha_A` | how it is obtained | vs Cloudy |
|---|---:|---|---:|
| Cloudy c25.00 | 4.3583e-13 | level-by-level Milne from its own cross sections, total cached | --- |
| MoCHII | 4.3719e-13 | Badnell RR total + DR | +0.3% |
| Wood `ionize` / CMacIonize | 4.2630e-13 | sum of the four BSS99 **effective** coefficients | -2.2% |
| MOCASSIN | **4.8394e-13** | Verner & Ferland 1996 total fit + dielectronic | **+11.0%** |

The MOCASSIN decode (`data/radrec.dat` in the MOCASSIN tree, read order
replicated from its `source/update_mod.f90:1707-1735`) is validated by two
controls: H I 4.1923e-13 against Cloudy's 4.1642e-13 (+0.7%) and He II
2.1880e-12 against 2.1817e-12 (+0.3%).

**No code truncates the recombination cascade to high levels.**  That was the
suspicion this comparison was made to test, and it is not what distinguishes
them.  Three treatments are in play: Cloudy sums explicitly over levels; MoCHII
and MOCASSIN use totals over all `n` (Badnell, VF96); Wood and CMacIonize use
*effective* coefficients, in which cascades from higher levels are already folded
into the `n = 2` terms.  All four therefore include every level.  MoCHII's
Badnell total agrees with Cloudy's explicit level sum to 0.3%, and the 20-level
list in Cloudy's `source/sanity_check.cpp` is a unit test of threshold cross
sections, not the extent of its recombination sum.

**The electron count settles that the spread is about atomic structure, not the
extent of the sum.**

| species | electrons | code-to-code spread in `alpha_A` |
|---|---:|---:|
| H I | 1 (hydrogenic) | 0.7% |
| He II | 1 (hydrogenic) | 0.3% |
| **He I** | **2** | **11-13%** |

The two one-electron species agree to a few tenths of a percent because a
hydrogenic system has an exact analytic solution and there is nothing for two
calculations to disagree about.  He I has no analytic solution: screening and the
singlet/triplet splitting enter, and that is exactly where the four codes part
company.  Cascade extent is common to all of them; the two-electron atomic data
are not.

**Within He I, three cluster and one is the outlier.**  Cloudy's level sum
4.3583e-13, Badnell 4.3719e-13 and BSS99 4.2630e-13 all lie within 2.5% of each
other; Verner & Ferland (1996) alone sits at 4.8394e-13.  MOCASSIN applies
`radRecFit` uniformly with no special case for H or He
(`source/update_mod.f90:1636-1639`, the routine's own comment naming "H-like,
He-like, Li-like, Na-like - Verner & Ferland, 1996, ApJS, 103, 467") and adds
dielectronic recombination only for `elem >= 3`, so its He I rate is VF96's
radiative value alone.  A plausible reason, offered as reasoning and not as a
verified claim: VF96's He-like fits lean on hydrogenic scaling for the excited
levels, and Z = 2 is where that approximation is worst, neutral helium being the
case in which screening and the singlet/triplet splitting matter most.

So **He I is where the codes disagree**: 11-13% across the four, against 0.7% for
H I and 0.3% for He II.  At the measured 0.80 pc per unit that spread is worth
0.09 pc, half the residual --- but MoCHII sits within 0.3% of Cloudy and with the
clustered three, so it is not MoCHII's share of the spread.

**What the spread is, and what it is not.**  The 11% is the **variance between
two data sets**, not a measured error in MoCHII, and for He I the variance is
attributable: Cloudy does *not* use Badnell for the He-like and H-like
sequences but derives them level by level from its own photoionization cross
sections through the Milne relation, which is an independent first-principles
determination.  Against it, Badnell (MoCHII) is +0.3% and VF96 (MOCASSIN) is
+11.0%.  **For He I, therefore, MoCHII is demonstrably right and the spread
belongs to VF96.**  That is a measurement, not a preference among references.

**A caution for MOCASSIN comparisons.**  `CLAUDE.md` lists MOCASSIN's 11-case
regression suite as a 3D validation reference.  Any helium comparison against
MOCASSIN inherits a known **+11% in `alpha_A(He I)` that has been traced to
VF96** by the first-principles determination above, so that difference should be
subtracted before a disagreement in a helium front position or a helium line
intensity is attributed to transport or geometry.

**He I excited levels cannot supply the difference.**  MoCHII's band floor is
`par%eion_min = 13.598 eV`, while every excited He I level has an ionization
potential below it (2^3S 4.767, 2^1S 3.971, 2^3P 3.623, 2^1P 3.369 eV), so
MoCHII structurally cannot photoionize them --- this is a band definition, not a
missing rate.  The size of what is missed: 2^3S dominates the excited population
(2^1S decays by two-photon emission at 51.3 s^-1, leaving it 10^7 times less
populated), with `n(2^3S)/n(He^0)` = 6.2e-5 at r = 3 pc and
`Q(>4.767 eV)/Q(>24.587 eV)` = 33.8 for the 40 kK blackbody, so the extra
ionization is **0.2-0.6%**.  Cloudy's own 2^3S threshold cross section,
5.48e-18 cm^2, confirms the value used in that estimate, and collisional
ionization out of 2^3S gives a comparable ~0.6%.

**The metals are a different case: no arbiter.**  The same comparison was run
across the whole metal registry (radiative recombination only, dielectronic
excluded on both sides because MOCASSIN adds it separately).  Eleven ions: 7
within 10%, 10 within 20%, median 0.971, with the spread concentrated in the
rates that form **neutral atoms** --- C I 0.861, O I 0.837, Ne I 0.898,
Mg II 0.791, against 0.999-1.016 for the rates forming ions.  Here Cloudy
cannot adjudicate: it ships `data/badnell_rr.dat` and `data/badnell_dr.dat`
**md5-identical to MoCHII's copies** and evaluates them with
`RR_Badnell_rate_coef` / `DR_Badnell_rate_coef` in `source/ion_recomb.cpp`, so
Cloudy and MoCHII agree on metal recombination by construction.  With only two
determinations in play, which one is right for C I, O I, Ne I or Mg II is
**undetermined**.  One may reason by analogy from the He I result --- Badnell
(2006) is the later and more careful calculation, VF96 (1996) leans on
simplified treatments of many-electron systems, and He I is where that was
actually demonstrated --- but that is inference, not measurement.  Settling it
needs an independent determination for each ion: R-matrix values for C I, O I
and Ne I recombination, or CHIANTI totals that do not derive from Badnell.

**Wording withdrawn.**  An earlier phrasing in this investigation, *"a 10-20%
uncertainty in MoCHII's neutral-forming recombination rates"*, is **withdrawn**
as imprecise: it reads a spread between two data sets as an error bar on one of
them.  The two accurate statements are

- **He I:** MoCHII agrees with an independent first-principles calculation to
  0.3%, and the VF96 set MOCASSIN uses is 11% high.  Measured.
- **Metal neutrals:** Badnell and VF96 differ by 10-16%, and this benchmark
  configuration has no third source able to say which is right.  A code-to-code
  spread, not a MoCHII error bar.

**The registry-wide audit itself.**  The `alpha_1` defect prompted a sweep of
all 32 registry entries.  It found **no transcription error** (the 22 entries
that take Badnell coefficients reproduce `data/atomic/badnell_{rr,dr}.dat` at
10^4 K to a ratio of 1.0000, with the registry `TRANSITION it` of element Z
mapping to the Badnell row `(Z, N = Z - it, M = 1)`), and it found that the
older `RR2`/`DR2` forms **cannot be migrated**: every one of them sits at a
(Z, N) outside Badnell's coverage (N <= 15, plus N = 18 for Z >= 20), and where
Badnell has nothing the pipeline takes CHIANTI v11.0.2 in exactly the form
CHIANTI declares, so the Shull & Van Steenberg form is what CHIANTI still ships
for those ions rather than a stale choice made in preference to Badnell.  An
earlier suggestion to migrate them is **withdrawn** on both counts.  The full
audit, with the coverage and form-matching tables, is in
`docs/MoCHII_physics.pdf`, "Recombination data across the registry".  Its
verdict is that **no defect was found in MoCHII's recombination data**, which is
not the same statement as "the data are accurate" --- for the metals there is no
arbiter.

## The photon budget: source, and the He-ionizing balance in the nebula

The ionizing budget entering the calculation is verified directly against the
reference.  For the benchmark's 40 kK blackbody (`par%tstar=4.0e4`,
`par%luminosity=1.2734e39` over 13.598-100 eV):

| | MoCHII | Cloudy c25.00 | ratio |
| --- | ---: | ---: | ---: |
| `Q(>13.598 eV)` | 4.2656e49 s^-1 | 4.2653e49 (`Q > 1 Ryd`) | 1.00005 |
| `Q(>24.587 eV)` | 4.6082e48 s^-1 | 4.6081e48 (`Q > 1.8 Ryd`) | 1.00003 |

The two codes are illuminating the same nebula with the same number of photons
to 0.005%, and in particular the He-ionizing budget --- the quantity that would
have to rise by ~13% for a uniform explanation of the front --- agrees to
0.003%.  The source normalization, the spectral shape at the two thresholds, and
the band definition are all excluded.

**The He-ionizing balance inside the nebula also closes**, once the leaf volumes
are formed with the correct convention (see
[Common diagnostic](#common-diagnostic)): the summed gas volume matches the
geometric shell 0.972-4.729 pc to 0.25%, and the budget is

| HII40 | s^-1 |
|---|---:|
| stellar `Q(>24.587 eV)` | 4.6082e48 |
| diffuse He I ground continuum, recycled | 1.6135e48 |
| total He-ionizing supply | 6.2217e48 |
| He^0 photoionizations = case-A recombinations | 4.7008e48 |
| absorbed by H^0 or metals, or escaped | 1.5210e48 = **24.4% of the supply** |

Nothing is missing from the accounting.  What the closed budget shows is *where
the sensitivity lives*: **24.4% of the He-ionizing supply is taken by something
other than He^0**, and halving that loss fraction would move the front 0.211 pc
outward --- more than the whole residual.  That is the quantitative reason the
H^0/He^0 competition, and its energy dependence, is the axis on which a
few-percent difference can compound into 0.18 pc, while a few-percent change in
the He^0 rates themselves cannot.

Dust is excluded by construction rather than by measurement: the Lexington
benchmarks run with `par%dust_model='none'` and `par%ion_add_dust=.false.`, so
there is no grain opacity competing for ionizing photons in either calculation
being compared.

## Next decision

Excluded: diagnostic binning, Monte Carlo noise, spatial resolution, the H/He
recombination rate set, the emitted diffuse spectrum as a whole, the transport of
the diffuse photons, the He I ground-state photoionization cross section, the
He I recombination rate, `alpha_A`/`alpha_1`/`alpha_B`(He I) against Cloudy, the
metal opacity, metal recombination continua as a diffuse source, He I
excited-level photoionization, dust, the source photon budget, the He-ionizing
budget inside the nebula, the H ionization structure, and the resolution of the
reference calculation at the front.  The bounds are tabulated in
[What is excluded, and how tightly](#what-is-excluded-and-how-tightly).

**No single parameter remains.**  Every candidate above has been given its own
paired run or its own external reference, and each one is short of 0.17750 pc by
between one and three orders of magnitude.  The list is a set of measured
levers, not a proof of completeness.

What remains, in the order the evidence supports:

1. **A combination of several few-percent differences, compounding with depth.**
   This is now the leading hypothesis, and it is the one thing the whole
   experimental program so far cannot address: every run varied one quantity.
   The mechanism is already measured --- He^0 carries 84-96% of the opacity above
   24.6 eV and 24.4% of the He-ionizing supply is lost to other absorbers, so a
   difference in the H^0/He^0 competition of a few percent at 1.6 pc is 20% at
   3.8 pc ([The radial signature](#the-radial-signature)).  Testing it requires
   perturbing several quantities together, in the directions their individual
   external uncertainties allow, rather than looking for one more lever.
2. **The He ionization network itself** --- specifically, what Cloudy's resolved
   He I model atom contains that MoCHII's effective treatment does not.  The
   code-by-code comparison narrowed this: it is not the extent of the
   recombination cascade, and it is not the excited-level photoionization, both
   of which are now measured.  What has not been read is Cloudy's He I
   implementation channel by channel.  That is a **separate task**.
3. **Collisional redistribution among the n = 2 terms**, the most concrete
   instance of item 2 that MoCHII could already act on.  The branch fractions are
   effective coefficients in the low-density limit and `p584` is an
   escape-probability closure; `hei_metastable_mod` already carries the
   `q_31a`/`q_31b` rates that would couple the terms at finite density.  It is
   density dependent, so it is unlikely to dominate at the benchmark densities,
   but it is testable without leaving the code.
4. **A Cloudy-side modeling choice not yet identified.**  The metal opacity, the
   remaining candidate outside the He treatment when this list was last written,
   has since been measured and excluded, and the one thing Cloudy demonstrably
   does that MoCHII does not --- photoionizing metals out of excited levels ---
   is of order 1e-4 of the ground-level opacity.

The quantitative statement to keep in view while pursuing any of these: closing
0.17750 pc means moving the front out by 4.1% in radius, hence 12.9% in volume,
which a uniform explanation would have to supply as `Q(He0)` higher by ~13% or
`alpha_B(He I)` lower by a comparable amount.  `Q(He0)` agrees with Cloudy to
0.003%, and although the code-to-code spread in `alpha_A(He I)` does reach 11-13%,
MoCHII sits within 0.3% of Cloudy and with the clustered three of the four
determinations.  **Neither is supported by the data on MoCHII's side.**  The
explanation, if there is one there, has to be progressive in depth rather than a
change of scale.

## The G2a thermal claim: bisected, and it was never a regression

The manuscript carried a median `|dTe|` of 0.14% and a mean bias of +0.006% for
the G2a pure H/He thermal gate. The gate reports 0.97% and +1.02%. The gap was
traced by rebuilding the code at `29b4a46` (2026-07-26), the commit whose
manuscript first carried the 0.14%, in a separate worktree and rerunning the same
case with the period's own gate script:

| code | median `|dev|` | mean bias | max `|dev|` |
|---|---|---|---|---|
| `29b4a46` (2026-07-26) | 0.994% | +1.051% | 3.337% |
| current | 0.969% | +1.016% | 3.108% |

So the gate has reported about one percent since at least 2026-07-26, the
0.14% / +0.006% pair was never produced by it, and the alpha_1 correction made
today **improved** the agreement by 0.025 percentage points rather than degrading
it. There is no regression to find; the manuscript carried a number the gate does
not measure.

What the ~1% is made of, measured rather than assumed:

- **0.19 percentage points** is a data-source mismatch in the gate itself. Its 1D
  reference uses the Hui & Gnedin (1997) case-B recombination coefficient while
  the run uses the Badnell total with the Milne ground-level rate, and the two
  differ by 1.2% at the 19--24 kK of this test. Rerunning the same case with
  `par%recomb_model = 'hui_gnedin'`, matching the reference, gives a median of
  0.776%.
- **The remainder** is the reference's own free-free Gaunt approximation, which
  `check_g2a.py` documents in its header as 1--2% of a cooling term worth about
  10% of the total.

The comparison is therefore a bound on the 1D reference, not a resolution of the
Monte Carlo solution, and the manuscript now says so. Tightening it would mean
giving the reference the same recombination set and a proper Gaunt factor; that
is a gate-side change and is recorded as an open item rather than done here,
because it would move every G2a number ever published.

The same class of mismatch exists in the G1 Stromgren gate, whose 1D reference
also uses Hui & Gnedin while the run uses Badnell (0.86% apart at 10^4 K, which
maps to 0.29% in radius against a reported deviation of +0.60%).

## Aligning the analytic gate references with the code

The G2a finding above -- that a fifth of its deviation was the gate's reference
using a different recombination set -- turned out to be one of three such
mismatches. The gates and the paper's figure generators each carried their own
copy of the code's physics, and the copies had drifted:

| axis | reference carried | the run uses | size |
|---|---|---|---|
| recombination | Hui & Gnedin (1997) case B | Badnell total minus the Milne ground-level rate | 0.86% at 10^4 K, 1.2% at 19-24 kK |
| free-free Gaunt | frequency-averaged fit, no charge dependence | `gbar_ff(T, Z)` with Z = 1 and Z = 2 separated | 6.5% between Z = 1 and Z = 2 at 10^4 K |
| photoionization | Verner et al. (1996) fits for H I and He II | the exact hydrogenic expression | 0.775% at the H I threshold |

All three are gate-side. The Fortran is correct in each case --
`cooling_mod.f90:257` separates the charge states exactly as its own header
states, and `photo_xsec.f90` uses the exact hydrogenic form for both hydrogenic
species.

`tools/python/mochii_rates.py` now holds all three, taken from the code's own
routines, and the seven files that had copies import from it.
`tests/rates/check_python_rates.py` compiles a driver against the production
objects and checks the module against the Fortran: recombination to 8.9e-16,
`gbar_ff` to 2.2e-16, cross sections to 1.3e-15.

### What each alignment was worth

| stage | G0 `Gamma_HI` bias | G1 level 6 | G2a median | dusty full / l09 |
|---|---|---|---|---|
| before | +0.004% | +0.60% | 0.969% | +0.31% / +0.57% |
| + recombination | --- | +1.00% | 0.765% | +0.66% / +0.98% |
| + Gaunt charge dependence | --- | --- | 0.762% | --- |
| + cross sections | +0.077% | +1.0011% | 0.773% | +0.66% / +0.98% |

**Only the recombination mattered.** The Gaunt charge dependence moved the G2a
median by 0.003 percentage points, because a 40 kK blackbody leaves almost no
photons above 54.4 eV and the Z = 2 term is therefore negligible; the cross
sections moved it 0.011, and in the wrong direction. The earlier attribution of
the G2a remainder to the reference's Gaunt approximation is **withdrawn**: that
approximation is accurate to 0.06-0.6% for Z = 1, measured against the code's
own tabulated values.

### The G1 criterion moved to level 7

Aligning the reference pulled its R_eff from 3.0365 to 3.0243 pc while the
level-6 Monte Carlo value stayed at 3.0546 pc, so the level-6 deviation went
from +0.60% to +1.0011% and crossed the gate's documented 1% criterion. The
former margin was two errors cancelling.

The +1.00% is not an error but the finite front width on a 0.147 pc cell against
a front at 3.02 pc, about a fifth of a cell. A uniform level-7 run
(`tests/g1_stromgren/g1_stromgren_l7.in`, 2 097 152 leaves) gives
R_eff = 3.0448 pc, **+0.6789%**, and the criterion now sits there with level 6
reported alongside as the resolution it is. The convergence from +1.00% to
+0.68% on halving the cell is what identifies the effect as discretization.

The gate also only printed its result; it now computes a verdict and exits
nonzero on failure, like the others.

The G0 case is the clean counter-example to G1. Its rate integral is exactly
what a cross-section difference biases, and there the alignment showed a bias of
+0.077% where the mismatched reference had reported +0.004% -- the old figure was
the reference's 0.073% overestimate cancelling against the true bias. Nothing
about the run changed; only what it was compared against.

### Completing the alignment, and what it revealed about RQMC

Four more files carried copies: `tests/qmc/check_v3.py`, `paper/python/make_qmc.py`,
`tests/peel/check_peel.py` and `tests/peel_ext/check_peel_ext.py`. Three used the
Verner et al. (1996) fit for H I; `peel_ext` already had its own exact hydrogenic
form and its output is bit-for-bit unchanged, which is the control that the
module reproduces what a correct copy did.

The two transport tests were expected to be insensitive. They were not.

| | before | after |
|---|---|---|
| `peel` analytic-attenuation ratio - 1 | +0.1507% | **+0.0195%** |
| `qmc` RQMC mean \|bias\|, N = 2^20 | 0.062% | **0.011%** |
| `qmc` radial-profile RQMC rms, N = 2^20 | 0.119% | **0.032%** |
| `qmc` radial-profile gain, N = 2^20 | 4.12x | **30.17x** |
| `qmc` cell-wise gain, N = 2^20 | 222.5x | **231.1x** |

The peel deviation was almost entirely the fit's 0.775% error at the threshold,
where the optical depth is largest -- an order of magnitude of apparent transport
error that was the reference's.

The RQMC case is the more interesting one. The scrambled-Sobol launch had been
reported as carrying a residual bias of about 0.06%, and that bias was the
reference, not the sampler: with the cross sections aligned it is 0.011%. The
radial-profile variance-reduction gain had been floored by it at 4.1x and is
30.2x once it is gone. **The variance reduction was understated because the thing
it was measured against was wrong.** The manuscript's gains are corrected to
28.8, 66.7, 143 and 231 times at N = 2^17 to 2^20.

`sig_HI_pow` in `peel_ext` was kept and renamed `hydrogenic_power_law`: it is a
deliberate foil, showing that the measured optical-depth ratio follows the exact
hydrogenic shape to 0.00% and a textbook E^-3 power law only to 9.65%. A
reference the code is *not* meant to match has to be labelled as such, or the
next reader aligns it too.

### The recombination-cooling and collisional-ionization coefficients

`betaA_HII`, `betaB_HII`, `beta_HeII` and `beta_HeIII` are now exported from
`cooling_mod`, and the Voronov collisional-ionization rates were already public,
so all seven live in the shared module and
`tests/rates/check_python_rates.py` holds them to the Fortran along with
everything else -- 23 quantities, all to machine precision. Those expressions
agreed before the move, so the G2a numbers are unchanged by it; what changes is
that they can no longer drift apart silently, which is how the three mismatches
above reached the published gate numbers in the first place.
