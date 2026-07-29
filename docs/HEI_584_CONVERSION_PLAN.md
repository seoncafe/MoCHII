# Plan: the fate of the He I 584 A (2^1P) photon in the diffuse field

- Date: 2026-07-30
- Target: channel 4 of `src/diffuse_mod.f90` (`par%hei_diffuse`), and
  `gamma_2q_HeI` in `src/nebcont_mod.f90`
- Status: **implemented, gated and measured (2026-07-30)**.
  `hei_584_ionizes_H` in `src/hei_cascade_mod.f90`, folded into
  `hei_branch_energies` in `src/diffuse_mod.f90` and into `gamma_2q_HeI` in
  `src/nebcont_mod.f90`; gate `tests/hei_cascade/check_hei_584.f90`.
  Measured effect and the prediction grading are in section 7 below.
- Follows: [HEI_CASCADE_DIFFUSE_PLAN.md](HEI_CASCADE_DIFFUSE_PLAN.md), whose
  step 1 (the Drake two-photon spectrum) and step 2 (temperature-dependent
  branch fractions) are implemented and gated. This is step 3, and it is the
  larger of the three by an order of magnitude.

## 1. The defect

Roughly one sixth of every case-B He II recombination reaches the ground
state through 2^1P, radiating He I Ly-alpha at 584 A = 21.22 eV. MoCHII emits
that photon into transport as a 21.22 eV H-ionizing packet and lets the grid
decide where it is absorbed.

That is one of two limits, and it is the wrong one. The 584 A photon is a
**resonance** photon: the medium it is emitted into is optically thick to it
in He^0, so it scatters many times before it does anything irreversible. Two
irreversible outcomes compete during that random walk:

1. absorption by a hydrogen atom, which photoionizes it (21.22 > 13.60 eV) and
   deposits 7.62 eV of heat;
2. collisional transfer 2^1P -> 2^1S in the He atom that last absorbed it,
   after which the decay is the two-photon continuum instead.

Which one wins is set by the competition between the H^0 and He^0 opacities
seen by a trapped 584 A photon, so it is a **local, density- and
ionization-dependent branching**, not a constant. MoCHII currently assumes
outcome 1 with probability one.

`nebcont_mod` independently assumes the opposite limit --- it uses the total
singlet recombination coefficient `6.23e-14 (T/1e4)^-0.827` for the He I
two-photon nebular continuum, which is correct only if *every* 2^1P decay
converts. So the two modules disagree about the fate of the same decay, and
each is right only at an opposite extreme.

## 2. What the reference codes do

This was settled by reading the codes, not by argument. All three of the
Monte Carlo photoionization codes in the reference trees treat 2^1P the same
way, and MoCHII is the outlier.

### K. Wood, `ionize` (`~/RT_Codes/Wood_Codes/ionize/`)

`probset.f` sets four effective recombination coefficients from BSS99 Table 1
+ J. S. Mathis and normalizes by their sum:

```fortran
      alpha_1_He =1.54e-13*(Telec/1.e4)**(-0.486)
      alpha_e_2tS=2.1e-13 *(Telec/1.e4)**(-0.381)
      alpha_e_2sS=2.06e-14*(Telec/1.e4)**(-0.451)
      alpha_e_2sP=4.17e-14*(Telec/1.e4)**(-0.695)
c**** Set alphaHe=sum(alpha).  Do this because of different T-dependences
c**** in fitting formulae
      alphaHe=alpha_1_He+alpha_e_2tS+alpha_e_2sS+alpha_e_2sP
```

`mcionize.f` then handles the 2^1P branch **without transporting it**:

```fortran
c***************** HeI Ly-alpha.  Do not propagate this with Monte Carlo.
c***************** Test if it is absorbed by H OTS or converted to HeI
c***************** two-photon continuum.
                  pHots=1./(1.+0.77/(sqrt(temp(...))/1.e2)
     +                              *nfracHe(...)/nfracH(...))
                  if(ran2(iseed).lt.pHots) then  !absorbed by H OTS
                     ... H Lyman continuum or terminate ...
                   else
                     if(ran2(iseed).lt.0.56) then
c****************** prob of ionizing photon from He 2-photon continuum is 0.56
                       pnu=nusetHe2q(pHe2q,nu,iseed)
```

The same `pHots` enters the ionization balance (`h0find.f`) and the heating
(`ioneng.f`, where the converted-away fraction deposits
`1.2196e-11` erg = 7.62 eV = 21.22 - 13.60).

### CMacIonize (`~/RT_Codes/CMacIonize/CMacIonize-master/src/`)

`PhysicalDiffuseReemissionHandler` is a line-for-line reimplementation of the
above, and cites its source: *"Wood, Mathis & Ercolano (2004), section 3.3"*.
Same four coefficients, same sum normalization, same cumulative channel
structure (`HELIUM_LYC`, `HELIUM_NPEEV`, `HELIUM_TPC`, `HELIUM_LYA`), and the
same 584 A branching:

```cpp
    // helium Lyman alpha, is either absorbed on the spot, or converted to
    // helium two-photon emission
    const double pHots = sqrtTnH0 / (sqrtTnH0 + 77.*ionic_fraction(ION_He_n));
    if (x < pHots) { /* absorbed on the spot by hydrogen */ }
    else           { /* helium two-photon continuum */ }
```

Its form is algebraically identical to Wood's: dividing through by
`sqrt(T) x(H^0)` gives `1/(1 + 77 x(He^0)/(sqrt(T) x(H^0)))`, which is Wood's
`1/(1 + 0.77 (100/sqrt(T)) x(He^0)/x(H^0))`.

CMacIonize is also where MoCHII's constant `0.56` came from, with its
derivation in the comment: *"two photons are emitted; a fraction of 28% of
these photons has an energy high enough to ionize hydrogen (and since there
are two, the probability is 56%)"*. MoCHII step 1 replaced the constant by
0.55644 integrated from the Drake table.

### MOCASSIN (`~/RT_Codes/MOCASSIN/.../source/emission_mod.f90`)

MOCASSIN does not resolve the branching at all; it takes the opposite limit
and says so:

```fortran
        ! assume all HeI singlets finally end up in the 2^1S
        ! use total recombination cofficient to all singlets
        alphaEff21SHeI = 6.23*((TeUsed/10000.)**(-0.827))
```

`6.23e-14` is the singlet total: `2.06e-14 + 4.17e-14 = 6.23e-14` exactly, the
2^1S and 2^1P coefficients of BSS99 added together. (The exponent `-0.827` is
a separate fit to the sum, not the sum of the two exponents.) This is the
`pHots -> 0` limit made into a constant.

MoCHII inherited this expression for `gamma_2q_HeI` when the nebular continuum
was ported, which is why the two MoCHII modules sit at opposite limits.

### Summary

| code | 2^1S branch | 2^1P (584 A) fate | transported? |
|---|---|---|---|
| Wood `ionize` | `2.06e-14 T4^-0.451` | `pHots` -> ionizes H; `1-pHots` -> 2-photon | no, on the spot |
| CMacIonize | `2.06e-14 T4^-0.451` | identical to Wood | no, on the spot |
| MOCASSIN (continuum) | `6.23e-14 T4^-0.827` = singlet total | all converted to 2^1S | n/a, emissivity only |
| **MoCHII now** | `2.06e-14` ratio (step 2) | **all ionize H** | yes, full transport |
| **MoCHII after step 3** | same | `pHots` split, then transported | yes, full transport |

MoCHII differs from Wood and CMacIonize on one axis in its own favor: those
codes destroy the 584 A photon on the spot, while MoCHII transports whatever
comes out of the branch through the grid, so the *deposition* is resolved
rather than assumed local. That advantage is worth nothing while the branch
content itself is wrong, which is what this step fixes.

## 3. How large this is

`pHots` evaluated on the converged production grids, using the same formula
the reference codes use, weighted by the 2^1P emission rate `n_e n(He^+)`:

| benchmark | median `x(He^0)/x(H^0)` | emission-weighted `<pHots>` | MoCHII assumes |
|---|---|---|---|
| HII40 | 8.0 | **0.117** | 1.0 |
| HII20 | 282 | **0.0044** | 1.0 |

So 88% of the 2^1P branch in HII40, and 99.6% in HII20, should convert to the
two-photon continuum and does not. HII20 is the more extreme case because its
softer 20 kK spectrum leaves helium far more neutral relative to hydrogen.

Consequences for what channel 4 emits per case-B He II recombination, using
the step-2 branch fractions and the step-1 two-photon moments:

| | HII40 (`pHots` = 0.117) | HII20 (`pHots` = 0.0044) |
|---|---|---|
| energy now | 19.22 eV | 19.24 eV |
| energy after | 17.47 eV (**-9.1%**) | 17.12 eV (**-11.0%**) |
| photons now | 0.966 | 0.966 |
| photons after | 0.903 (**-6.6%**) | 0.890 (**-7.9%**) |

and the spectrum softens: a 21.22 eV line is replaced by a continuum spread
over 13.6-20.6 eV with mean 16.1 eV.

**Honest scaling against the open residual.** Enabling channel 4 at all moves
`r(He^0 = 0.5)` outward by 0.103 pc. Removing ~10% of its content scales
naively to ~0.01 pc. That is twenty times larger than every correction
measured so far (Milne ground continua +0.00063 pc, Drake two-photon
-0.00030 pc, step-2 branch fractions ~+0.0004 expected) but still about ten
times smaller than the 0.18 pc Cloudy residual. This step should be expected
to be the largest remaining item inside the diffuse field, and should **not**
be expected to close the residual on its own. The spectral softening acts on
top of the photon-count change and is not captured by the linear scaling, so
the measured value may exceed the estimate; that is a reason to measure, not
a reason to predict a larger number now.

## 4. Implementation

### 4.1 The branching probability

Add to `hei_cascade_mod`, since it is a property of the cascade and belongs
with the branch fractions:

```fortran
  !--- Escape-probability branching of the trapped He I 584 A resonance
  !--- photon: absorbed by hydrogen (photoionizing it) versus collisionally
  !--- converted 2^1P -> 2^1S and re-emitted as the two-photon continuum.
  !--- Wood, Mathis & Ercolano (2004) sect. 3.3; the same expression is in
  !--- K. Wood's ionize (mcionize.f, h0find.f) and CMacIonize.
  real(kind=wp) function hei_584_ionizes_H(T, xHI, xHeI) result(p)
     p = 1.0_wp/(1.0_wp + 77.0_wp/sqrt(T)*xHeI/max(xHI, tinest))
  end function
```

with `xHI -> 0` giving `p -> 0` (fully converted) and `xHeI -> 0` giving
`p -> 1` (nothing left to convert), both correct limits.

### 4.2 The branch energies

No new Sobol dimension is needed, and none may be added: the dimension
assignment is fixed (7 = channel, 8 = branch, 9 = two-photon energy) and
changing it would invalidate every stored RQMC comparison. The split folds
into the branch weights at build time, because `pHots` depends only on the
leaf's own state:

```fortran
  subroutine hei_branch_energies(T_K, xHI, xHeI, e3s, e1p, e1s)
    call hei_branch_fractions(T_K, f3s, f1s, f1p)
    p   = hei_584_ionizes_H(T_K, xHI, xHeI)
    e3s = f3s*HEI_E3S
    e1p = f1p*p*HEI_E1P                                  ! survives as 584 A
    e1s = (f1s + f1p*(1.0_wp - p))*hei_2ph_energy_per_decay()
  end subroutine
```

`diffuse_build` already has `gas_xHI(il)` and `gas_xHeI(il)` in hand, and both
`gen_diffuse_photon` and `gen_diffuse_photon_qmc` already call
`hei_branch_energies` with the leaf index available, so the change is
confined to one signature and its three call sites.

**Energy bookkeeping.** A converted decay radiates 20.62 eV as the two-photon
pair, not 21.22 eV. The difference, 0.60 eV, leaves as the 2^1P -> 2^1S
transition at 2.06 micron, which is non-ionizing and correctly absent from the
band. This must be stated at the code site, since dropping 0.6 eV without
explanation looks like a leak.

### 4.3 The nebular continuum

Once step 3 lands, `gamma_2q_HeI` must use the same split:

```
alpha_2q(T, xHI, xHeI) = alpha(2^1S) + [1 - pHots] alpha(2^1P)
```

which reduces to MOCASSIN's `6.23e-14` exactly in the `pHots -> 0` limit and
to `2.06e-14` in the other. Only then do the transported packets and the
written continuum describe the same nebula. Doing this **before** step 3 would
make the two modules consistent at the wrong limit, which is why the earlier
proposal to simply replace `6.23e-14` by `2.06e-14` is withdrawn.

The exponent question resolves itself: with the split evaluated from the two
BSS99 coefficients, the fitted `-0.827` sum is no longer used at all.

## 5. Validation

1. **Unit** (`tests/hei_cascade/check_hei_584.f90`, new): `pHots` reproduces
   Wood's and CMacIonize's algebraically distinct forms to round-off on a grid
   of `(T, xHI, xHeI)`; the limits `xHI -> 0` and `xHeI -> 0` are exact; the
   branch energies still sum to a physically bounded per-recombination energy
   and the photon count is monotonically decreasing in `1 - pHots`.
2. **Regression**: with `pHots` forced to 1, channel 4 must reproduce the
   step-2 output bit-for-bit. This is the check that the refactor of
   `hei_branch_energies` changed nothing but the branching.
3. **Paired production runs**, HII40 and HII20 at 68 ranks with the same Sobol
   seed, measured against the 0.002% seed-to-seed noise floor:
   `r(He^0=0.5)`, `r(C^2+=0.5)`, `r(H^0=0.1)`, `<T[NpNe]>`, `L(Hbeta)` and the
   Lexington line ratios.
4. **Prediction to record before running**: the diffuse field loses 7-8% of
   its channel-4 photons and softens, so the He I front should move **inward**,
   by of order 0.01 pc, more in HII20 than in HII40. Recording the sign in
   advance is required: the Milne step's sign prediction was wrong, and the
   only reason that was caught was that it had been written down first.
5. **Nebular continuum**: the He I two-photon component of `_nebcont.txt`
   rises where `pHots` is small, and the total must remain finite and smooth
   across the He^+ front.

## 6. Out of scope

- Explicit He I 584 A resonance-line transfer. `pHots` is an
  escape-probability closure, and this step adopts it as such. Its domain of
  validity --- an optically thick, static, dust-free 584 A line --- must be
  documented at the code site. Dust competes for the trapped photon and is not
  in this closure; at the benchmark dust content that is a small correction,
  but at higher dust-to-gas it is not, and the code comment must say so.
- The 2^3S branch, which stays a 19.82 eV packet in transport as in every
  reference code.
- Collisional redistribution among the `n = 2` terms at high density, still
  deferred from [HEI_CASCADE_DIFFUSE_PLAN.md](HEI_CASCADE_DIFFUSE_PLAN.md) 8.

## 7. Outcome

### 7.1 What was measured

Paired 68-rank runs at the same Sobol seed, against step 2, front positions
from the 0.01-pc volume-weighted shell diagnostic:

| | HII40 | HII20 |
|---|---:|---:|
| dr(He0 = 0.5) | **-0.00363 pc** | -0.00022 pc |
| dr(C2+ = 0.5) | -0.00332 pc | -0.00001 pc |
| dr(H0 = 0.1) | **-0.01002 pc** | (H is density bounded) |
| dL(Hbeta) | -0.1997% | -0.0704% |
| mean Te(x_HI < 0.95) | 9080.2 -> 9053.4 K | 8523.0 -> 8522.3 K |
| diffuse band luminosity | -1.56% | -0.22% |

All three cascade gates pass, and
`tests/continuous_energy/check_hii_production.py` passes on both benchmarks.
The gate confirms this module's expression agrees with Wood's and
CMacIonize's algebraically distinct forms to 2.2e-16.

This is the largest single correction of the four made to what the diffuse
field emits --- six times the Milne ground continua and twelve times the
step-2 branch fractions.

### 7.2 The predictions, graded

Section 5.4 recorded three predictions before the runs. Two held and one
failed, and the failure is the informative one.

- **Sign: correct.** The fronts moved inward, as the loss of channel-4
  photons requires.
- **Magnitude of order 0.01 pc: correct.** The HII40 H0 front moved
  -0.01002 pc; the He0 front moved -0.00363 pc.
- **"More in HII20 than in HII40": WRONG.** HII40 moved sixteen times
  further. The prediction was read off `p_Hots` alone, which measures only
  what fraction of the 2^1P branch converts --- 99.6% in HII20 against 88% in
  HII40. Front displacement instead scales with the branch's *absolute*
  contribution to the ionizing budget. Inverting the measured luminosity
  drops against the per-recombination change gives channel 4 as **17.6% of
  the HII40 diffuse field and 2.0% of the HII20 one**: the 20 kK spectrum
  makes few photons above 24.6 eV, so He+ is scarce and the channel that
  feeds on its recombinations is small. A larger converted fraction of a much
  smaller channel moves less.

  The lesson generalizes: a local branching probability predicts what happens
  per event, never how far a front moves. The second factor is the channel's
  share of the budget, and it is set by the source hardness.

### 7.3 What this says about the Cloudy residual

r(He0 = 0.5) against the Cloudy c25 reference at 4.30531 pc:

| | r(He0=0.5) | residual |
|---|---:|---:|
| step 2 | 4.12333 | -0.18199 |
| step 3 | 4.11970 | **-0.18561** |

**Step 3 moves MoCHII away from Cloudy.** Summed over all four corrections to
the emitted diffuse spectrum --- Milne ground continua +0.00063, Drake
two-photon -0.00030, Wood branch fractions +0.00024, 584 A conversion
-0.00363 --- the net is **-0.00306 pc**, i.e. getting the emission physics
right makes the agreement marginally worse.

That is not an argument against any of the four. Each is independently
correct, each matches an external reference implementation or an exact
relation, and `docs/../CLAUDE.md`'s acceptance criterion is physical
correctness, not agreement with a particular reference calculation. What the
result does settle is the diagnosis: **the 0.18 pc residual is not in what
the diffuse field emits.** Four independent corrections to the emitted
spectrum, spanning its shape, its branch weights and the fate of its
strongest line, together amount to 1.7% of the gap and in the wrong
direction.

What remains, in order of remaining plausibility:

1. **The transport of the diffuse photons**, as opposed to their emission ---
   the one axis never varied in this investigation.
2. **He I ground-continuum atomic data.** `milne_setup` reports alpha_1 from
   the Milne integral disagreeing with the fitted alpha_A - alpha_B by
   2--3.4% for He I, against 0.5% for H I and 0.4% for He II, in exactly the
   continuum whose front is discrepant.
3. **Collisional redistribution among the n = 2 terms**, still out of scope
   here and density dependent.
4. **A difference outside the diffuse field entirely** --- the metal opacity
   at the front, or a Cloudy-side modeling choice not yet identified.
