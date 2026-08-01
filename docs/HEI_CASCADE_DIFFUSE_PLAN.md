# Plan: level-resolved He I recombination cascade in the diffuse field

- Date: 2026-07-29
- Target: channel 4 of `src/diffuse_mod.f90` (`par%hei_diffuse`)
- Status: **both steps implemented, gated and measured (2026-07-30)**, and so
  is step 3, so the physics of channel 4 is now closed end to end: the
  branch weights, the two-photon shape, and the fate of the 584 A photon all
  come from one consistent set of data, and `nebcont_mod` and `diffuse_mod`
  no longer disagree about any of them.
  Step 1 `src/hei_twophoton_mod.f90`, step 2 `src/hei_cascade_mod.f90`,
  step 3 `hei_584_ionizes_H` in the same module folded into
  `hei_branch_energies` and `gamma_2q_HeI`, gate `tests/hei_cascade/`
  (three checks).  Measured effect and the conclusion are in
  [HII40_CLOUDY_RESIDUAL_INVESTIGATION.md](HII40_CLOUDY_RESIDUAL_INVESTIGATION.md):
  the two steps here move `r(He0=0.5)` by -0.00030 and +0.00024 pc and step 3
  by -0.00363 pc, so like the Milne ground continua all three are correctness
  fixes and not the explanation of the Cloudy residual --- the four corrections
  to the emitted diffuse spectrum sum to -0.00306 pc against a 0.18 pc gap,
  which is what exhausted that axis.

  Section 5 is followed as written, including its prohibition on deriving
  `alpha_eff(2^1P)` from closure: `hei_cascade_mod` carries an independent
  2^1P coefficient and the three branches are normalized by their own sum.

  **What was implemented first, and why it was wrong.**  A first version did
  derive 2^1P from closure, reading `nebcont_mod`'s
  `6.23e-14 (T/1e4)^-0.827` as `alpha_eff(2^1S)`.  That gave branch fractions
  0.746 : 0.224 : 0.031 at 8000 K and an argument that case-B trapping of the
  584 A line inverts the singlet split.  The premise was false:
  `2.06e-14 + 4.17e-14 = 6.23e-14` exactly, so `6.23e-14` is the singlet
  *total*, which is what MOCASSIN says at the line it comes from
  (`source/emission_mod.f90`: "assume all HeI singlets finally end up in the
  2^1S / use total recombination cofficient to all singlets").  Section 3
  below propagated the same misreading in its table and is corrected there.
  The independent 2^1P coefficient the plan asked for did exist in the
  reference trees all along, in K. Wood's `ionize`.  The closure derivation
  and the inversion argument are withdrawn; the physical observation behind
  the argument is real and became step 3
  ([HEI_584_CONVERSION_PLAN.md](HEI_584_CONVERSION_PLAN.md)).
- Follows: the Milne ground-continuum sampler, which
  [MILNE_DIFFUSE_LYC_IMPLEMENTATION_PLAN.md](MILNE_DIFFUSE_LYC_IMPLEMENTATION_PLAN.md)
  deferred this work from
- Followed by: [HEI_584_CONVERSION_PLAN.md](HEI_584_CONVERSION_PLAN.md), step
  3 --- the local split of the 2^1P branch between hydrogen absorption and
  collisional conversion to 2^1S, now implemented, gated and measured.
  Section 8 of this plan treated the 584 A photon as ionizing hydrogen on
  emission; that is one of two limits and it is the wrong one, and the step
  turned out to be twelve times larger than step 2 here.

## 1. Why this is the next step

[HII40_CLOUDY_RESIDUAL_INVESTIGATION.md](HII40_CLOUDY_RESIDUAL_INVESTIGATION.md)
excludes diagnostic binning, Monte Carlo noise, spatial resolution, the H/He
recombination rate set, and now the ground-continuum spectral shape. The He I
excited cascade is what the remaining evidence points at:

- **It is the largest measured single effect.** Enabling channel 4 moves
  `r(He0 = 0.5)` outward by 0.103 pc, the same order as the unexplained
  0.18 pc residual. The Milne ground-continuum correction moved it by
  0.0006 pc.
- **It is the crudest physics left in the diffuse field.** The whole channel
  is four constants, none of which depends on temperature.
- **The residual is in the He I front**, and the `alpha_1` cross-check added
  by `milne_setup` shows He I disagreeing by 2-3.4% between the Milne
  integral and the fitted `alpha_A - alpha_B`, while H I and He II agree to
  0.5%. That points at the He I data specifically.

## 2. What the code does now

`src/diffuse_mod.f90` represents the entire case-B He I cascade with

```fortran
HEI_F3S = 0.75                       ! 2^3S branch
HEI_F1P = 0.25*(2/3), HEI_E1P = 21.22 ! 2^1P, 584 A
HEI_F1S = 0.25*(1/3), HEI_P2PH = 0.56 ! 2^1S two-photon
HEI_E1S = 0.5*(13.598 + 20.62)       ! mean of a FLAT distribution
```

Three approximations are stacked here:

1. **Fixed statistical branching.** The 3/4 : 1/4 triplet-singlet split and
   the 2/3 : 1/3 split within the singlets are the low-density statistical
   weights of \S4.6 of Osterbrock & Ferland (2006). They carry no temperature
   dependence and no level resolution.
2. **A flat two-photon spectrum.** The 2^1S continuum is sampled uniformly
   over [13.598, 20.62] eV. The real 2^1S -> 1^1S distribution is strongly
   peaked at half the transition energy.
3. **A fixed H-ionizing fraction**, `HEI_P2PH = 0.56`, which should instead
   follow from the spectral shape above the H I edge.

## 3. The data are already in the repository

This is the decisive point: the same physics is represented twice in MoCHII,
correctly in one place and by a surrogate in the other.

| quantity | already present | used by | not used by |
|---|---|---|---|
| He I 2^1S two-photon shape `A(y)` | `data/HeI2phot.dat`, Drake, Victor & Dalgarno (1969), 41 points symmetric about `y = 0.5` | `nebcont_mod` (nebular continuum) | `diffuse_mod` (flat surrogate) |
| `alpha_eff(2^1S + 2^1P, T)` --- the singlet **total**, not the 2^1S branch | `6.23\times10^{-14}(T/10^4)^{-0.827}` | `nebcont_mod`, where it assumes every 2^1P decay converts to 2^1S | `diffuse_mod` (fixed 1/12 for 2^1S) |
| `alpha_eff(2^3S, T)` | `hm_alpha3 = 2.10\times10^{-13}(T/10^4)^{-0.778}`, overridable from `data/atomic/hei_metastable.txt` | `hei_metastable_mod` | `diffuse_mod` (fixed 3/4) |
| He I case-B line emissivities | `data/atomic/hei_porter_caseB.txt`, Porter et al. (2012, 2013), density-dependent | the line output | `diffuse_mod` |

So the nebular continuum MoCHII writes and the diffuse packets MoCHII
transports currently disagree about the shape of the same 2^1S two-photon
emission. Removing that inconsistency is most of this task.

The one coefficient not already present is `alpha_eff(2^1P, T)`, the 584 A
resonance channel. Section 5 covers how to obtain it.

*Correction (2026-07-30).* This table originally listed the `6.23e-14` fit as
`alpha_eff(2^1S, T)`. It is the singlet total, and reading it as the 2^1S
branch is what produced the withdrawn closure derivation recorded in the
status block. The two modules also disagree about more than the two-photon
*shape*: `nebcont_mod`'s use of the singlet total assumes every 2^1P decay
converts to 2^1S, while `diffuse_mod` assumes none does. That is step 3.

## 4. What must be conserved

Two invariants have to survive the change, and both are testable.

**Photon rate.** Channel 4 emits one packet per case-B He II recombination
that reaches an H-ionizing decay. With equal-energy packets counted from an
energy luminosity, the sampled energies must follow `E n(E)`, not `n(E)`, and
the luminosity must use the number-weighted mean --- the same pairing the
Milne module establishes and `tests/milne` checks. Reuse that reasoning; do
not re-derive it per branch.

**Total case-B rate.** The sum of the level-resolved `alpha_eff` over the
branches must reproduce `alpha_B(He II)` to within the atomic data's own
consistency, and any shortfall must be reported rather than renormalized
away. If the branch coefficients sum to something other than `alpha_B`, that
is an atomic-data statement worth printing, exactly as `milne_setup` prints
the `alpha_1` comparison.

## 5. Implementation, in two separately measurable steps

The two changes must land separately. Together they are expected to move the
He I front by an amount comparable to the residual, and a single combined
commit would make the source of that motion unattributable.

### Step 1 --- the physical two-photon spectrum

Replace the flat draw with the `A(y)` table already in `data/HeI2phot.dat`.

- Add the shape to a small module, or extend `milne_recomb_spectrum_mod`,
  building the CDF of the **energy-luminosity** density `E n(E)` restricted
  to `E >= E_th(H I)`, since only that part ionizes hydrogen.
- Derive `HEI_P2PH` from the same table instead of carrying 0.56 as a
  constant: it is the fraction of the two-photon *photon* emission above
  13.598 eV, and the mean energy that `diffuse_build` uses must be the mean
  of the same truncated distribution.
- `A(y)` is symmetric about `y = 0.5` and each decay emits two photons whose
  energies sum to 20.62 eV; the sampler must draw one photon, not the pair,
  and must not double-count the normalization.
- Keep the Sobol dimension assignment unchanged: dimension 9 remains the
  two-photon energy variate, dimension 8 the branch selector.

Expected sign: the true distribution puts more weight near 10.3 eV, below
the H I edge, than a flat distribution does, so the H-ionizing fraction
should **fall** relative to 0.56 and this step should move the He I front
inward. That prediction should be recorded before the run, and checked.

### Step 2 --- temperature-dependent branch weights

Replace `HEI_F3S`, `HEI_F1P`, `HEI_F1S` with `alpha_eff(T)` per branch.

- The three coefficients must come from **one** source. A fraction assembled
  from fits with unrelated temperature dependences is not a fraction of
  anything, and this requirement decides every choice below.
- Report the branch sum against `alpha_B(He II)` at setup.

*As implemented (2026-07-30).* `src/hei_cascade_mod.f90` takes all three from
Benjamin, Skillman & Smits (1999) Table 1 with J. S. Mathis, in the form
tabulated in K. Wood's `ionize` (`probset.f`) and used with the same sum
normalization by CMacIonize, which cites Wood, Mathis & Ercolano (2004)
section 3.3:

    alpha(2^3S) = 2.10e-13 (T/1e4)^-0.381
    alpha(2^1S) = 2.06e-14 (T/1e4)^-0.451
    alpha(2^1P) = 4.17e-14 (T/1e4)^-0.695     [cm^3 s^-1]

giving 0.762 : 0.076 : 0.162 at 8000 K against the statistical
0.750 : 0.083 : 0.167. Two instructions above were changed on contact with
that requirement:

- **2^3S is deliberately not shared with `hei_metastable_mod`.** Its
  Oklopcic & Hirata (2018) fit `2.10e-13 (T/1e4)^-0.778` is the same
  prefactor with a steeper exponent, 9% higher at 8000 K. It is an absolute
  production rate for the metastable population; this one is a member of a
  set whose ratios are the branch fractions. Sharing it would have removed
  exactly the internal consistency the fractions depend on, so the two fits
  serve their own purposes and the module says why at the code site.
- **The 2^1P coefficient did exist in the reference trees**, so nothing had
  to be fitted and closure stayed forbidden as section 4 requires.

The fractions are normalized by the sum of the three coefficients rather than
against `alpha_B(He II)`, which is Wood's own choice, made for the reason his
source states --- "different T-dependences in fitting formulae". The sum is
0.767, 0.896, 0.970 and 1.134 of MoCHII's `alpha_B(He II)` at 5000, 8000,
10^4 and 1.5e4 K; `hei_cascade_setup` prints that comparison instead of
absorbing it, exactly as section 4 asks. The fits are nebular and are clamped
outside 3e3-2e4 K rather than extrapolated.

## 6. Validation

1. **Unit, before any transport.** The two-photon CDF integrates to one,
   reproduces the tabulated `A(y)` quantiles, and its mean over
   `E >= 13.598` eV matches quadrature. The pseudorandom and Sobol paths
   return identical energies for identical uniforms.
2. **Photon-rate identity**, as in `tests/milne`: the represented photon rate
   equals the emitting recombination rate.
3. **Branch closure**: the level-resolved fractions sum to one, so the
   channel emits exactly the recombinations `alpha_B(He II)` feeds it, and
   the sum of the coefficients is compared with `alpha_B(He II)` and printed.
   *As gated* (`tests/hei_cascade/check_hei_branches.f90`): closure and
   positivity, the fractions against the source coefficients recomputed
   independently, the ordering `2^3S > 2^1P > 2^1S`, a bound on the distance
   from the statistical weights, the sign of the temperature trend, and
   clamping outside 3e3-2e4 K. Restoring the withdrawn `6.23e-14` reading
   fails three of those. The `alpha_B` comparison is printed and not
   enforced, since a shortfall is a statement about two independent data
   sets. `tests/hei_cascade/run_check.sh` runs this with the two-photon gate.
4. **Paired production runs**, HII20 and HII40, at 68 ranks with the same
   Sobol seed, after step 1 and again after step 2, so each step's effect on
   `r(He0=0.5)`, `r(C2+=0.5)`, `r(H0=0.1)` and the line ratios is separately
   measured against the 0.002% seed-to-seed noise floor.
5. **Nebular-continuum consistency**: the He I two-photon component of
   `_nebcont.txt` and the diffuse packets must now derive from one table.
   Confirm the continuum output is unchanged, since `nebcont_mod` already
   used the correct shape.

## 7. Deliverables

- The plan's two steps as two commits, each with its measured effect.
- `tests/hei_cascade/` gate covering section 6 items 1-3.
- Updated `HII40_CLOUDY_RESIDUAL_INVESTIGATION.md` with both rows.
- If the production configuration changes, regenerate
  `results/continuous_energy/production/`, Tables 4/5 and Figures 14/15, and
  update `MoCHII_paper/mochii.tex`, `docs/MoCHII_physics.tex`,
  `docs/MoCHII_UserGuide.tex` and
  `docs/PHOTON_ENERGY_SAMPLING_CODE_COMPARISON.tex` --- the same sweep the
  Milne change required.

## 8. Out of scope

- Collisional redistribution among the `n = 2` terms at densities above the
  benchmark values; `hei_metastable_mod` already carries `q_31a`/`q_31b` for
  the metastable diagnostic and coupling them into the diffuse branching is a
  separate density-dependent change.
- He I resonance-line transfer. The 584 A photon is treated here as escaping
  or ionizing hydrogen on emission, as now. **Superseded and now done**: that
  is the `p -> 1` limit of a local branching that Wood's `ionize` and
  CMacIonize both resolve, and on the benchmark grids the correct value is
  0.117 (HII40) and 0.0044 (HII20). Step 3,
  [HEI_584_CONVERSION_PLAN.md](HEI_584_CONVERSION_PLAN.md), implemented that
  branching as `hei_584_ionizes_H` and measured it: `r(He0=0.5)` moves inward
  by 0.00363 pc in HII40 and 0.00022 pc in HII20, the largest of the four
  corrections to what the diffuse field emits, and it also removed the
  `nebcont_mod` / `diffuse_mod` disagreement about the same decay. Full
  resonance-line transfer remains out of scope there too, since the split is
  an escape-probability closure whose domain of validity --- an optically
  thick, static, dust-free 584 A line --- is documented at the code site.
- Any change to `alpha_B(He II)` itself.
