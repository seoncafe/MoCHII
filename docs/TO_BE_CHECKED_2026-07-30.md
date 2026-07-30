# Open items after the He I alpha_1 correction

- Date: 2026-07-30
- Context: the session that corrected the ground-level recombination rate
  (`tools/fitting/fit_alpha1_milne.py`), audited the recombination data across
  the registry, and aligned every analytic gate reference with the code
- Full audit trails: [HII40_CLOUDY_RESIDUAL_INVESTIGATION.md](HII40_CLOUDY_RESIDUAL_INVESTIGATION.md),
  [HEI_584_CONVERSION_PLAN.md](HEI_584_CONVERSION_PLAN.md),
  [COOLING_LOCAL_NE_PLAN.md](COOLING_LOCAL_NE_PLAN.md)

Items are separated by what is known about them. Everything in section A is
backed by a measurement made in this session; the numbers are the measured ones,
not estimates.

## A. Items this session created

### A1. The HII40 Cloudy front residual, now -0.17750 pc, unexplained

Every **single-parameter** candidate is closed by measurement:

| candidate | measured lever | needed to close | verdict |
|---|---|---|---|
| diffuse emission (four corrections) | net -0.00306 pc | --- | wrong direction |
| diffuse transport | 0.062 pc between explicit and OTS | 31-77% in diffuse band luminosity | excluded |
| `sigma_pi(He I)` | 0.27 pc per unit | ~70% | excluded |
| `alpha(He I)` | 0.80 pc per unit | ~23% | excluded; code-to-code spread is 11% |
| metal opacity | 0.103 pc per unit upward, saturating | ~180% | excluded |
| metal recombination continua as a diffuse source | 0.29% of channel 2 | --- | ~0.001 pc |
| He I excited-level ionization | 0.2-0.6% | --- | outside the 13.598 eV band floor |
| source photon budget | agrees with Cloudy to 0.005% | --- | excluded |
| dust | the benchmark runs dust free | --- | excluded by construction |
| `alpha_A` / `alpha_1` / `alpha_B` (He I) vs Cloudy | 0.3 / 0.34 / 1.9% | --- | excluded |
| Cloudy's own front resolution | 0.0106 pc median zone width | --- | excluded |

**The surviving hypothesis** is several few-percent differences compounding
together through the He^0 optical depth. No single-parameter run can produce it
or exclude it, which is why the list above closes without explaining anything.

The signature to explain: x(H^0) agrees with Cloudy to 1-2% through the whole
ionized zone while the x(He^0) ratio grows monotonically outward -- 1.007 at
1.6 pc, 1.015, 1.031, 1.059, 1.109, 1.196 at 3.8 pc. Not a normalization offset,
which would be radius-independent, but helium-specific and cumulative, as it must
be when He^0 carries 84-96% of the opacity above 24.6 eV and 24.4% of the
He-ionizing supply is lost to other absorbers.

**Recommendation: do not open this next.** A multi-parameter hypothesis is hard
to design a test for, and A2 and A3 below will narrow the candidate set before it
is worth trying.

### A2. Recombination forming neutral atoms: a 10-16% code-to-code spread

Radiative recombination only, at 10^4 K, MoCHII (Badnell/CHIANTI) against the
Verner & Ferland (1996) fits MOCASSIN uses:

| forming | ratio |
|---|---|
| C II, C III, N II, Ne II | 0.999-1.016 |
| N I | 0.971 |
| Si IV (Na-like) | 0.932 |
| Ne I | 0.898 |
| C I | 0.861 |
| O I | 0.837 |
| Mg II (Na-like) | 0.791 |

The spread is concentrated in the rates that form **neutral** atoms, which is
where He I also sits at 11%.

**This is a variance between two datasets, not a MoCHII error bar.** For the
metals there is no third arbiter: Cloudy ships `data/badnell_rr.dat` and
`badnell_dr.dat` **md5-identical to MoCHII's** and `source/ion_recomb.cpp` uses
`RR_Badnell_rate_coef` / `DR_Badnell_rate_coef`, so Cloudy and MoCHII agree by
construction and cannot adjudicate. Which one is right for C I, O I, Ne I or
Mg II is undetermined here.

He I is the exception and there MoCHII is demonstrably right: Cloudy does not use
Badnell for the He-like and H-like sequences but derives them level by level from
its own photoionization cross sections through the Milne relation, an independent
first-principles determination. Against it Badnell (MoCHII) is +0.3% and VF96
(MOCASSIN) is +11.0%.

One may reason by analogy that the metal neutrals behave the same way -- Badnell
(2006) is the later calculation, VF96 leans on simplified treatments for
many-electron systems, and He I is where that was actually demonstrated -- but
that is inference, not measurement.

**What would settle it**: an independent determination for each ion. R-matrix
values for C I, O I and Ne I recombination (Nahar & Pradhan and successors), or
CHIANTI totals that do not derive from Badnell. **One ion is enough to start**,
because the direction generalizes: if the R-matrix value sits with Badnell for
C I, the same pattern as He I is established.

**Recommendation: do this first.** It is the only item with a clear next
experiment, and the He I sensitivity of 0.80 pc per unit is already measured, so
the verdict feeds straight into A1.

**Caution that follows from the measured case**: `CLAUDE.md` lists MOCASSIN's
11-case regression suite as a 3D validation reference. Helium comparisons against
MOCASSIN inherit a +11% in `alpha_A(He I)` that has been traced to VF96, so
subtract that before attributing a discrepancy to transport or geometry. It does
**not** mean MoCHII may be wrong there.

### A3. The G2a thermal deviation of 0.773% is not atomic data

After aligning the reference on all three axes (recombination, free-free Gaunt
charge dependence, photoionization cross sections) the pure H/He thermal gate
still reports a median `|dTe|` of 0.773% and a mean bias of +0.790%. Two
candidates were checked and excluded:

- `par%cooling_model='local_ne'` versus the reference's n_e -> 0 Tier-1 fits: the
  code's own printout gives `Ssup(1e4 K, 1e4 cm^-3) = 1.0000`, so H I carries no
  density suppression at a hundred times this test's density.
- Voronov collisional ionization: `ci_dere_ratio` returns 1.0 at the default
  `ci_model='voronov'`, so gate and code were already identical.

The likely remainder is the difference between a 1D exact solve and a
level-6 octree Monte Carlo solution, which is what the gate is meant to measure.
**It is not confirmed.**

**What would settle it**: a level-7 g2a run. The same method resolved G1, where
halving the cell took the deviation from +1.00% to +0.68% and identified it as
discretization.

**Recommendation: second priority.** One run, and it either closes the item or
turns it into something specific.

### A4. Numbered here so the section reads continuously; see C and D

The item that stood here -- `par%recomb_model = 'badnell_mao'` is a historical
name, left alone because renaming it would break input files -- is closed. The
name is now `'badnell_milne'` (section C) and the reason given for keeping it was
false (section D, item 6).

### A5. Thirty-seven namelist parameters are read nowhere

`par` now carries 229 fields, and **36 of them are never referenced outside
`src/define.f90`**; `out_bitpix` brings the count to 37, appearing only inside a
comment in `src/fitsio_mod.f90`. A run reads a value for each from the input
file and does nothing with it. That is worse than a switch that refuses an
unsupported value, because the run completes and reports success while silently
ignoring what was asked for.

They are not transcription mistakes but MoCafe features whose implementing code
was never ported into MoCHII:

| group | count | parameters |
|---|---|---|
| stochastic dust-emission alternative path | 8 | `use_dustemis`, `dust_emission_method`, `dust_niter`, `dust_no_photons`, `dust_tol`, `dust_fast_table`, `dust_nU`, `dust_single_teq` |
| source spatial profiles (Sersic / disk) | 9 | `src_spectrum`, `src_zscale`, `src_rscale`, `src_sersic_index`, `src_reff`, `src_axial_ratio`, `src_boxiness`, `source_zscale`, `source_rscale` |
| all-sky HEALPix interior observer | 5 | `allsky`, `allsky_nside`, `allsky_x`, `allsky_y`, `allsky_z` |
| MPI master-slave scheduling | 3 | `num_send_at_once`, `use_master_slave`, `use_reduced_wgt` |
| polarization and the scattering matrix | 2 | `use_stokes`, `scatt_mat_file` |
| Modified Random Walk | 2 | `use_mrw`, `mrw_gamma` |
| FITS output precision | 2 | `out_bitpix`, `out_bitpix_force` |
| unrelated leftovers | 6 | `sightline_tau`, `nr`, `use_amr_grid`, `cone_opening`, `output_normalization`, `radiation_angular_PDF_file` |

The 34 fields referenced on exactly one line were all checked and each does real
work, so nothing of this kind hides among them. One further field is dead but is
not a user knob: `nscatt_tot`, marked in `define.f90` as not an input parameter,
is zeroed in `setup.f90` and then never accumulated or read.

**Deliberately left in place for now**, unlike the four groups removed in
section C. Nothing in the repository references any of them, so deletion is
trivial, but each would have to be written back when its feature is ported.

**What would settle it**: reject them in `read_input` instead of deleting them.
A parameter that an input file sets but MoCHII does not implement then aborts
with that reason, the silent no-op is gone, and porting the feature later means
deleting one abort rather than reconstructing a namelist entry. This is the
recommended shape whenever the item is opened.

## B. Items that predate this session, unchanged

- **Hard-spectrum (PN/AGN) extension** -- the main open development track.
  `docs/HARD_SPECTRUM_PLAN.md`: band headroom to ~500 eV, registry up to
  O/Ne VI, inner-shell and Auger, MOCASSIN PN and Cloudy gates. Its
  density-dependent Lambda(T, n_e) prerequisite is done.
- **He I 10833 below the Porter n_e floor** -- over-estimated at the ten-percent
  level for n_e < 10. Refuse below the floor, or feed from the explicit 2^3S
  population. Diffuse-medium runs only.
- **SH95 line-table interpolation** -- bilinear-in-log is systematically low by
  0.06-0.13% on L(Hbeta). Replace with a cubic spline or the Storey-Hummer form.
- **Cold-gas cooling boundary (T < 1e3 K)** -- the Tier-1 fits collapse there.
  Extend them or call the n-level solve, and add a [C I] H-impact channel.
- **He I 2^3S metastable, remaining gates** -- H2b Cloudy anchor, H4 EXHALE
  triplet profile. `docs/HEI_10833_PLAN.md`.
- **Plane-parallel slab leftovers** -- composing slab illumination with internal
  sources; the Chandrasekhar H-function gate.
- Broad MPI collective and I/O status wrapper. Further ions as science needs them;
  P and K would reuse the VY95 path.

## C. Closed in this session, no longer open

| item | outcome |
|---|---|
| He I `alpha_1` | **corrected** -- from the Milne relation on the transport cross sections, within 0.35% of Cloudy on all three channels, gated by `tests/milne/check_alpha1` |
| registry-wide recombination audit | **no defect found** -- zero transcription errors, and the older `RR2`/`DR2` forms cannot be migrated because every one sits outside Badnell's coverage and is what CHIANTI v11.0.2 still ships |
| analytic gate references | **aligned** -- three mismatched axes removed, 23 quantities held to the Fortran to machine precision by `tests/rates/check_python_rates.py`, 11 gates passing |
| the metal-opacity "negligible above 13.6 eV" comments | **withdrawn** -- removing it moves the He front 0.615 pc |
| the G1 criterion | **moved to level 7** at +0.68%; the level-6 +1.00% is the discretization it is, and the gate now exits nonzero on failure |
| the G2a 0.14% manuscript claim | **was never a regression** -- rebuilding `29b4a46` gives 0.994%, so the gate always read about one percent |
| the RQMC residual bias near 0.06% | **was the reference** -- 0.011% once aligned, and the radial-profile variance gain rises from 4.1x to 30.2x |
| `par%use_ion_band` | **removed** -- MoCHII has no mode that runs without the ionizing band, so the switch had exactly one legal value and `read_input` aborted on the other; gone from `define.f90`, `setup.f90`, `grid_mod_amr.f90`, `ion_band_mod.f90` and 199 input files |
| `par%use_shared_memory` | **removed** -- never read. A global switch was tried and dropped because ranks accumulate into the tally arrays and one shared window races; sharing is fixed at the allocation site, where read-only structure calls `create_shared_mem` directly. The 177 bare occurrences of the name are *local* variables in `memory_mod_mpi.f90` fed by its optional `shared_memory` argument, which is what made the knob look live to a grep |
| the (a, g) and tau scan knobs | **removed** -- `use_ag_list`, `albedo_list`, `hgg_list`, `use_tau_list`, `tau_list`, the constant `MAX_SCAN`, and the `observer_type` arrays that existed only for them (`scatt_ag`, `scatt_agt`, `direc_t`, never allocated) |
| `par%recomb_model = 'badnell_mao'` | **renamed to `'badnell_milne'`** -- only the fit FORM of `alpha_1` is still Mao & Kaastra's, so `recomb_mod` keeps `rr_mao` while the namelist value now names the two determinations it combines, Badnell's total and the Milne ground level. The name also records that `alpha_B` here is a difference rather than a fit, which is where the He I defect entered. Twenty mentions, no namelist migration |
| the `recomb_model` dispatch | **made exhaustive** -- the six coefficient functions read `if hui_gnedin ... else <Badnell>`, so any other string was served the Badnell rate. They now `select case` on both names with no catch-all: a third model added to `read_input`'s whitelist but not to these functions returns NaN instead of silently getting Badnell's |
| the clumpy-medium knobs | **removed** -- 22 `clump_*` parameters plus `use_clump_medium`, `save_clump_info` and `photon_type%icell_clump`. `read_input` accepts only `grid_type` = `car` or `amr`, so a clumpy setup silently produced a uniform medium; the `grid_type` comment claiming `clump` as a legal value and a backward-compatible mapping of the legacy booleans was false and is corrected |

## D. Retractions recorded during this session

Kept as an audit trail; each was written into the documents before being
measured, and each is now marked withdrawn where it appeared:

1. "The Milne integral is the outlier among three `alpha_1(He I)` determinations."
   Wrong -- Cloudy's independent level-resolved value lands within 0.34% of it.
   The error was taking a majority vote between two fits that happened to be
   close, neither derived from first principles.
2. "A ~30% error in the metal opacity would be worth ~0.18 pc." Wrong -- 30% is
   worth 0.031-0.072 pc. The error was extrapolating an on/off bound as if the
   response were symmetric; it is not, 0.615 pc per unit removed against 0.103
   per unit added.
3. "A 10-20% uncertainty in MoCHII's neutral-forming recombination rates."
   Imprecise -- it is a code-to-code spread with no arbiter for the metals, not a
   MoCHII error bar.
4. "Migrate the legacy `RR2`/`DR2` forms to Badnell." Impossible on two counts:
   Badnell has no entry at those (Z, N), and the forms are what CHIANTI's current
   data uses.
5. "The G2a remainder is the reference's free-free Gaunt approximation." Wrong --
   that approximation is accurate to 0.06-0.6% for Z = 1 against the code's own
   tabulated values, and adding the charge dependence moved the median by 0.003
   percentage points.
6. "Renaming `par%recomb_model = 'badnell_mao'` would break existing input files
   and the examples, so the stale name stays." Wrong on the premise -- **no input
   file in the repository sets `recomb_model` at all**. It was reached as a
   default, so the rename to `'badnell_milne'` touched twenty mentions, only two
   of them executable, and needed no namelist migration. The error was carrying
   over the cost of removing `use_ion_band`, which 199 input files did set,
   without checking whether this parameter was set anywhere.
