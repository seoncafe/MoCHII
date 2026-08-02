# Where the MoCHII namelist departs from MoCafe's

Date: 2026-08-02.  Status: audit result.  Two of the names it listed
(`use_master_slave`, `num_send_at_once`) have since been implemented and are
excluded from the count below; nothing else here has changed.

MoCHII's `par` was copied from MoCafe v2.00 and then grew, shrank and shifted
underneath the same names.  Two hazards follow, and this document records both,
because both are silent by default:

1. **Declared but never read** (section 1) — a MoCafe feature whose implementing
   code was not ported.  Setting it does nothing and says nothing.
2. **Same name, different meaning** (section 2) — the name exists in both codes
   and does something in both, but not the same something.  Worse than (1): the
   run does work, just not the work the input asked for.

Section 2 is the reason the file is named for the divergence rather than for the
unread list alone.

---

# 1. Declared but never read

**Supersedes the count in `TO_BE_CHECKED_2026-07-30.md` section A5**, which
found 37.  The list is now 33, and the difference is accounted for:

- `dust_niter`, `dust_no_photons`, `dust_tol` (−3) are **implemented** since
  that audit; they drive the infrared transport passes in
  `src/dust_emis_transport_mod.f90`.
- `save_jlam` (+1) was missed by A5.  It has been comment-only for as long as
  the current history records, so this is a correction to that audit rather
  than a regression.
- `use_master_slave`, `num_send_at_once` (−2) are **implemented** since this
  audit was first written, in `src/photon_schedule_mod.f90`: rank 0 serves
  batches of `num_send_at_once` photon indices on demand instead of every rank
  walking a fixed cyclic sequence, for the stellar loop, the diffuse LyC loop
  and the peel-off imaging pass.  One divergence from MoCafe survives and is
  deliberate — MoCafe defaults `use_master_slave = .true.`, MoCHII `.false.`,
  because with the default `launch_sequence = 'random'` the schedule also
  selects which rank draws which packet, i.e. the realization
  (`MPI_LOAD_BALANCE_2026-08-02.md`).

37 − 3 + 1 − 2 = 33.  A5's own resolution — reject an unimplemented parameter in
`read_input` rather than deleting it — remains the recommended shape and is
repeated under "Open decision" below.

`par` is one derived type and the input file is `namelist /parameters/ par`, so
**every** component is settable from an input file whether or not any code reads
it.  A name in this list can be set in a run, will be accepted without an error,
and will do nothing.  That is the failure this document exists to record: it is
silent, and it is silent in exactly the direction that looks like success.

All 33 of them **are** read by MoCafe v2.00, which is where they came from.  So
the reader who sets one is not inventing a parameter — they are using a MoCafe
parameter in a code that dropped, renamed or restructured the feature behind it.

## Method

Reproducible, and the reason the counts below can be trusted:

1. Extract every field of `type params_type` in `src/define.f90` (230 fields),
   handling multi-name declarations, array specifications and initializers.
2. Strip Fortran comments from every `src/*.f90` — respecting quoted strings, so
   a `!` inside a string literal does not truncate the line.
3. For each field, search the stripped code for `par % <name>` with a word
   boundary (case-insensitive, optional whitespace around `%`), skipping array
   subscripts when deciding whether an occurrence is an assignment.
4. Classify: no occurrence outside `define.f90` = never referenced; occurrences
   that are all assignment targets = written but never read.
5. Repeat steps 2-3 against `MoCafe_v2.00/src/*.f90` for the cross-reference.

`src/` is the whole build: the Makefile compiles `$(SRCDIR)/%.f90` with
`SRCDIR = src` and nothing else, so a name absent from `src/` is absent from the
executable.

Two names, `save_jlam` and `out_bitpix`, do appear once each outside
`define.f90` — both inside comments (`src/jtally_mod.f90:24`,
`src/fitsio_mod.f90:5`), which is why the comment-stripped pass classifies them
as unreferenced.

## Never referenced outside `define.f90` (33)

| group | parameters |
|---|---|
| transport / parallel control (4) | `use_reduced_wgt`, `use_mrw`, `mrw_gamma`, `nr` |
| polarization / scattering matrix (2) | `use_stokes`, `scatt_mat_file` |
| source geometry (9) | `src_reff`, `src_rscale`, `src_zscale`, `src_sersic_index`, `src_axial_ratio`, `src_boxiness`, `source_rscale`, `source_zscale`, `cone_opening` |
| SED mode (2) | `save_jlam`, `src_spectrum` |
| dust emission (5) | `use_dustemis`, `dust_emission_method`, `dust_single_teq`, `dust_fast_table`, `dust_nU` |
| output / observation (7) | `allsky`, `allsky_nside`, `allsky_x`, `allsky_y`, `allsky_z`, `sightline_tau`, `radiation_angular_PDF_file` |
| file I/O (3) | `out_bitpix`, `out_bitpix_force`, `output_normalization` |
| grid (1) | `use_amr_grid` |

## Written but never read (2)

| parameter | assigned at | note |
|---|---|---|
| `lambda0` | `src/sed_mod.f90:94` | inside `setup_sed`, which nothing calls |
| `nscatt_tot` | `src/setup.f90:73` | reset at startup, never consulted afterward |

## The three that will bite hardest

- **`use_amr_grid`** — MoCafe maps it to `grid_type = 'amr'` (`setup.f90:115`).
  MoCHII ignores it, and a run that sets it lands on the Cartesian grid instead
  of the octree without saying so.  `grid_type = 'amr'` is the working spelling.
- **`use_dustemis`** — MoCafe's master switch for dust thermal emission.  MoCHII
  split the feature into `dust_sed` (the SEDust solve) and
  `dust_emis_transport` (transporting the reradiated infrared) and reads
  neither name.  A MoCafe input carried over therefore runs with dust emission
  off.
- **`dust_single_teq`, `dust_fast_table`, `dust_nU`, `dust_emission_method`** —
  described in the User Guide, read by no code.  Documentation and behavior
  disagree here, and the documentation is the optimistic one.

`use_master_slave` was the fourth entry here — a performance switch, so a run
that set it and got plausible timings had no signal that it did nothing.  It is
implemented now, and the remaining hazard is only its default, which is `.false.`
in MoCHII against `.true.` in MoCafe.

## Open decision

Three routes, not mutually exclusive:

**(a) Remove them.**  A parameter that does nothing should not be advertised.
Costs: input files that set them stop parsing (a namelist read fails on an
unknown name), which is loud but breaks old inputs.

**(b) Keep them and refuse them.**  This is A5's recommendation and the one to
follow: in `read_input`, compare each against its default and abort on any a run
actually set, naming it as unimplemented.  The silent no-op is gone, MoCafe
inputs still parse (nothing is deleted from the namelist), and porting a feature
later means deleting one abort rather than reconstructing a namelist entry.  A
warning instead of an abort is the weaker version of the same idea and is worth
considering only for parameters whose neglect cannot change a result.

**(c) Implement them.**  Some are real gaps rather than dead weight:
`use_stokes` / `scatt_mat_file` (polarized scattering), the `allsky_*` block
(all-sky imaging), `use_mrw` / `mrw_gamma` (modified random walk for optically
thick cells).  Each is its own piece of work and should be decided on its own
merits, not as part of this cleanup.

Recommended order: **(b) first** — the silent ignore is the actual defect, and
the warning is cheap and reversible — then (a) or (c) per item.

---

# 2. Same name, different meaning

These names exist in both codes and are read by both.  Carrying an input file
across therefore produces a run that works and is wrong, which no error message
will announce.  Every row below was read out of both source trees, not inferred
from the declarations.

| name | MoCafe v2.00 | MoCHII | how it fails |
|---|---|---|---|
| `rmax`, `rmin` | **geometry**.  `rmax > 0` declares the system a sphere and *overwrites* `xmax/ymax/zmax` (`setup.f90:86-96`), sets the source sampling radius (`gen_photon_car.f90:102`) and the clump placement shell, `rmin` being its inner cavity (`clump_mod.f90:474`) | **a density cut**.  A leaf whose center radius falls outside `[rmin, rmax]` gets `nH = 0` (`grid_mod_amr.f90:142-146`); the grid keeps its size.  `rmax` additionally serves as the entry-sphere radius that `ext_geometry='sph'` requires (`setup.f90:366`) | silent.  A MoCafe input meaning "make the box a sphere of this radius" instead empties the gas outside it on an unchanged box |
| `luminosity` | source luminosity over the whole SED grid | the **ionizing-band** luminosity `[eion_min, eion_max]`; with `add_fuv` the FUV segment rides on top (`define.f90:205-209`, `ion_band_mod.f90:452`) | silent.  The same number means a different total, so the photon budget shifts by whatever fraction of the SED lies outside the band |
| `density_file` | cube values become `grid%rhokap` directly, scaled by `cext_dust` only (`grid_mod_car.f90:141-143`) | cube values are **nH [cm^-3]** per leaf (`grid_mod_amr.f90:109-127`); dust follows through `DGR * cext_dust`, and the gas density is set as well | silent unless `DGR = 1`, and even then the cube now also drives the gas |
| `ext_intensity` | `-999` sentinel; the external field is on when it is set *or* `ext_spectrum` is given (`setup.f90:632`) | `> 0` is itself the on switch (default `-1`), and it is interpreted on the **ionizing segment** of the field (`ion_band_mod.f90:190, 514`) | silent for the on/off edge cases; the normalization differs whenever the field extends past the ionizing segment |
| `ext_spectrum` | a file path, or blank | a file path **or a built-in ISRF preset name** — `'draine'`, `'habing'`, `'mathis'`, scaled by `ext_scale` (`ion_band_mod.f90:530, 891`) | not a hazard in the MoCafe→MoCHII direction (a path stays a path); it is an extension, recorded so the two lists can be compared |
| `launch_sequence = 'sobol'` | the whole launch draw is quasi-random | **only** the frequency and direction of a single internal point source; transport, diffuse and scattering draws stay on the Mersenne Twister (`define.f90:193-199`) | silent.  The variance reduction a MoCafe user expects is narrower here |
| `grid_type` | `'car'` (face-based, non-cubic cells allowed, `geometry`/`density_*` density models apply) / `'clump'` / `'amr'` | `'car'` (octree leaf list + DDA walk, **cubic cells enforced**, `geometry`/`density_*` not applied) / `'amr'`; `'uniform'` is an accepted alias | `'clump'` **aborts** (`setup.f90:90-93`), so the value set is safe; `'car'` is the silent one, because the name survives while its cell model does not |
| `source_geometry` | the source's spatial distribution: `'point'`, `'uniform'`, `'exp_spiral'`, `'sech'`, … (`gen_photon_car.f90:100`) | a **legacy alias only**: `'point'` is inert, `'external*'` is shorthand for the external field, `'slab'` is the one value that selects real behavior.  Internal sources are always point sources | **loud**: any other value aborts at `setup.f90:408`, so this one cannot bite silently |
| `ext_geometry` | `'sph'` / `'cyl'` / `'rec'` | `'rec'` (default) / `'sph'` only | **loud**: `'cyl'` aborts at `setup.f90:361` |
| `dust_niter` | Lucy iteration for dust self-absorption in the SED band (`lucy_mod.f90:67`), default 1 | passes of the infrared transport (`dust_emis_transport_mod.f90:242`), default 2 | same physical meaning (self-absorption iteration), different loop and default; listed for completeness, not as a hazard |

Two names that *were* collisions were resolved on 2026-08-02 rather than
documented:

- `dust_model` meant the SEDust grain model in MoCafe and the leaf dust density
  law in MoCHII.  MoCHII now uses MoCafe's spelling for both: `dust_model` is
  the grain model, `dust_density_law` is the density law.
- `spectrum_type = 'shape'` meant a `lambda [um]` abscissa in MoCafe and an
  `E [eV]` abscissa in MoCHII, with `'shape'` the default on both sides — the
  most dangerous item found, since it read a different spectrum without a word.
  MoCHII's value is now `'shape_ev'`, so a MoCafe input carrying `'shape'` is
  rejected by name.  The `per_ev` / `per_hz` / `per_ang` / `per_um` values carry
  their abscissa in the name and mean the same thing in both codes.

## Reading the two directions

The rows above are written from the MoCafe→MoCHII direction because that is how
inputs actually travel here.  The reverse direction is not symmetric: a MoCHII
input handed to MoCafe fails on the names MoCafe does not declare (`eion_min`,
`abund_*`, `dust_density_law`, …), and a namelist read fails on an unknown name,
so it is loud.

## What this audit did not check

- Only `src/*.f90` was analyzed.  Standalone programs under `tests/` and
  `tools/` have their own `par` assignments and were not swept.
- Whether a parameter is *used correctly* where it is read was not examined;
  this audit only distinguishes referenced from unreferenced.
- The MoCafe cross-reference of section 1 establishes that each name is read
  there.  It does not establish that MoCafe's meaning is the one MoCHII would
  want; that is what section 2 is for.
- Section 2 is not proven exhaustive.  It was assembled by comparing the two
  `params_type` declarations name by name and then reading the use sites of
  every shared name whose declaration, default or documented value set differed.
  A pair that agrees in all three and still diverges deep inside a routine would
  not have been caught.
