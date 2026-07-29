# Plane-parallel (xy-periodic slab) illumination and transport — plan

Date: 2026-07-21.  Status: PP1-PP6 and the follow-ons (QMC slab launch,
scattered-escape boundary tally) are IMPLEMENTED and committed (ff8a4e8).
A code-review pass then hardened the mode (2026-07-24): the boundary output
now labels the always-valid L_escape_top / L_escape_bottom and only splits
them into reflected/transmitted when exactly one face is lit; a near-grazing
ray (kz -> 0) is guarded on both the source side (beam theta strictly < 90,
isotropic mu floored away from 0) and the transport side (a max-step counter
in every periodic ion walk terminates a wrapping packet safely and reports
once); the log no longer mislabels a slab as an external isotropic field;
setup validates each lit face's angle and source strength, slab_nmu >= 1,
the namelist grid dimensions, and rejects combining a slab with internal or
external sources; and tests/slab/check_slab.py auto-checks the analytic gates.
Remaining OPTIONAL items: composing slab illumination WITH internal point
sources in one run (needs a source-model decision); the Chandrasekhar
conservative-scattering boundary gate G-c (now possible with the scattered-
escape tally).
Goal: let MoCHII model a plane-parallel slab whose density/temperature vary
only in z, periodic in x and y, illuminated by a collimated beam (at a chosen
incidence angle) and/or an isotropic field entering the top and/or bottom
face, with each face independently specified — plus internal or external
sources inside the same periodic slab.  Natural observable: the emergent
intensity I(mu) at the two boundaries, and z-profiles of the gas state.

Reference: `~/LaRT/combine/LaRT_v2.00` and MoCafe both implement xy-periodic
slabs, but **only for internal sources** (the user confirmed this).  So the
transport wrapping is a port; the **external plane-parallel illumination into
a periodic slab is genuinely new** and has no template in those trees.

## 1. What already exists in MoCHII (inventory)

- `par%xy_periodic`, `par%xyz_symmetry`, `par%z_symmetry` are declared
  (inherited from MoCafe), but the periodic/symmetry wrap is wired **only**
  into the grey-dust routines `raytrace_to_tau_amr` / `raytrace_to_edge_amr`
  (the MoCafe SED/dust path).  On an x/y face exit (iface 1-4) they translate
  the position by `+/- xrange` / `+/- yrange` and re-find the leaf; z faces
  escape.  This is the proven pattern to copy.
- The **ionizing-band transport does NOT honor periodicity at all**:
  `raytrace_ion_to_edge_amr`, `raytrace_ion_to_tau_amr`,
  `raytrace_ion_tau_only_amr`, their DDA Cartesian twins
  (`raytrace_ion_*_car`), `transport_ion_packet`, the diffuse launch, and the
  peel walks all simply terminate when a packet leaves the box.  Wiring
  periodicity into these is the core of the project.
- Box shape is already selectable in `grid_mod_car`: `z_symmetry` (or
  `xyz_symmetry`) sets `zmin = 0` (the `(0, zmax)` half-slab); the default is
  the full `(-zmax, zmax)`.  Same for x/y.
- The Cartesian ('car') grid is `nx x ny x nz` with the medium from
  `density_file` (a 3D cube), `nH_const`, or the namelist.  A z-only medium is
  trivial: a horizontally uniform column.
- External field: `emit_external_rec` injects a Lambert (mu = sqrt(xi))
  isotropic field over all six faces, area-weighted; `emit_external_sph` over a
  bounding sphere.  Each carries `from_external` + the surface normal
  `snx/sny/snz` for the peel Lambert weight.  Entering power L = pi I A, with
  the specific intensity I taken from `ext_intensity` (F = pi I).
- Composable sources already exist: `par%nsource` internal points + the
  external field, split by a component CDF (`ion_band_mod`), with QMC twins.

## 2. Reference behavior (LaRT / MoCafe, internal-source only)

- `par%xy_periodic` == slab.  Transport wrap: on an x/y face exit, translate
  the position by `+/- xrange` / `+/- yrange` and re-find the leaf; z faces
  escape.  Identical to MoCHII's grey-path block.
- The slab volume is stored in full (periodicity does NOT fold it, unlike the
  reflecting symmetries that fold by 2/4/8).
- Slab launch (`plane_atmosphere`): a collimated beam straight down from
  `z = zmax` (cost = -1), entry at x = y = 0 (any xy point is equivalent under
  periodicity for a horizontally uniform medium).  This is an INTERNAL-source
  style injection at the boundary; there is no external-field face integral.
- Slab output (`output_sum_rect`): the emergent field is a histogram of
  escaping photons at the two z boundaries — 2 faces, solid angle 2 pi, no
  1/d^2 or area factor.  Luminosity is normalized as "1 photon/cm^2" times the
  slab area.

## 3. Gaps to close (new work)

1. **Transport wrapping in the ionizing band** (port + DDA extension).  The
   octree walks copy the grey-path block; the DDA Cartesian walks need the
   wrap folded into the Amanatides-Woo stepper (re-init at the wrapped face).
2. **External plane-parallel illumination into a periodic slab** (new).  A
   collimated beam at incidence angle theta and/or an isotropic field, entering
   the top and/or bottom face, each face independent.  Not present in
   LaRT/MoCafe.
3. **Slab source geometry + parameters** and its composition with internal
   sources.
4. **Plane-parallel output** — the emergent I(mu) boundary tally.  Point-
   observer peel-off is ill-defined under periodicity (infinitely many
   periodic images), so the boundary angular histogram replaces it.

## 4. Design

### 4.1 Transport wrapping (the enabler)

**Octree ion walks** (`raytrace_ion_to_edge_amr`, `_to_tau_amr`,
`_tau_only_amr`): insert the same periodic block the grey path already uses —
on `il_new <= 0` with `par%xy_periodic` and `iface <= 4`, translate the
position by `+/- xrange` / `+/- yrange` and `il_new = amr_find_leaf(...)`.  z
faces (iface 5/6) escape.  This is a direct copy of the validated pattern.

**DDA Cartesian walks** (`raytrace_ion_to_edge_car`, `_to_tau_car`,
`_tau_only_car`): the incremental Amanatides-Woo walk terminates at a face.
For a UNIFORM grid the DDA parametric coordinates (`tmax`, `tdel`, `t_cur`) are
translation-invariant, so a periodic crossing needs **no re-initialization**:
when the `ix`/`iy` index steps out of `[0, n)`, simply wrap it with
`modulo(ix, nx)` (or `iy`) and continue — `tmax`/`tdel` already point at the
next face one cell on, which is the correct parametric distance in the periodic
tile.  z faces still escape.  The only extra step is in `_to_tau_car`, which
sets a final position from the accumulated (un-wrapped) path length: fold that
`photon%x`, `photon%y` back into the box with the `floor` wrap so it matches the
wrapped index cell.  This is exact for any lateral resolution (`nx = ny = 1` a
crossing returns the ray to the same column; `nx, ny > 1` enters the opposite-
side column) and cheaper than a re-init.

**Scattering walk** (`transport_ion_packet`): it calls `raytrace_ion_to_tau`;
verify the post-step position update and leaf re-find survive the wrap (the
photon's stored `x,y` become the wrapped values).

**Diffuse packets and peel-tau-only walks**: launched/traced inside the slab,
they must wrap too (otherwise the diffuse field and any boundary tally leak at
the x/y faces).

**z boundaries**: both z faces are always OPEN (a packet reaching `+/- zmax`
escapes and is tallied at that boundary).  No reflecting base — an albedo = 1
floor is too artificial for an astrophysical slab (decision in section 7).
`z_symmetry` is used, if at all, only to place the box at `(0, zmax)`; it does
not turn the base into a mirror in a slab run.

**Edge case**: a ray with `kz` exactly 0 (grazing) would wrap in x/y forever.
Guard `dda_init` / the walk with a wrap-count cap (or reject exactly-horizontal
directions at launch); document.  The analytic edge walk must also cap the
accumulated optical depth (reuse the existing `tau > 500` early skip) so a
non-absorbing periodic column cannot loop.

### 4.1a Internal sources under periodicity: interpretation, not extra logic

The transport wrap of 4.1 is **source-agnostic**.  An internal source emits a
packet at its position; that packet then wraps at the x/y faces exactly like
any other packet.  There is **no separate "periodic source" code** — the
source is emitted once, and periodicity is purely a transport boundary
condition.  So internal sources add no implementation complexity beyond PP1.

The only thing to be clear about is the **physical interpretation**, which is
the user's modeling choice:

- `xy_periodic` + an internal point source == a horizontally infinite,
  periodic **array** of identical sources (a disk / ISM slab with a
  statistically uniform source distribution).  This is a valid model and is
  exactly how LaRT/MoCafe use it.
- If instead a **single isolated** internal source is intended, `xy_periodic`
  is the wrong boundary condition (it turns the one source into a lattice);
  use a wide box or a reflecting symmetry instead.

So the recommendation matches the user's intuition: for the clean, unambiguous
plane-parallel case, the illumination is the **external** field of 4.3 (a
beam / isotropic field entering the faces), and `xy_periodic` is the correct
boundary condition for it.  Internal-source periodicity is available for free
once PP1 exists, but should be used only when the periodic-array
interpretation is the one intended.

### 4.2 Slab geometry and the medium

Box shape: `(0, zmax)` via `z_symmetry` (already sets `zmin = 0`);
`(-zmax, zmax)` is the default.  No new box code needed; only clarify the
parameter meaning in the UserGuide.

Two medium regimes share the identical transport (the wrap of 4.1 does not
care which):

- **True 1D plane-parallel** — the medium varies only in z.  Then
  `nx = ny = 1`, `nz = N` is the minimal grid: xy-periodicity makes the lateral
  extent irrelevant, so a single column suffices.  Medium from a
  `density_file` `(1, 1, N)` cube or a small analytic z-profile (a
  `tauhomo`-style uniform column, already partly in `grid_mod_car`).  This is
  the classic plane-parallel atmosphere.
- **General xy-periodic slab** — an ARBITRARY 3D medium
  `nH(x, y, z)` (turbulent / clumpy / multiphase) with `nx, ny > 1`, tiled
  periodically in x and y.  Here `xy_periodic` is a genuine horizontal boundary
  condition, not a statement that the medium is uniform: the single simulated
  tile `[xmin, xmax) x [ymin, ymax)` repeats to a horizontally infinite slab.
  Medium from a full `(nx, ny, nz)` `density_file` cube.  This is the more
  general and often more interesting case (a periodic box of ISM under
  external or internal illumination), and it uses exactly the same wrap,
  sources, and output as the 1D case.

So `nx = ny = 1` is the principle for a z-only medium, but it is a special
case: the plane-parallel machinery is written for general `nx, ny` and the
lateral resolution is simply the user's choice of how finely the periodic tile
is resolved.  The emergent `I(mu)` of 4.4 is then the tile-averaged (i.e.
horizontally averaged over one period) plane-parallel intensity.

**File-input medium under periodicity (a concrete blocker to remove).**  The
medium can already be loaded from a file on both backends — an AMR leaf list
(`par%amr_file`, grid_type = 'amr') or a car density cube
(`par%density_file`, grid_type = 'car', FITS/HDF5 `(nx, ny, nz)`) — and both
must work under `xy_periodic`.  The file read itself is fine (the nH load and
the box bounds `xmax/ymax/zmax -> xrange` used by the wrap are already set on
every path).  The single obstacle is that the one grid builder,
`grid_mod_amr::grid_create_amr` (called by `main.f90` for BOTH grid types),
**force-resets `par%xy_periodic = .false.`** (together with the two
symmetries) near the end of the build.  So today periodicity is silently
disabled for every grid, file or namelist.  The fix is to stop forcing it off
and honor the namelist `par%xy_periodic`:
- For the **car cube**, `nx, ny > 1` with the wrap gives a general periodic 3D
  tile; `nx = ny = 1` gives the z-only column — no other change.
- For the **AMR leaf list**, the wrap relocates the position and
  `amr_find_leaf` finds whatever leaf sits at the wrapped point (exactly how
  the grey-dust path and LaRT/MoCafe already do AMR periodic).  The user is
  responsible for supplying a laterally periodic-consistent refinement (the
  leaf layout at `x = xmin` should tile with `x = xmax`); the code only
  relocates and re-finds, it does not enforce tiling.
This reset removal is part of PP1/PP2; the z faces stay open (the symmetry
resets can remain, since a slab uses no reflecting boundary).

### 4.3 Slab illumination (the new external source)

A new emission routine `emit_slab` (with a QMC twin), modeled on
`emit_external_rec` but restricted to the z faces and with a collimated option:

Parameters (names to be finalized with the author):
- `par%slab_faces` in {`top`, `bottom`, `both`}  (which z face(s) are lit).
- `par%slab_top_mode`, `par%slab_bot_mode` in {`beam`, `isotropic`}.
- `par%slab_top_theta`, `par%slab_bot_theta` [deg from the inward normal] for
  the beam (0 = normal incidence).
- `par%slab_top_phi`, `par%slab_bot_phi` [deg] optional beam azimuth (default:
  random, i.e. an azimuthally symmetric cone of one ray direction — for a true
  single beam, fix it).
- `par%slab_top_source`, `par%slab_bot_source` — the source strength of the
  face, interpreted per mode (see the bookkeeping below):
  isotropic -> the specific intensity I [erg/s/cm^2/sr];
  beam -> the beam flux F_n normal to the beam [erg/s/cm^2].

Emission:
- **beam**: fixed `mu = cos(theta)`; direction downward (`kz = -mu`) for the
  top face, upward for the bottom; `kx, ky` from theta/phi; entry `(x, y)`
  uniform over the face (or the single column when `nx = ny = 1`).  The face
  flux is `F_z = mu F_n` from the input beam flux F_n (see the bookkeeping).
- **isotropic**: Lambert `mu = sqrt(xi)` into the inward hemisphere, exactly as
  `emit_external_rec` but only for the chosen z face(s).
- `from_external = .true.`; `snx/sny/snz` = the face inward normal.

Entering-power bookkeeping — **the input is a source strength; the flux
crossing the horizontal face is the derived normalization that sets the
absorbed energy** (state this explicitly; it is where plane-parallel codes
differ from a total-luminosity code):

- The single quantity that sets all physics is `F_z`, the energy flux crossing
  the (horizontal) illuminated face per area per time [erg/s/cm^2], obtained
  from the input source strength by the cosine-weighted solid-angle integral
  `F_z = integral_hemisphere I_in cos(theta) dOmega` of the incident field:
  - **isotropic**: constant specific intensity I over the inward hemisphere,
    so `F_z = I * 2*pi * integral_0^(pi/2) cos sin dtheta = pi I`.
  - **beam**: a single direction at angle theta (mu = cos theta) with beam
    flux F_n normal to the beam, so `F_z = mu F_n = F_n cos(theta)`.
  So the source strength is input per mode (I or F_n), and the geometric
  factor (pi for the hemisphere, cos theta for the beam) produces F_z.  The
  cos(theta) is physics, not a hidden convention: at fixed beam flux F_n, a
  more grazing beam deposits less energy per horizontal area.
- `F_z` alone fixes the physics: it is the normalization that sets the total
  absorbed energy, the ionization/heating structure, and the emergent
  intensity.  The total simulated luminosity `L = F_z * A_tile`
  (A_tile = the illuminated face area of the ONE periodic tile =
  `xrange * yrange`) is an **internal MC bookkeeping quantity only** — it sets
  the packet weight `Lpacket = L / N_photon = F_z * A_tile / N_photon`.  It is
  NOT a physical input: the physical slab is horizontally infinite, so any
  total luminosity is meaningless; only F_z per unit area is.
- **A_tile cancels in every physical output.**  A leaf column of lateral area
  `a` receives energy `F_z * a` (its share of the entering flux), carried by
  `N_photon * a / A_tile` packets of weight `F_z * A_tile / N_photon`; the
  product is `F_z * a`, independent of A_tile.  The volumetric absorption is
  then `F_z * a / (a * dz) = F_z / dz` — set by F_z, independent of the tile
  size and lateral resolution.  So every per-volume / per-atom result (Gamma,
  heating, x_HI(z), Te(z)) scales linearly with F_z and is invariant under
  changing `nx, ny, xmax, ymax` at fixed source strength.  This invariance is
  a validation gate (G-g below).
- Relation to the existing external field: the isotropic case with
  `F_z = pi I` is exactly the existing `ext_geometry = 'rec'` one-face
  convention (its `pi * ext_intensity * A`, i.e. `ext_intensity` plays the
  role of the specific intensity I).

Composition: a lit face is one more component in the existing source CDF, so
it can be combined with internal `nsource` points and with the other face
(independent flux/angle).  Both faces lit = two components; the result is the
superposition of two single-face runs (a validation gate).

### 4.4 Output: emergent intensity I(mu)

- **Boundary angular tally** (replaces peel under periodicity): bin each
  escaping packet by (top / bottom, `mu = |kz|`) into `I(mu)`
  [erg/s/cm^2/sr] — 2 faces, hemispherical `mu in (0,1]`, no `1/d^2` or area
  factor (the LaRT `output_sum_rect` slab convention).  ALLREDUCE and write
  `<base>_slab_Imu.txt` (mu grid + I_top(mu) + I_bot(mu)).
- **Integral diagnostics, normalized by F_z**: reflected flux (escaping the
  illuminated face) / F_z = albedo, transmitted flux (the far face) / F_z, and
  absorbed flux / F_z, all per unit area and independent of A_tile.  Energy
  conservation per unit area: `F_reflected + F_transmitted + F_absorbed = F_z`
  (the budget must close, and it is stated in units of F_z so the horizontal
  extent drops out).
- **z-profiles**: already available on the car `nz` axis; the Python readers
  (`RatesData`, `EmisData`) slice x_HI, Te, n_e, J, and line emissivity along
  z with the existing `slice`/`plot_field_slice` methods.
- **Peel-off**: deferred.  The point-observer peel is disabled when
  `xy_periodic` is on (a setup guard + note), because the periodic images need
  a careful special-case design; peeling a periodic slab is probably feasible
  but is out of scope here.  The boundary `I(mu)` tally is the plane-parallel
  observable.

## 5. Validation gates

- **G-a  collimated beam, pure absorption**: normal or angled beam into an
  absorbing slab; transmitted flux vs analytic `F exp(-tau/mu)`; I(mu) at the
  far face a spike at the beam mu attenuated by `exp(-tau/mu)`.
- **G-b  isotropic illumination, pure absorption**: emergent transmitted field
  vs the analytic slab (exponential-integral `E_3(tau)` two-stream), or a 1D
  reference.
- **G-c  conservative isotropic scattering (albedo = 1)**: reflection and
  transmission I(mu) of a finite slab vs the Chandrasekhar H-function /
  published plane-parallel benchmark — the gold-standard RT gate.
- **G-d  periodic self-consistency**: an internal source in a periodic slab
  reproduces the horizontally averaged result of a wide non-periodic box
  (validates the wrap).
- **G-e  photoionization slab**: a 1D ionization front illuminated from one
  face vs the existing 1D photoionization reference (ties to the gas physics
  already gated).
- **G-f  two-face superposition**: both faces lit at independent angles/flux
  equals the sum of the two single-face runs.
- **G-g  flux-normalization invariance**: at fixed source strength (I or F_n,
  hence fixed F_z), the per-volume gas structure (x_HI(z), Te(z), Gamma(z))
  and the F_z-normalized emergent I(mu)/F_z and albedo are invariant under
  changing the tile area and lateral resolution (`nx, ny, xmax, ymax`), and
  scale linearly with F_z.  This is the direct test that F_z, not the total
  luminosity, is the normalization.

## 6. Staged order

- **PP1**  [DONE 2026-07-24]  Transport wrapping in the ion band: the octree
  walks (`raytrace_ion_to_edge/tau/tau_only_amr`) got the xy-periodic block, the
  DDA car walks got the index `modulo` wrap + final-position fold, and the
  forced `xy_periodic = .false.` reset in `grid_create_amr` was removed (namelist
  value honored).  Verified: (a) tall thin periodic column ionizes fully
  (`<x_HII> = 0.9999`) vs the isolated control (`0.7072`) — the wrap fires
  decisively; (b) shared vs dda under periodicity agree to `8.5e-17` (both wraps
  correct); (c) off-path (non-periodic) unchanged, dda vs shared `1.2e-15`.
  Source-agnostic (see 4.1a).  Still to do in a later pass: the diffuse/peel
  walks under periodicity (only matter once slab sources / dusty scattering are
  used) and the `kz = 0` grazing guard.
- **PP2**  [PARTIAL 2026-07-24]  Slab box + z-only medium + minimal grid.
  DONE: the `par%nx < 2` setup guard was relaxed to allow `nx = 1` under
  `xy_periodic`; **file input verified under periodicity on BOTH backends** —
  a z-varying car density cube (`density_file`, nH read correctly, wrap fires:
  periodic more ionized than the non-periodic control) and an AMR leaf list
  (`amr_file`, octree walk, small-box wrap fires: `<x_HII>` 0.9998 vs 0.9889).
  RESOLVED and FIXED (via the PP3 analytic gate): the `nx = 1` dda-vs-shared
  discrepancy was an **octree periodic-wrap FP photon loss** (`x -/+ xrange`
  landing outside the box -> `amr_find_leaf` drops the packet), now fixed by the
  `slab_wrap_xy` exact-face re-entry (see the PP3 entry).  Both walks now match
  the analytic at nx=1 (dda == shared to 1.3e-18); `grid_type='amr'` and
  `car_walk='amr'` are also correct.  Box `(0,zmax)` vs `(-zmax,zmax)`: `grid_create_amr`
  always builds the symmetric `(-zmax,zmax)`; `(0,zmax)` is a cosmetic origin
  shift, deferred (physics identical).
- **PP3**  [DONE 2026-07-24, core]  `emit_slab` + `slab_setup` in
  `ion_band_mod`, `source_geometry = 'slab'`, and the `par%slab_*` parameters
  (faces top/bottom/both; per-face mode beam/isotropic, theta, phi, source).
  Normalization implemented as designed: input a source strength (isotropic I /
  beam F_n), the derived face flux `F_z = pi*I` or `mu*F_n` sets
  `L_slab = F_z * A_zface` and the packet weight.  Verified: (a) **beam analytic
  gate G-a** — frozen neutral slab (`gas_niter=0`), `Gamma_HI(z)` is a PERFECT
  single exponential (R^2 = 1.00000), and at `theta = 60` (mu = 0.5) the decay
  slope is EXACTLY `slope0 / mu` (ratio 1.0000) — the `exp(-tau/mu)` law and the
  cos-theta convention are correct; (b) isotropic mode gives `F_z = pi*I`
  exactly and a shallower monotone profile (angular spread); (c) both-faces
  superposition is symmetric (`L_slab = (F_z_top + F_z_bot) A`).
  **nx=1 octree bug FOUND AND FIXED (was the PP2 open item):** the analytic gate
  first showed the octree walk (`car_walk='amr'` / `grid_type='amr'`) ~10x too
  steep at nx=1 while the DDA matched exactly.  Root cause (found by instrumenting
  the edge walk): the octree periodic wrap did `x -/+ xrange`, which can land a
  rounding-epsilon OUTSIDE `[xmin,xmax]`, and `amr_find_leaf`'s strict bounds
  check then DROPS the packet — a photon loss that, at nx=1 where every lateral
  crossing wraps, depletes the deep flux and over-attenuates ~10x.  Fix: a
  `slab_wrap_xy` helper re-enters the packet at the OPPOSITE face's EXACT boundary
  value (`xmin`/`xmax`/`ymin`/`ymax`) instead of `x -/+ xrange` (used in all three
  octree ion walks).  After the fix: nx=1 shared matches the analytic (ratio
  0.9998, R^2 1.0, was -2.386 -> -0.2389); the nx=1 point-source discrepancy is
  gone (dda 0.3766 == shared 0.3766 to 1.3e-18); nx=8 periodic unchanged
  (shared vs dda 8.5e-17); non-periodic untouched (wrap guarded by xy_periodic).
  So the octree is now correct at nx=1 too.  Deferred: the QMC twin
  (`emit_slab_qmc`; `sobol` is guarded off for slab in setup), and the isotropic
  `E_3(tau)` analytic (the smoke already shows the right behavior).
- **PP4**  [DONE 2026-07-24, core]  Boundary `I(mu)` output + energy budget.
  A rank-local `slab_Iesc(nmu, 2)` tally in `jtally_mod` accumulates each
  packet's surviving luminosity (`expo*wl`) at a z-face escape (top +z / bottom
  -z, binned by `mu = |kz|`), added in both octree and DDA edge walks;
  ALLREDUCEd and written to `<base>_slab_Imu.txt` (mu, I_top, I_bot in
  erg/s/cm^2/sr with `F = 2 pi int I mu dmu`) with the reflected/transmitted/
  absorbed budget in the header (`par%slab_nmu`, default 20).  Verified against
  the analytic transmission: normal beam transmits `exp(-tau_perp) = 0.8874`
  EXACTLY (budget closes 0.88739 + 0.11261 = 1.0), I_bot a single spike at
  mu = 1, I_top = 0 (no scattering); `theta = 60` beam transmits
  `exp(-tau/0.5) = 0.7874` exact, I_bot spike at mu = 0.5; isotropic gives a
  broad I(mu) over all bins (mean mu 0.583, forward-shifted by the `exp(-tau/mu)`
  weighting), budget closes.  Reflection (I_top) stays 0 until dust scattering
  (PP5+); the top-face tally already handles it.  Deferred: gate G-c
  (Chandrasekhar conservative-scattering slab) needs the dust-scattering slab.
- **PP5**  [PARTIAL 2026-07-24]  DONE: **peel guard under periodicity** — setup
  now aborts on `par%ion_peel .and. par%xy_periodic` (the point-observer peel is
  ill-defined under the periodic images; the boundary I(mu) tally is the
  observable); gate **G-f** (two-face superposition) already passed in PP3
  (symmetric, `L = (F_z_top + F_z_bot) A`); internal-source-in-a-periodic-medium
  already works (PP1/PP2, transport wrap is source-agnostic).  OPEN/OPTIONAL:
  composing the slab illumination WITH internal `nsource` point sources in ONE
  run (a slab both externally illuminated and with an embedded source).  The
  multi-component builder (`setup_multi_components`) is clean and extensible, but
  adding the slab as a component needs a source-model decision (a new
  `par%add_slab`-style toggle vs a `source_geometry` restructure) plus QMC
  guarding — deferred pending demand, since both cases work independently.
- **PP6**  [DONE 2026-07-24]  Photoionization slab gate G-e + docs.  G-e: a
  full ionization solve (case B OTS) builds a plane-parallel front at the
  Stromgren depth `d = Phi/(nH^2 alpha_B)` — front position within ~5-7% (alpha_B
  + finite-front-width + convergence) and `d` scales with `F_z` (ratio 1.94 for a
  x2 change, the deficit the ~constant front width).  Docs: UserGuide gained the
  slab parameter table (`slab_faces/_mode/_theta/_phi/_source/_nmu`), the
  `source_geometry='slab'` note, and the `<base>_slab_Imu.txt` output row;
  physics gained algorithm sec 4.12 (xy-periodic wrap incl. the exact-face fix,
  the F_z normalization equation, boundary I(mu)) and validation sec 5.13 (all
  the PP1-PP4 gate numbers).  Both docs switched to my_memo.cls (per the user)
  and rebuilt clean (UserGuide 22 pp, physics 31 pp; no overfull, refs resolved,
  new pages inspected).  Test inputs kept under `tests/slab/`.
- **QMC for the slab (`emit_slab_qmc`)**  [DONE 2026-07-24]  the sobol launch
  now covers the slab in both modes (6-dim layout: frequency, face, xy, mu, phi;
  `ion_qmc_ndim` returns 6, a slab branch in `gen_ion_photon_qmc`, the setup
  sobol guard removed).  Verified: isotropic slab transmission random 0.80503 vs
  sobol 0.80511 (0.01%), beam identical (0.88739) --- same physics, sobol enabled.
- **Scattered-escape boundary tally**  [DONE 2026-07-24]  `slab_bnd_add` is now
  also called at the z-face escape of the scattered flights (`raytrace_ion_to_tau`
  octree + DDA), so a dust-scattering slab's diffuse/reflected component enters
  the boundary $I(\mu)$ (the edge walk = direct, to_tau = scattered, no double
  count; the peel tau-only walk is untouched).  Verified on a pure dust-scatter
  slab (transparent ionized gas + global_dgr dust): reflection went from exactly
  0 to nonzero and the refl/tran/abs budget closes to 1.0000 (previously the
  scattered escapes were misattributed to absorption); non-scatter gates
  unaffected (beam transmission 0.88739 unchanged, to_tau not called without
  dust scattering).  This is the prerequisite for the Chandrasekhar G-c gate.

PP1 is the prerequisite for everything and delivers value on its own; PP3-PP4
are the genuinely new plane-parallel-illumination contribution.

## 7. Decisions

All decided (2026-07-21):
- **Normalization convention** — input a source strength (isotropic: specific
  intensity I; beam: beam flux F_n normal to the beam); the derived horizontal-
  face flux `F_z` (`= pi I` isotropic, `= mu F_n` beam) is the sole physical
  normalization and sets the absorbed energy.  The `cos(theta)` for the beam is
  explicit physics.  See 4.3.
- **Parameter surface**: `ext_geometry = 'slab'` + `slab_*`, composing with the
  existing external-field plumbing (CDF, QMC, spectrum resolution).
- **z boundaries are always OPEN (escape), both faces.**  No reflecting base:
  an albedo = 1 floor is too artificial for an astrophysical slab.  So a slab
  run does not use `z_symmetry` for a mirror; the box is `(-zmax, zmax)` or,
  with `z_symmetry` used only to place the box at `(0, zmax)`, the z faces
  still escape.  (This removes the reflecting-base option from 4.1.)
- **Default plane-parallel grid**: `nx = ny = 1` (a z-only medium column); a
  general 3D periodic tile sets its own lateral resolution (4.2).
- **Observable**: the boundary `I(mu)` tally.  Peel-off under periodicity is
  deferred — it is probably feasible but needs careful design for the periodic
  images (a special case), so it is out of scope for this plan.
