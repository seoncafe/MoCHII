# MoCHII — founding design (photoionization on the MoCafe AMR engine)

Written 2026-07-11; package name **MoCHII** (MOnte Carlo for H II regions)
confirmed the same day, developed as a separate package at
`~/MoCHII/MoCHII_v1.00/` (this tree). Direction agreed with K. Seon: instead of
retrofitting an octree into MOCASSIN (see
`RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/OCTREE_PLAN.md` and the
`P1.11-feasibility` analysis of why its traversal resists structural
change), build H/He(+alpha) photoionization **on top of MoCafe v2's
existing, validated AMR dust-RT engine** — by copying the needed modules
into this package (see `docs/PORTING.md`), never modifying the MoCafe tree —
with modern analytic atomic data following the EXHALE pattern
(`RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/docs/`).

Scope: H/He ionization structure and (Stage G2) temperature on AMR grids —
Stromgren-type problems, I-fronts, post-processing of TNG/RAMSES snapshots —
with dust + PAH emission coupled throughout (section 7), and metal lines
added incrementally (section 8).

---

## 1. What MoCafe v2 already provides (the decisive asset list)

| Need for photoionization | Already present | Where |
|---|---|---|
| Octree AMR grid + builder + neighbor lists | yes, validated (AMR ↔ Cartesian sphere) | `octree_mod.f90` (536 lines) |
| Event-driven traversal (interaction-or-crossing, no re-init) | yes — the design MOCASSIN's `pathSegment` lacks | `raytrace_amr.f90` |
| Multi-wavelength transport + per-leaf J_lambda tally | yes (SED mode, Stage 1 complete) | `define.f90`, `jtally_mod.f90` |
| Iterate-to-convergence driver | yes (Lucy dust-emission iteration, Camps+15 validated) | `lucy_mod.f90` |
| Variance reduction | peel-off (NEE) + analytic first-flight + MRW | `peelingoff_mod`, `raytrace_amr`, `mrw_mod` |
| Gas fields from simulation snapshots | `nH`, `xHI`, `Z` columns read (currently transient) | `read_generic_amr.f90`, `grid_mod_amr.f90` |
| MPI, HDF5/FITS, converters (TNG/RAMSES) | yes | `memory_mod_mpi`, `hdf5io_mod`, tools |

The expensive parts of an "AMR photoionization code" — grid, traversal,
frequency-resolved transport, the convergence loop, I/O, parallelism — already
exist and are validated. **The new work is the gas microphysics and its
feedback on opacity**, which is exactly the part that is small and analytic
for H/He (the EXHALE approach).

## 2. What is missing (the actual new work)

1. **Ionizing frequency bins.** The SED grid covers dust wavelengths (um).
   Add an ionizing band nu >= 13.6 eV (~10–20 bins reach percent-level rate
   integrals for H/He; EXHALE integrates the same way). Source spectrum
   (blackbody or table) sampled in the band; the grey `s_ext` rescaling trick
   does not apply — the gas absorption coefficient varies with both cell and
   bin, so the ionizing band carries its own opacity array.
2. **Leaf gas state.** Persist `nH`, `x_HI`, `x_HeI`, `x_HeII`, `T_e`, `n_e`
   as leaf arrays (today `nH`/`xHI` are deallocated after the dust-density
   step). Initialize from the snapshot columns or fully neutral.
3. **Atomic rates (EXHALE pattern, all analytic fits + verification scripts):**
   - photoionization cross sections: Verner et al. (1996) — HI, HeI, HeII;
   - radiative + dielectronic recombination: Badnell (2006);
   - collisional ionization: Voronov (1997);
   - case A/B switch; recombination emissivities later from Storey & Hummer
     (1995) when emission maps are wanted;
   - **nebular continuum: port MOCASSIN's `fb_ff()` block as one unit.**
     Three components, all living together in `emission_mod::fb_ff` +
     `ph_mod::initLogGamma*`; their currency differs (K. Seon, 2026-07-11):
     * free-bound — **confirmed current**: Ercolano & Storey (2006, MNRAS
       372, 1875), by MOCASSIN's own authors, still the definitive
       free-bound reference. HI/HeI/HeII gamma_nu(T) tables
       (`data/gammaHI.dat`, `gammaHeI.dat`, `gammaHeII.dat`, with
       ionization-edge structure) + the log-interpolation in `fb_ff`.
       On port, verify the compact tables in this tree (~21 T x ~22 nu
       points) are the ES06 data; if they predate ES06, swap in the ES06
       machine-readable tables at the same interface (drop-in);
     * free-free — `gaunt.f90::getGauntFF`, Hummer (1988) Chebyshev fits of
       Karzas & Latter (1961), max error 0.7%. **Not the current standard**
       — that is van Hoof et al. (2014, MNRAS 444, 420; adopted by Cloudy).
       In the HII-region domain (Te ~ 5e3-2e4 K, radio-optical) the two
       agree to ~<=1%, so: port `gaunt.f90` verbatim now (free-standing
       ~155 lines, zero dependencies) and schedule a **van Hoof+2014 table
       swap as a small follow-up** — the getGauntFF interface
       (z, log10Te, log10nu array -> g array) matches a vH14 lookup
       exactly, so it is a drop-in; gate the swap with a Hummer-vs-vH14
       difference plot over the HII domain;
     * two-photon — HI 2s-1s and HeI (`data/HeI2phot.dat`) continua.
     Port note: `fb_ff` allocates its logGamma scratch arrays on every call
     (a known MOCASSIN inefficiency flagged in its Phase-B follow-on list) —
     allocate once at setup in the ported version.
     Used twice: free-free/recombination cooling in the G2 thermal balance,
     and the full nebular continuum emissivity for output maps (radio/mm
     free-free and the Balmer/Paschen jumps are primary HII-region
     observables alongside the PAH/IR bands and recombination lines).
4. **Rate integrals from the existing J tally.** For each leaf:
   Gamma_i = integral of 4 pi J_nu sigma_i(nu) / (h nu) dnu over the ionizing
   bins; photoheating H_i uses the same integrand times (h nu - h nu_i).
5. **Ionization equilibrium solve.** Three ionization states (HII; HeII,
   HeIII) + n_e closure: photo + collisional ionization vs. recombination.
   A small nonlinear solve (Newton/bisection on n_e) in each cell — EXHALE has
   this worked out and documented.
6. **Thermal balance (Stage G2).** Photoheating vs. H/He line cooling,
   recombination cooling, free-free (EXHALE `cooling_formulas`); solve for
   T_e together with the ionization state.
7. **Opacity feedback + iteration.** Ionizing-band absorption coefficient
   n_HI sigma_HI(nu) + n_HeI sigma_HeI(nu) + n_HeII sigma_HeII(nu) (+ dust);
   re-transport and repeat until x and T_e converge. This is the Lucy loop
   with one difference: gas opacity changes between iterations, so the
   stellar tally must be recomputed each iteration (no stellar-tally reuse).
8. **Diffuse field.** Stage G1 uses the on-the-spot approximation (case B).
   Stage G3 emits recombination-continuum packets explicitly from each leaf —
   structurally identical to the dust re-emission that `lucy_mod` already does.

## 3. Iteration structure (Lucy loop, generalized)

```
read grid -> leaf state (nH, x, Te) -> ionizing-band opacities
do iter = 1, gas_niter
    transport stellar (+ Stage G3: diffuse) packets  -> J_nu per leaf
    for each leaf: Gamma, H from J_nu
                   solve ionization (+ Stage G2: thermal) equilibrium
    update ionizing-band opacities from new x
    converged when max |delta x_HII| and |delta Te|/Te below tolerance
end do
outputs: x_HI/x_HeI/x_HeII/Te leaf maps, I-front radius, (later) emission maps
```

Under-relaxation on x between iterations if oscillation appears near sharp
I-fronts (known failure mode; MOCASSIN's own iteration never protects against
it and its disabled 1D mode collapses exactly there).

## 4. AMR-specific design points

- **Refinement where it matters:** the I-front is the one place needing
  resolution. Start with the snapshot's static tree; a later stage can
  re-refine on the x_HI gradient between iterations (`amr_build_tree` is
  cheap to re-run; leaf state maps over by position).
- **Memory:** the ionizing band adds (nbin_ion x nleaf) reals for J and for
  opacity — with 10–20 bins this is small next to the existing SED J tally;
  the ROADMAP chunking item applies unchanged if it ever grows.
- **MRW off in the ionizing band.** MRW accelerates scattering-dominated
  diffusion; the ionizing band is absorption-dominated, so MRW must be
  restricted to the dust band.
- **Peel-off still applies** to dust-scattered ionizing photons for images;
  rate integrals come from the volume tally, not the peel.
- **Cartesian (car) grids** (DONE, 2026-07-12/13; grid_type='car',
  'uniform' accepted as alias): three layers.  (a) Geometry is fully
  analytic — cell centers/half-widths and the leaf<->cell maps are
  accessor functions of the raster index (octree_mod: cell_cx/cy/cz/ch,
  leaf_cell, leaf_cx/cy/cz/half); the car build allocates NO tree
  arrays at all (no children/neighbor AND no
  cx/cy/cz/ch/ileaf/icell_of_leaf — the whole ~10 GB at 512^3).  The
  octree path keeps the arrays; all physics modules read geometry only
  through the accessors.  (b) Traversal selects via par%car_walk (not a
  physics choice — both agree to rounding): 'dda' (default) = dedicated
  incremental Amanatides-Woo walks (per-ray tMax/tDelta, comparisons and
  additions only), the fastest car walk — transport x1.64 over shared,
  x2.0 over octree (total Stromgren 13.0 -> 6.5 min); 'shared'
  (verification) = the octree walks with the two geometry primitives
  branched (amr_find_leaf -> integer divisions, amr_next_leaf -> integer
  face arithmetic; amr_cell_exit common), BIT-IDENTICAL to the
  single-level octree (max|dx_HI| = 0 — the cross-check DDA cannot give,
  differing at rounding 2.7e-15 by design; x1.25 over octree, slower
  than DDA).  refine_front stays octree-only.  (c) The grid can
  be built from the namelist without a file: for grid_type='car',
  amr_file empty -> par%nx/ny/nz over the box [-xmax,xmax]... (cubic
  cells).  A gas-density model applies on the leaf centers for BOTH
  car and amr: par%nH_const (uniform density override) and
  par%rmax/rmin (zero nH outside a sphere / inside a shell) — so a
  sharp sphere or shell edge is resolved by the octree while the
  density stays uniform.  Gate: namelist car (nH_const=100, rmax=4)
  reproduces the file-built car sphere bit-identically.

## 5. Staged plan with validation gates

| Stage | Content | Gate |
|---|---|---|
| G0 (3–5 d) | ionizing bins + Verner cross sections + Gamma/H integrals from J (no feedback) | Gamma(r) vs analytic point-source attenuation; cross-section curves vs EXHALE's `compare_photoion_cross_sections.py` |
| G1 (1–2 wk) | equilibrium solve + opacity feedback + iteration (case B, fixed Te) | **Stromgren sphere vs analytic R_S** (iteration- and resolution-convergent — the test MOCASSIN's 1D fails); AMR ↔ uniform-grid match (existing methodology) |
| G2 (2–4 wk) | thermal balance + the generic trace-metal frame (registry, cascade, n-level solver — see section 8) with the first coolants O II/III, N II | Te profile vs MOCASSIN 3D HII20/HII40 (regression cases 01/02, freshly optimized reference — **requires O/N cooling to match**) and vs EXHALE 1D; emissivities of each ion vs PyNeb |
| G3 (~1 wk) | explicit diffuse field (recombination continuum packets); case A | case A vs case B bracketing; energy conservation |
| G4 (1–2 wk) | I-front refinement; TNG/RAMSES post-processing demo | refined vs uniform high-res run |
| G5 (open-ended) | further ions one at a time (S, Ne, C, Ar, ...) — data operations per section 8 | each ion: emissivity ratios vs PyNeb; line fluxes vs MOCASSIN where applicable |

Total: **~9–12 weeks** to a validated gas+dust+PAH HII-region code with the
first diagnostic ion set, each stage independently testable. Compare with
the MOCASSIN octree route (realistically 5–7 months including the traversal
rewrite and spatial decomposition it would eventually require — and without
peel-off, MRW, or a modern dust-emission stage).

## 6. Role of the existing codes

- **MOCASSIN** (freshly optimized: spline cache ~36% gas, hinted DDA ~31%
  dust; 11-case regression suite): the 3D reference for G1/G2 gates. No
  further structural work planned there.
- **EXHALE**: the atomic-data source pattern and 1D reference profiles, with
  ready-made comparison scripts.
- **LaRT**: conventions already shared (distance_unit, memory layout,
  random_mt); resonance-line transfer stays in LaRT.

## 7. Dust + PAH coupling in HII regions (why this combination is the point)

Confirmed direction (K. Seon, 2026-07-11): HII regions carry dust, and PAH
emission is a primary observable — the goal is ionization **added to** the
existing dust RT, not a gas-only code. MoCafe's SEDust stage already provides
the emission side (DL07 = the PAH model, plus astrodust/Zubko; stochastic
heating validated against Camps+15). The couplings to design in:

1. **One radiation field, two consumers.** The same J_nu tally feeds the gas
   rate integrals (Gamma, H) and the grain/PAH heating integral. Extend the
   dust-heating U integral into the ionizing band — EUV heating of grains is
   significant inside HII regions and currently the SED band stops short.
2. **Dust competes for Lyman-continuum photons.** Dust opacity in the
   ionizing band absorbs a fraction of the ionizing photons and shrinks the
   Stromgren radius (the classic dusty-Stromgren correction). This falls out
   automatically once the ionizing band carries gas + dust opacity, and is a
   good quantitative gate: reproduce the analytic dusty R_S reduction.
3. **PAH abundance tied to the ionization state.** PAHs are observed to be
   destroyed inside ionized gas. `laursen09_ndust(nH, xHI, ...)` with
   `par%f_ion_dust` already expresses dust density as a function of the
   ionization fraction — generalize it to use the **computed** x_HI from the
   G1 solve (updated each iteration) instead of the static snapshot column,
   optionally with a separate PAH-specific ionized-gas survival factor.
4. **Science outputs of the combined code:** x/Te structure, Halpha/Hbeta
   recombination-line maps (Storey & Hummer emissivities, peel-off imaging
   already in place), free-free radio/mm continuum maps (Gaunt factors from
   the ported `gaunt.f90`), and PAH/IR band maps from SEDust — gas and dust
   observables from a single self-consistent run.

Stage placement: (1) and (2) belong to G1/G2 (they are just band/opacity
bookkeeping); (3) is a small G2/G4 refinement; the recombination-line maps
in (4) slot in after G2.

## 8. Metal lines — incremental by design (K. Seon, 2026-07-11)

Metal lines are **in scope**, added one ion at a time. Two facts shape the
design:

1. **Physics requires it anyway.** HII-region Te (~8,000 K) is set by metal
   forbidden-line cooling ([O III], [O II], [N II], ...). A pure H/He
   thermal balance overshoots to ~15-20 kK, and the G2 gate (Te vs MOCASSIN
   HII20/40, whose abundance files include C/N/O/Ne/S) cannot pass without
   at least O and N cooling. So the metal frame enters **with G2**, not
   after it.
2. **Trace-species approximation enables incrementality.** At HII abundances
   metals barely perturb n_e or the H/He opacity, so each metal element can
   be solved *after* the H/He+thermal iteration state of each cell:
   consume J_nu and (n_e, T_e), return its ionization fractions, cooling,
   and line emissivities. Adding an ion changes no existing results paths —
   it adds a cooling term and new output lines.

**Design rule: adding an ion is a data operation, not a code operation.**
Three generic pieces, written once:

- **Species registry.** A table of (element, ion, data file, abundance).
  Leaf state grows a compact ion-fraction block indexed by the registry.
  Abundances from a file (MOCASSIN-style abundance input is a fine format).
- **Element ionization cascade.** For each registered element: photo rates
  (Verner+96 covers all ions) + Badnell RR/DR + Voronov collisional +
  **charge exchange with H/He** (Kingdon & Ferland 1996 — essential for O;
  EXHALE's `charge_exchange_table4.md` already tabulates this) -> tridiagonal
  stage solve. Generic over the registry.
- **n-level atom solver.** Build the rate matrix from (A-values, effective
  collision strengths Upsilon(T), level energies/weights), solve for level
  populations, emit cooling + line emissivities. One routine serves every
  ion. (MOCASSIN's `equilibrium` is the reference implementation; a clean
  rewrite is ~150 lines and avoids its known quirks, e.g. the spline issue
  fixed this session.)

**Atomic data: CHIANTI-derived formula fits, built on EXHALE's existing
pipeline (K. Seon, 2026-07-11).** EXHALE already contains the core tooling
in `EXHALE/cooling_data/`:
- `chianti_cooling.py` — a dependency-light direct reader of CHIANTI v11
  ASCII files (.elvlc/.wgfa/.scups, no ChiantiPy) with Burgess-Tully
  descaling, written as "the auditable source for every cooling coefficient
  that gets ported into the Fortran";
- fitted coronal cooling functions Lambda(T) = T^(-1/2) sum_i A_i exp(-T_i/T)
  with coefficients already produced for C I-II, N I-II, O I-II
  (`cno_formula_coefficients.txt`, fit range 1e3-1e5 K), plus comparison
  notebooks; CHIANTI v11.0.2 lives locally at `~/RT_Codes/CHIANTI/dbase`.

The fitting workstream extends this pipeline rather than starting over.
Products, two tiers (the tier split matters because the EXHALE fits are
coronal/low-density-limit — correct for cooling below the critical
densities, but by construction unable to give density-dependent line
ratios like [S II] 6717/6731):

- **Tier 1 — cooling fits (hot loop).** Lambda_ion(T) in the EXHALE form
  for every registry ion (extend from C/N/O to S, Ne, Ar, ...). Fast,
  smooth, Newton-friendly inside the G2 thermal iteration. Where an HII
  case approaches a line's critical density, add a fitted density
  suppression factor (two-level correction) or fall through to Tier 2 for
  that ion.
- **Tier 2 — collision-strength fits + n-level solve (diagnostics).**
  Fit Upsilon(T) per transition in Burgess-Tully scaled space (low-order
  Chebyshev; a handful of coefficients per transition, smooth and
  machine-accurate where the underlying data are) and feed the generic
  n-level solver. Evaluated per cell only at output time (converged state),
  giving correct n_e-dependent emissivities and the diagnostic ratios.
  This replaces the raw-table converter idea from
  `RT_CODES/MOCASSIN/.../CHIANTI11_UPGRADE_PLAN.md` Phase A — same B-T
  descaling decisions, but the product is fitted coefficients in the
  registry files (Fortran stays table-free), each with a provenance header
  (CHIANTI version, fit range, max fit error).
- **Verification per ion:** fit vs direct CHIANTI evaluation (the existing
  `cno_cooling_comparison` notebook pattern) + PyNeb emissivity-ratio
  cross-check; both become part of the G5 gate on each ion.

**Per-ion verification gate:** emissivity ratios vs (n_e, T_e) against PyNeb
(and line fluxes vs MOCASSIN HII20/40 for ions in its abundance set).
Diagnostic pairs come first: [O III] 5007/4363, [S II] 6717/6731,
[N II] 6584/5755, [O II] 3726/3729.

**Suggested order of addition** (cooling importance + diagnostics):
O (II/III) -> N (II) -> S (II/III) -> Ne (II/III) -> C -> Ar -> Fe later.

Staging impact: the three generic pieces + O/N land inside G2 (adds ~1-2
weeks to the earlier G2 estimate); further ions are G5, each hours-to-days
(convert data, register, verify vs PyNeb). Realistic total for a validated
gas+dust+PAH HII code with the first diagnostic ion set: **~9-12 weeks**.

## 9. Non-goals (kept out deliberately)

- Metals feeding back on n_e / opacity (trace approximation) — revisit only
  if super-solar or dust-free regimes demand it.
  > **Retired by measurement, 2026-07-31; the text above is kept as the
  > founding intent.** Metal opacity has been on by default for some time
  > (removing it moves the HII40 He front 0.615 pc), and metal electrons and
  > metal photoheating were made default-on once their cost was measured:
  > −21.3 K (0.26%) at HII40, −5.8 K (0.084%) at HII20, almost all of it the
  > photoheating (`docs/METAL_COUPLING_PLAN.md`, gate G-M2). Charge exchange
  > is now two-way, so the charge a reaction moves is conserved in the
  > hydrogen-plus-metal pair. What survives of the trace approximation is the
  > *solver structure* of item 2 above — each element still an independent
  > chain evaluated on the H/He and thermal state — and a source survey found
  > that structure is what MOCASSIN, CMacIonize, TORUS and Wood all use, and
  > that even Cloudy is one matrix per element. So the founding design was
  > right about the method and wrong only about which feedbacks to leave out.
- Completeness as a prerequisite: ions are added as science needs them, not
  Z <= 30 up front.
- Dust-gas thermal coupling beyond photoheating (later, if needed).
- Velocity fields / resonance-line transfer (LaRT's domain).
- Time-dependent (non-equilibrium) ionization — equilibrium only at first.
