# Hard-Spectrum Plan: PN and AGN Fields in MoCHII

Kwang-il Seon — drafted 2026-07-18 (Claude session; reviewed plan, not yet
started).  Companion to `docs/PLAN.md`; same staged-gate discipline.
A more detailed version (verified benchmark parameters, module-level work
items, quenching-order numbers) is the "Hard-spectrum (PN and AGN)
extension plan" section of `docs/MoCHII_vs_M3.tex/pdf`.

## 1. Scope and the two tracks

"Hard spectrum" splits into two physically distinct regimes, and the plan
keeps them as separate tracks because they need different physics:

- **PN track** — thermal central stars, T_eff ~ 75–200 kK.  A 150 kK
  blackbody has kT = 12.9 eV; its flux at 300 eV (~23 kT) is negligible,
  and the lowest K-shell edge of an abundant metal (C, 290 eV) is never
  reached.  So the PN track needs **no inner-shell/Auger physics and no
  Compton terms** — it is band headroom + higher ion stages + the He II
  diffuse channels + (at PN densities) density-dependent cooling.
- **AGN track** — a power-law/big-bump continuum extending to keV X-rays.
  Above ~0.3 keV inner-shell absorption dominates the metal opacity and
  Auger multi-electron ejection couples non-adjacent ion stages; above
  ~0.5 keV Compton heating and Thomson scattering become non-negligible.
  First increment caps the band at ~500 eV (photoabsorption-dominated) and
  treats Compton as out of scope (Section 9).

Current baseline (v1.00-dev, 2026-07-18): band 13.6–100 eV (+FUV option),
registry stages C 4 / N 3 / O 3 / Ne 3 / S 4 / Ar 4 / Mg 3 / Fe 4 / Si 5 /
Cl 5 / Ca 5, te_max = 5e4 K, gbar_ff tabulated for Z_ion = 1–4, secondary
ionization (`use_sec_ion`, SvS85) available but default off, RQMC launch
available.  The He II SH95 lines, He I Porter lines, and the case-A/B
switch already cover the He II 4686 diagnostic side.

## 2. Reference assets (verified present)

| Asset | Location | Use |
|---|---|---|
| MOCASSIN PN150 + PN75 benchmarks | `~/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/benchmarks/gas/{PN150,PN75}` | primary 3D MC gates (run MOCASSIN itself CONVERGED — the HII40 lesson) |
| Cloudy c23.01 tsuite | `~/CLOUDY/c23.01/tsuite/auto/pn_paris*.in`, `nlr_paris.in`, `agn_lex00_u{0,1,m1}.in` | 1D cross-code gates: Meudon PN, AGN NLR, Lexington-2000 AGN slabs at three ionization parameters |
| Cloudy AGN SED | `~/CLOUDY/c23.01/source/parse_agn.cpp` | the standard AGN continuum parameterization (big-bump T, alpha_ox, alpha_uv, alpha_x; alpha_ox anchored at 2 keV/2500 A) for a MoCHII preset |
| Kaastra & Mewe (1993) Auger yields | `~/CLOUDY/c23.01/data/mewe_nelectron.dat` (+ `mewe_fluor.dat`) | electrons ejected per inner-shell vacancy, shell-resolved — machine-readable provenance for the Auger stage |
| Shell-resolved cross sections | Verner & Yakovlev (1995) — `sigma_vy95` already in `photo_xsec.f90` (Cl path); phfit2.f partials | inner-shell photoionization above K/L edges |
| MOCASSIN 3.x | `~/RT_Codes/MOCASSIN/MOCASSIN-3.0/` | the only 3D MC code in reach that transports hard photons: Auger-aware ion balance, Compton, `nstages` to 31, and the Ercolano et al. (2007) X-ray slab benchmark inputs. **Read section 8b first — one defect blocks quantitative use above 100 eV** |
| CHIANTI v11.0.2 | `~/RT_Codes/CHIANTI/dbase` | levels/collision strengths for the new ions (Ne V, O IV, N IV, Ar V ... all present) |
| Badnell 2023 RR/DR | `data/atomic/badnell_{rr,dr}.dat` (already in-tree) | recombination for the new stages (isoelectronic coverage is complete there) |

## 3. Stage HS1 — band headroom (small code work, mostly validation)

1. Validate `eion_max` up to 500 eV: VFKY96 fits and the exact hydrogenic
   forms are valid far beyond that; the work is checking the aligned-edge
   builder, the FUV+ionizing two-segment bookkeeping, and the SvS85
   integrals with many high-energy bins.  Recommend `nnu_ion >= 64` when
   the band is widened (more thresholds to align; the largest-remainder
   allocator already degrades gracefully).
2. Extend `gbar_ff` from Z_ion 1–4 to 1–8 (one parameter, NGZ; the vh14
   table already covers the gamma^2 range).
3. `te_max` as a hard-field option: raise to 1e5–1e6 K for X-ray-heated
   zones; the Tier-1 cooling fits are fitted over 1e3–1e5 K, so any ion
   used with a higher te_max needs its fit range extended by the pipeline
   (a regeneration flag, not new code).
4. Guidance defaults for hard fields: `use_sec_ion = .true.` (plentiful
   E0 > 40 eV photoelectrons; exactly the SvS85 regime) and
   `launch_sequence = 'sobol'` (the partially neutral shielded zones where
   secondaries matter are also where launch noise hurts most).

**Gate HS1**: G0-style fixed-state analytic attenuation with a 150 kK
blackbody on the widened band — Gamma_HI/HeI/HeII vs the analytic
integrals; Q(>13.6)/Q(>24.6)/Q(>54.4) ratios vs exact Planck integrals to
<0.1% with aligned edges.

## 4. Stage HS2 — registry upward in ion stage (the data operation)

Target stages, thresholds, and the diagnostics they unlock:

| Element | now → target | new thresholds [eV] | key new lines |
|---|---|---|---|
| C | 4 → 5 | C IV→V 64.5 | C IV 1548/1550 (resonance — flag optically thin caveat) |
| N | 3 → 5 | N III→IV 47.4, IV→V 77.5 | N IV] 1486, N V 1240 |
| O | 3 → 6 | O III→IV 54.9, IV→V 77.4, V→VI 113.9 | [O IV] 25.9 um, O IV] 1400, O V] 1218 |
| Ne | 3 → 6 | Ne III→IV 63.4, IV→V 97.1, V→VI 126.2 | [Ne IV] 2422/2425, [Ne V] 3346/3426 + 14.3/24.3 um, [Ne VI] 7.65 um |
| S | 4 → 6 | S IV→V 47.2, V→VI 72.6 | [S IV] 10.5 um (needs the S IV n-level file), [S V]? UV |
| Ar | 4 → 6 | Ar IV→V 59.8, V→VI 75.0 | [Ar V] 7.90/13.1 um, [Ar VI] 4.53 um |
| Mg/Fe/Si/Cl/Ca | hold | — | extend only on diagnostic demand (memo rule); coronal [Fe VII+] deferred (Section 9) |

Work items, all through the existing pipeline (`tools/fitting/`):
1. Element files: VFKY96 rows exist for every stage; Voronov CI all
   stages; Badnell RR+DR cover these sequences (verified isoelectronic
   coverage); **charge exchange**: Kingdon & Ferland fits stop at stage
   ~4 — for stages >= 5 adopt the standard fast-rate scaling used by
   Cloudy (Dalgarno-type ~1.92e-9 z cm^3/s recombination-only) with an
   explicit provenance note, or document neglect; decide once, apply
   uniformly.
2. Tier-1 cooling fits for each new stage (fit range matched to the
   chosen te_max).
3. Tier-2 n-level files for the diagnostic ions above (CHIANTI v11 has
   all of them; the 40-level/600-transition caps are already in place).
4. PN abundance set: the MOCASSIN PN150/PN75 input abundances as a
   documented example file (they differ from the HII20 defaults).

**Gate HS2**: PyNeb line-ratio checks ion by ion ([Ne V] 3426/24.3 um
temperature pair, [Ne IV] and [Ar IV] density pairs, [O IV] 25.9 um) at
the established database-spread tolerance; end-to-end smoke showing the
new lines in `_lines.txt` with sane ratios.

## 5. Stage HS3 — He II excited-channel diffuse photons

The PN He III zone is large, and case-B He II recombinations emit
H/He-ionizing photons beyond the ground continuum already carried:
He II Ly-alpha (40.8 eV — ionizes H I AND He I), the He II Balmer
continuum (>13.6 eV), and the He II two-photon continuum (partly
ionizing).  Add a fifth diffuse channel following the `hei_diffuse`
pattern exactly (branching ratios from the same AGN3/Osterbrock tables;
energies sampled analogously).  This is the same "correct energy
bookkeeping" class as `hei_diffuse` — measurable in a PN model where
He III fills much of the nebula, unlike HII40 where it moved Te by ~1%.

**Gate HS3**: energy conservation of the new channel; PN150 He I/He II
structure and Te shift recorded with the channel on/off.

## 6. Stage HS4 — PN validation gates

1. **MOCASSIN PN150 and PN75**, run CONVERGED on both sides (the HII40
   re-run lesson: the stored reduced baselines are not equilibrium
   states).  Compare Te(H II), the He I/II/III structure, ionization
   fractions of the new high stages, and the line table (He II 4686/Hbeta,
   [Ne V]/[Ne III], [O IV] 25.9 um, [Ar V]).
2. **Cloudy `pn_paris`** (same SED, abundances, density, dust off): 1D
   sphere vs a MoCHII 3D sphere — compare zone-by-zone Te and ion
   fractions along the radius plus ~20 integrated line ratios.  This
   separates atomic-data differences from transport differences, with the
   Lexington-published spread as the acceptance band.

**Acceptance**: land inside (or bracket) the published Lexington PN150
code spread the way HII20/HII40 do today; attribute any residual with the
HII40 methodology (binning, diffuse treatment, data vintage — one lever
at a time).

## 7. Stage HS5 — AGN SED input

1. `ion_spectrum = 'agn'` preset with the Cloudy parameterization
   (parse_agn.cpp): big-bump temperature, alpha_ox (anchored 2500 A–2 keV),
   alpha_uv, alpha_x, and the IR/X-ray cutoffs; normalized by the existing
   luminosity/derive-or-rescale rules.  (The tabulated-file route via
   `spectrum_type='per_ev'` already works today and stays the reference
   path; the preset is convenience + reproducibility.)
2. An ionization-parameter helper in the log output (U computed from
   Q_ion, n_H, and the source distance) so AGN-slab setups can be matched
   to Cloudy inputs directly.

**Gate HS5**: bin-integrated preset spectrum vs Cloudy's `save continuum`
output for the same parameters (shape to ~1%); Q ratios exact.

## 8. Stage HS6 — inner-shell absorption + Auger (AGN track only)

The one genuine framework extension:

1. **Shell-resolved opacity**: above each K/L edge add the inner-shell
   photoionization cross sections (VY95 partials; the `sigma_vy95`
   path generalizes — the outer-shell VFKY96 total stays below the first
   inner edge).
2. **Auger yields**: Kaastra & Mewe (1993) electrons-per-vacancy from the
   Cloudy data file (shell-resolved, machine-readable) — an inner-shell
   ionization of stage i lands at stage i + n_e(shell), not i + 1.
3. **Cascade generalization**: `species_mod` currently solves a strict
   adjacent-stage chain (product form).  Auger couples non-adjacent
   stages, so the cascade becomes a (small) transition-matrix solve for
   each element — reuse the n-level partial-pivot solver on the stage
   populations.  Keep the chain path as the default when no inner-shell
   channel is active (bit-identical off-path, the house rule).
4. Auger photoelectrons are fast: they enter the existing SvS85 secondary
   partition automatically through the per-absorber excess energy — no
   new heating physics needed at this stage.

**Gate HS6**: X-ray-illuminated slab vs Cloudy `agn_lex00_u*` (three
ionization parameters): ion-fraction profiles of C/N/O/Ne through the
partially ionized zone, Te profile, and the He II/H recombination lines.
Acceptance = the Lexington-2000 published cross-code spread.

## 8b. What MOCASSIN 3.x does about hard spectra, and what to take from it

Read 2026-07-31 in `~/RT_Codes/MOCASSIN/MOCASSIN-3.0/`
(`github.com/mocassin/MOCASSIN-3.0`).  It is the closest thing to a worked
example of HS6 in a 3D Monte Carlo code, and worth reading before implementing.
It is **not** usable as a quantitative reference without patching, for the
reason in 8b.5.

**On which repository to read.**  `github.com/rwesson/mocassin_xray` (one
commit, 2017-09-18, tagged `mocassin.3.04.00`) was downloaded and compared file
by file.  It is **the same code, not a fork**: of its 23 Fortran sources, 22 are
byte-identical to MOCASSIN-3.0, and the two that differ, `photon_mod.f90` and
`update_mod.f90`, differ by exactly one line each -- `use inf_nan_detection`.
MOCASSIN-3.0 adds `infnan.f90` to supply portably the `isnan` that both trees
already call (four cooling guards in `update_mod`, one `costheta` guard in
`photon_mod`); a portability fix, not physics.  Its only later commits, in 2021,
fold and then unfold the same two continuation lines in `iteration_mod.f90`, so
there is **no physics change between 2017 and 2021**.  The older snapshot was
deleted; read MOCASSIN-3.0.

That matters for how the rest of this section reads: what follows describes
**mainline MOCASSIN 3.x**, not one author's experimental branch.  The design
choices and the defects below are what the reference code has shipped for eight
years.  Line numbers are MOCASSIN-3.0's; they run one higher than the 2017
snapshot's in `update_mod.f90` and `photon_mod.f90`, and are identical
elsewhere.

### 8b.1 What it implements that stock does not

- **Inner-shell cross sections**, Verner & Yakovlev (1995) for every inner
  shell and VFKY96 for the outer, up to 7 subshells per ion, each with its own
  energy window carried to `kShellLimit = 7.35e4` Ryd
  (`source/ph_mod.f90:586-664`, `source/hydro_mod.f90:238-309`,
  `source/constants_mod.f90:42`).  Every subshell of every ion enters the cell
  opacity (`source/ionization_mod.f90:429-437`), so packets are genuinely
  attenuated by inner shells.  This is the same data path HS6 item 1 plans
  (`sigma_vy95` is already in `photo_xsec.f90`), which is a useful confirmation
  that the plan's choice of source is the one a working code made.
- **Auger yields** from Kaastra & Mewe (1993), `data/auger.dat` (Z = 3-30, 1696
  records, 351 with multiplicity >= 3, neutral Fe K peaking at 6-8 ejected
  electrons), read by `makeAugerData` (`source/hydro_mod.f90:35-53`).  Note this
  is the same determination HS6 item 2 already names from the Cloudy data file,
  so the two references agree on provenance.
- **Compton scattering wired into transport**: Klein-Nishina cross sections and
  an angle CDF (`ph_mod.f90:1883-1922`), added to the cell opacity
  (`ionization_mod.f90:405-410`), packets scattered and redshifted
  (`photon_mod.f90:2509-2551`, `:6806-6845`), with Tarter exchange plus recoil in
  the thermal balance (`update_mod.f90:3835-3898`).  On by default.
- **`nstages` raised** from 7 to a default of 17 and a ceiling of 31
  (`set_input_mod.f90:75`, `:165-174`), so fully stripped Fe is representable.
- **Fluorescence as transported packets**, 14 lines including Fe K-alpha at
  6.4 keV, split into 20 sub-packets at the absorption site and traced
  (`photon_mod.f90:2652-2686`), with per-stage Fe yields from Krolik & Kallman
  (1987) (`photon_mod.f90:4123-4138`).  A separate warm-start driver
  `mocassinFluorescence` runs it.

### 8b.2 The lesson for HS6 item 3: it did NOT solve the stage matrix

This is the most useful thing in MOCASSIN 3.x, and it is a negative result.  Faced
with exactly the problem HS6 item 3 identifies -- Auger moves a stage by more
than one, which an adjacent-stage chain cannot represent -- MOCASSIN **kept the
chain** and patched the flux.  The live routine is `ionBalance`
(`source/update_mod.f90:1576`, called at `:318`); `ionBalance2` at `:1265` and
`ionBalance3` at `:1995` are dead code, which is worth knowing before quoting
line numbers out of this file.  The Auger block of the live routine is
`:1840-1856`:

```
maxim = min( ion+nelec-1, min(nstages-1,elem))
if (grid%ionDen(cellP,elementXref(elem),maxim) > 1e-30 ) then
   ratio = ionDen(...,ion) / ionDen(...,maxim)
else
   ratio = 1.
endif
out(maxim) = out(maxim) + photoIon(1,nshell,elem,ion)*auger(elem,ion,nshell,nelec)*ratio
```

and the stage populations then come from an explicitly **bidiagonal** solve,
`mat(ion) = mat(ion-1)*out(ion-1)/in(ion-1)`, under the comment "invert and
solve bidiagonal matrix" (`:1869-1876`) -- adjacent stages only, by
construction.

The jump rate is charged to the *destination* stage, rescaled by the population
ratio, so that the net flux across the last boundary crossed comes out right
once `ratio` converges in the outer loop.  Two consequences follow directly, and
both argue for doing HS6 item 3 properly rather than copying this:

1. **Only the final boundary is fed.**  For `nelec >= 3` the intermediate
   boundaries `ion+1 ... ion+nelec-2` never receive the flux, so triple-and-higher
   Auger cascades underfill the high stages -- and those are exactly the Fe, Si
   and S K-shell cascades MOCASSIN's own data file says dominate.
2. **`ratio = 1.` when the destination is empty** (`update_mod.f90:1850`),
   which is the regime on the first iterations.

MoCHII's plan -- reuse the n-level partial-pivot solver on the stage populations
-- is the right call, and Cloudy independently confirms it: its ionization block
fills arbitrary super-diagonal entries precisely for Auger
(`ion_solver.cpp:727-736`, `:661-668`), and its own legacy nearest-neighbor chain
solver is dead code (`docs/METAL_COUPLING_PLAN.md` section 6).  **Two codes, one
conclusion: the chain has to go when Auger arrives.**

### 8b.3 Energy grid: a warning, not a model

MOCASSIN 3.x adds a separate grid constructor above `nuMax > 25` Ryd
(`source/grid_mod.f90:391-438`) reaching 7.3e6 Ryd, but it **abandons stock's
threshold-insertion scheme in doing so**: stock inserts every ionization
threshold as a triplet before log-filling (`grid_mod.f90:264-296`), and above
25 Ryd 3.x inserts nothing, so a K edge lands wherever a 0.008-dex (~1.9%)
mesh puts it.  Its own benchmarks compensate with `nbins 10000`
(`benchmarks/gas/X1/input/input.in`).  For HS1, keep MoCHII's edge-resolved
construction when the band is widened; brute-force bin counts are not the fix.

### 8b.4 Secondary ionization: heating only

MOCASSIN 3.x applies the Xu & McCray (1991) heating efficiency
`heatef = 0.9971*(1-(1-x_e^0.2663)^1.3163)` above 100 eV
(`update_mod.f90:3619-3627`, `:4023-4042`) but **never adds a secondary
ionization rate to `gamma`**.  MoCHII already has the Shull & van Steenberg
(1985) partition with both the heating fraction and the secondary ionizations
(`par%use_sec_ion`), so on this axis MoCHII is ahead, and HS6 item 4's claim
that Auger photoelectrons enter the existing partition automatically is sound.

### 8b.5 Defects that block quantitative use

Verified by reading, not inferred:

- **`update_mod.f90:4031-4032` uses the threshold bin for the whole X-ray
  band.**  In the secondary-ionization loop `do i = ilow, iup` the integrand
  reads `radField(IPnuP)`, where the low-energy loop twenty lines above
  correctly reads `radField(i)` (`:4006`).  Same structure, same names, one
  index different.  Every heavy-ion photoionization rate and photoheating rate
  above 100 eV is therefore computed from the threshold-bin mean intensity
  replicated across the entire band.  For a hard-spectrum calculation this is
  the dominant error path, and it means **MOCASSIN 3.x cannot be used as a
  quantitative gate above 100 eV unless this is patched first**.  It stands
  identically in the 2017 snapshot and the 2021 head, so it has been there at
  least eight years and will not be fixed upstream on our schedule.
- `ph_mod.f90:314` builds the Al K-alpha band from carbon's pointer
  (`elementP(6,...)`), so that line is emitted over the wrong band.
- The Compton scatter direction (`photon_mod.f90:6825-6832`) applies the
  KN-sampled angle to two components only rather than rotating about the
  incoming direction, which randomizes the phase function.
- `augerE` is built from `auger(...,1)`, a *probability*, used as a multiplier
  on the threshold energy to form a work function (`update_mod.f90:3990`).

### 8b.6 Benchmarks: setups yes, reference values no

`benchmarks/gas/` carries **X01, X1, X10** -- the Ercolano, Young, Drake &
Raymond (2007, ApJS) X-ray-illuminated constant-density slabs at ionization
parameter 0.1, 1 and 10, with a 5-point tabulated hard spectrum, 3 x 30 x 3
slab, `nbins 10000`, `nuMax 1.e6`, `nstages 31` -- plus an AGN **NLR** slab
(`contShape powerlaw 1.3`, `nstages 10`).  These are directly relevant to HS5
and HS6 and are worth adding to the gate list.  **But the `output/` directories
that `docs/README.3.00.tex:355-380` promises, holding the reference line fluxes,
are absent from the repository**: you get the model setup, not the published
numbers, so the comparison has to go through the ApJS tables.

### 8b.7 It is not an upgrade to the local 2.02.73.2 tree outside X-rays

Despite the higher version number, MOCASSIN 3.x is the narrower tree.  Comparing
Fortran sources: **no `.f90` exists in 3.x that is absent from the local
`mocassin-mocassin.2.02.73.2`**, while the local tree has fifteen that 3.x lacks
(`dda_mod`, `scattering_mod`, `echo_mod`, `multigrid_mod`, `reflection_mod`,
`shared_window_mod`, `grid_v2_mod`, `plane_ionization_mod`, `view_angles_mod`,
`subgrid_mod`, `mpi_utils_mod`, `readdata_mod`, `data_path_mod`, `gaunt.f90`,
`test_dda`).  Both already carry `fluorescence_mod` and an identical
`data/auger.dat`; what 3.x adds on top is Compton **wired into transport**
(the local tree has the Klein-Nishina routines in `ph_mod` but `lgCompton` is
declared and read at one site only), the raised `nstages`, the >25 Ryd grid, and
the Auger flux trick of 8b.2.  So: **mine it for the hard-spectrum pieces, and
do not treat it as a newer MOCASSIN for anything else.**

## 9. Stage HS7 — density-dependent metal cooling in the thermal loop

Required for PN densities (n_e ~ 1e3–1e6 cm^-3): the IR fine-structure
lines that dominate low-Te cooling have n_crit ~ 1e2–1e5 cm^-3 and quench
first, shifting the thermal balance toward the optical forbidden lines —
the low-density Tier-1 fits overestimate cooling there.  (This is the
comparison memo's top-value item and is PN-critical, so it belongs to the
PN track even though it is independent of the band work.)

Design (keeps "fitted offline, evaluated online"):
- **Primary**: extend Tier-1 to Lambda(T, n_e) tables — the pipeline
  already solves the n-level system; add a log-n_e grid (~6–8 nodes,
  1e1–1e7 cm^-3) and bilinear-in-log evaluation in `species_mod`.  The
  low-density node must reproduce today's fits (regression).
- **Verification**: an in-loop n-level solve for the dominant coolants on
  a small test, confirming the table interpolation error < 1% in Lambda.

**Gate HS7**: two-zone analytic check (low vs high n_e) against the exact
n-level cooling; PN150 Te with tables on/off; HII20/HII40 unchanged to
<0.1% (their densities sit in the low-density limit).

## 10. Deferred / out of scope (record once, revisit on demand)

- Compton heating/cooling and Thomson-scattering transport (> ~0.5 keV;
  Cloudy is the reference implementation).
- Fe K fluorescence; coronal ions beyond ~stage 8 ([Fe VII] 6087,
  [Fe X] 6374 need Fe 8–11 stages and a hotter thermal range).
- X-ray resonance-line transfer; dust sputtering/grain destruction in
  hard radiation fields (the PAH survival switches are a partial hook).
- Time dependence (finite source age) — MAPPINGS remains the 1D
  reference.

## 11. Suggested order and effort

| Order | Stage | Framework risk | Blocking? |
|---|---|---|---|
| 1 | HS1 band headroom | low (validation-heavy) | gates HS2+ |
| 2 | HS2 registry upward | none (data operation) | gates HS4 |
| 3 | HS3 He II diffuse channels | low (existing pattern) | sharpens HS4 |
| 4 | HS4 PN gates (MOCASSIN + Cloudy) | none | PN track done |
| 5 | HS7 density-dependent cooling | medium (pipeline + eval) | PN accuracy |
| 6 | HS5 AGN SED preset | low | gates HS6 |
| 7 | HS6 inner-shell/Auger | **high** (cascade matrix) | AGN track done |

HS7 can run in parallel with HS2–HS4 (independent code paths); doing the
first PN150 pass BEFORE HS7 is deliberate — it measures how much of the
PN residual the density-dependent cooling must explain, the same
one-lever-at-a-time attribution that worked for HII40.
