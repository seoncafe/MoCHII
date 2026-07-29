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
