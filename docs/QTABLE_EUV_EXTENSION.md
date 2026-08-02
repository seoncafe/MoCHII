# Extending the SEDust Q table into the EUV — feasibility review

Date: 2026-07-21.  Status when written: review only, no implementation.
**Status now (2026-08-02): implemented and in production.**  Both the
extension and the host adoption of section 7 happened, and one thing this
memo assumes has since changed: `dust_extinction` no longer performs the
size integral at call time.  It serves the model's *precomputed*
size-integrated curve — `SEDust/data/kext_astrodust_MW_euv.dat`,
`kext_dl07_MW_euv.dat`, `kext_zubko_BARE_GR_S.dat` — interpolated onto the
model's own wavelength grid `m%lam`; the size integral is now the separate
`size_integrated_extinction`, which is what wrote those tables.  Which table a
run serves is `par%sed_kext` (MoCafe's convention: blank takes SEDust's own
standard table for the named `par%dust_model`, a named file is mandatory).  The
gate `tests/grain_kext` holds the served curve against the size integral over
the model as built, so a table that does not belong to that model is caught.
The section 7 "decide the default" step (S6) resolved to the grain model: with
`par%ion_dust_kext` blank the transport optics come from `par%dust_model`, and
the dusty smoke tests (`sedust_smoke`, `dustemis_smoke`, `pahlive_smoke`) now
run that way — astrodust, so extinction and emission are one set of grains —
rather than through the D03 `ion_dust_kext` file they carried while this was
under review.  Everything below is the record of the 2026-07-21 review and is
left as written.

Context: the SEDust_v1.00 `dust_extinction` API (commit `765eff2`) lets an
RT host take its dust opacity from the same model object, on the same
wavelength grid `m%lam`, as its emission.  MoCHII cannot use it for the
ionizing band because the model grid starts at 0.0912 um = 13.6 eV exactly.
This memo reviews what it would take to extend the astrodust Q table (and
the PAH/graphite side) below 912 A, so that one optics set serves the whole
MoCHII band — absorption, scattering, Heat_dust, femit(T), and the SEDust
stochastic heating input.

Target range: the current MoCHII band reaches ~100 eV; the hard-spectrum
plan (`docs/HARD_SPECTRUM_PLAN.md`) needs headroom to ~500 eV (24.8 A).

## 1. Verdict up front

**Feasible with data already in the tree; no new external data is needed.**
The 0.0912-um floor is a choice of the SED pipeline wavelength grid
(`DH21_wave`), not a data limitation:

| Ingredient | File in `SEDust/data/dielectric/` | Range | Covers EUV? |
|---|---|---|---|
| Astrodust refractive index | `index_DH21Ad_P0.20_0.00_1.400` | 2.485e-5 eV – 12,398 eV (1 A – 5 cm) | **yes, to 1 A** |
| Graphite Q (for the PAH continuum) | `callqcomp_D16MGemt.gz` | 1e-3 um (10 A = 1.24 keV) – 1e4 um, 41 radii, includes g | **yes, to 10 A** |
| PAH cross sections (LD01/DL07 formulas in `qpah.f90`) | analytic | defined for lambda >= 0.058 um (E <= 21.4 eV); below that `qpah` already switches to **pure graphite** | yes, by construction |

The missing piece is purely computational: Q_ext/Q_abs/Q_sca/g for the
astrodust spheroid on the extension grid, where the T-matrix solver cannot
go, plus a regeneration of the graphite Q table on the extended grid.

References verified via ADS: Draine 2003, ApJ 598, 1026 (X-ray grain
optics; the method template); Draine & Hensley 2021, ApJ 909, 94 (DH21
dielectric function, tabulated to 1 A); Hensley & Draine 2023, ApJ 948, 55
(HD23 astrodust+PAH model).

## 2. Current generation pipeline (facts from the tree)

`tmatrix/driver/run_tmatrix.f90` sweeps the fixed grids
`DH21_aeff` (169 sizes, 3.16 A – 5.01 um, 40/decade) x
`DH21_wave` (1129 wavelengths, 0.0912 – 3.98e4 um, 200/decade),
interpolating the DH21 index and dispatching on the size parameter
x = 2 pi a_eff / lambda:

- x < 0.1: `rayleigh_limit` (Draine 1992 spheroid dipole, flag 10);
- 0.1 <= x <= 50: Mishchenko random-orientation T-matrix, oblate b/a = 1.4
  (flag 0; IERR failures redirect to the nearer limit);
- x > 50: `geometric_optics_limit` — Q_ext = 2, Q_abs = 1 - exp(-4 k x),
  g = 0 (flag 20).

Flag census of the current table (190,801 rows): 47,574 T-matrix,
140,320 Rayleigh, 2,907 geometric-optics.  The GO rows sit at the largest
sizes and shortest wavelengths; their crudeness is tolerated because the
FIR-dominated SED never samples them hard.  **That tolerance does not
transfer to the EUV** (section 3).

The loader chain is already grid-flexible: `q_table.f90` takes
`na_in`/`nw_in` and reads the sizes from the file; `sed_astrodust` sets
`NLAM = qt_n_lam` from the loaded table.  Only `run_tmatrix.f90` hardcodes
`NA = 169, NW = 1129`.

## 3. Physics regimes below 912 A

Astrodust refractive index m = n + ik from the DH21 file:

| E [eV] | lambda [um] | Re(n)-1 | Im(n) | \|m-1\| |
|---|---|---|---|---|
| 10.2 | 0.1216 | +6.78e-1 | 6.38e-1 | 0.93 |
| 13.6 | 0.0912 | +3.14e-1 | 7.02e-1 | 0.77 |
| 20 | 0.0620 | -1.32e-1 | 5.48e-1 | 0.56 |
| 25 | 0.0496 | -1.43e-1 | 2.75e-1 | 0.31 |
| 54.4 | 0.0228 | -5.66e-2 | 5.39e-2 | 0.078 |
| 100 | 0.0124 | -2.48e-2 | 1.83e-2 | 0.031 |
| 300 | 0.0041 | -4.42e-3 | 1.59e-3 | 0.0047 |
| 500 | 0.0025 | -1.61e-3 | 3.11e-4 | 0.0016 |

Two regimes:

1. **13.6–~50 eV (strong-contrast EUV)**: |m-1| = 0.1–0.9, grains still
   interact strongly.  The T-matrix converges only for x <= 50, i.e.
   a <= 50 lambda / 2 pi (at 20 eV: a <= 0.49 um; at 50 eV: a <= 0.20 um).
   Larger grains need another exact method.  **Mie theory for the
   equal-volume sphere** is the practical choice: BHMIE-class codes are
   stable to x ~ 1e4 in double precision, and the tree already carries one
   (`sed/src/mie.f90`, the Draine BHMIE wrap, returns g).  The shape
   (b/a = 1.4) correction abandoned by going to spheres is quantifiable in
   the overlap window (section 5, gate ii).

2. **>~50 eV (translucent grains)**: |m-1| < 0.08 and falling as E^-2.
   This is the anomalous-diffraction / Rayleigh-Gans regime (van de Hulst
   1957; Draine 2003): scattering is a narrow forward diffraction lobe
   (g -> 1) and absorption is volume absorption.  Shape effects are
   negligible; Mie spheres are essentially exact here too.

**Why the current x > 50 GO limit must not be used in the EUV**: it
assumes an opaque grain.  At 100 eV, a = 0.3 um: x = 152, and
4 k x = 4(0.0183)(152) = 11 — still opaque, marginal; but at 500 eV the
same grain has x = 760, 4 k x = 0.95: the grain is translucent and
Q_abs = 1 - exp(-4kx) applied on top of Q_ext = 2 misassigns the balance,
while g = 0 is qualitatively wrong (the D03 kext table, which does treat
this regime correctly, has <cos> = 0.71 already at 912 A, rising toward 1
in the EUV — MoCHII's band scattering currently samples g = 0.67–0.99).

## 4. Recommended method per (a, lambda) cell of the extension

Keep the existing dispatch and add one branch:

- x < 0.1: `rayleigh_limit` (unchanged; at 500 eV this covers only
  a < ~4 A, i.e. the first grid point);
- 0.1 <= x <= 50: T-matrix (unchanged — it still converges in the
  extension zone at these x; |m| <= 1.5 there);
- **x > 50: Mie equal-volume sphere** (new branch, new flag, e.g. 30),
  replacing the GO redirect for lambda below a threshold (~0.2 um);
  the GO limit remains for the long-wavelength large-grain corner where
  it is currently used.

Optionally, a closed-form `anomalous_diffraction_limit` in
`asymptotic_optics.f90` for the |m-1| < 0.03 corner — useful mainly as an
independent cross-check on the Mie branch (ADT and Mie must agree to a few
percent there); not required for production if Mie is the branch.

Grid extension: continue `DH21_wave` at 200/decade from 0.0912 um down.
Two candidate floors:

- **24.8 A (500 eV)**: +313 wavelengths (1129 -> 1442); covers the
  hard-spectrum plan exactly.
- **10 A (1.24 keV)**: +394 wavelengths (1129 -> 1523); matches the
  graphite `callqcomp` floor so the PAH side needs no extrapolation, and
  leaves headroom above 500 eV.  **Recommended.**

Cost: +169 x 394 = 66,586 rows (+~35% table size, ~+7 MB).  Only the new
columns are computed; the existing 1129 columns are untouched (append,
then verify the overlap column is byte-identical).  The T-matrix share of
the new zone is the moderate-x sliver (fast); the Mie share is trivial.

## 5. Validation gates

i.  **Continuity at 912 A**: Q(a, 0.0912 um) from the extension side
    (T-matrix/Mie) vs the existing first column — smooth join, no step
    beyond the method difference.
ii. **Shape error bound**: Mie-sphere vs T-matrix-spheroid over the
    overlap x in [10, 50] at 0.0912–0.2 um.  This measures exactly the
    error committed by the sphere branch at x > 50; expected a few
    percent (random orientation + forward-diffraction dominance).
iii. **High-E absorption sum**: size-integrated C_abs/H at E >~ 300 eV vs
    the atomic photoabsorption sum of the depleted elements (the Draine
    2003 consistency check — at high E, grain absorption approaches the
    sum of its atoms' cross sections).
iv. **Model-level comparison**: the extended astrodust kext columns vs
    the D03 WD01 table (`par%ion_dust_kext` default) over 13.6–500 eV.
    Differences are physical (different grain model), not errors: already
    at 13.6 eV, C_abs/H = 1.43e-21 (astrodust) vs 1.83e-21 (D03), a -22%
    model shift.
v.  **Emission regression**: rerun `dustemis_smoke` — the SED above 912 A
    must be unchanged (kappB Planck integrals gain nothing below 0.1 um
    for T < 3000 K); energy closure stays 1.000 by construction.

## 6. PAH / graphite side

`qpah.f90` already contains the EUV logic: for lambda <= 1/17.25 um
= 0.058 um (E >= 21.4 eV) the PAH cross section is **pure graphite**
(currently a dead branch, since the grid floor 0.0912 um never reaches
it); the LD01 blend covers 0.0912–0.058 um with defined formulas.  What
is missing is the graphite Q table on the extended grid: the current
`qlib_gra_D16MGemt_1.400` (169 x 1009) also floors at 0.0912 um, but
`callqcomp_D16MGemt.gz` holds Draine's D16 graphite Q down to 10 A
(41 radii x 3501 wavelengths, with g).  Regenerating/resampling qlib onto
the extended grid from that file is a data operation.  (The alternative —
atomic-carbon absorption above ~100 eV — is not needed while callqcomp
covers the range.)

## 7. Pipeline and host changes after the table exists

SEDust side (all downstream of the table, mostly automatic):
- `run_tmatrix.f90`: parametrize NW, add the Mie branch + flag; the
  existing `range` mode already supports partial parallel sweeps.
- Loader (`q_table.f90`) and `sed_astrodust` (NLAM = qt_n_lam): no
  structural change — grid comes from the file.
- kappB/enthalpy integrals: physically unchanged; slight cost increase.
- `calc_kext_astrodust.f90`: regenerate `kext_astrodust_MW.dat` on the
  extended grid (it reuses sed_init, so it inherits the extension).

MoCHII side (the payoff):
- `dust_extinction` then covers the whole band: band sabs/ssca/g, the
  C_ext(lambda_ref) normalization, and femit(T) can all come from one
  model object (option-gated; D03 kext stays the default until the
  d_dusty gates are re-baselined — expect ~1% R_eff shifts from the -22%
  absorption difference at the Lyman edge, a documented model change).
- **Stochastic-heating fix rides along**: `sedust_mod` currently deposits
  the EUV J below the grid minimum into the first (13.6 eV) bin —
  energy-conserving but quantum-wrong for stochastic heating (a 54-eV
  photon delivers 4x the temperature spike of a 13.6-eV photon on a small
  grain).  With the extended grid the band tally maps onto true EUV bins
  and single-photon spikes become correct for PAHs/small grains inside
  H II regions — a genuine physics improvement independent of the
  extinction motivation.

## 8. Staged plan (upstream SEDust_v1.00 work)

- S1. Regime-map script: (a, lambda) occupancy of Rayleigh/T-matrix/Mie
      zones on the extension grid; pick the Mie/GO lambda threshold.
- S2. Mie branch in `run_tmatrix` + gate (ii) overlap comparison; ADT
      cross-check at |m-1| < 0.03.
- S3. Extend `DH21_wave` to 10 A; regenerate the new columns; gates (i),
      (iii), (iv); append and re-verify the untouched columns.
- S4. Regenerate qlib graphite from `callqcomp` on the extended grid;
      PAH branch activates automatically.
- S5. Pipeline smoke (gate v) + regenerate `kext_astrodust_MW.dat`.
- S6. MoCHII host adoption: band optics from `dust_extinction`
      (option-gated), d_dusty re-gate, then decide the default.

Open decisions for the author: extension floor (24.8 A vs 10 A —
recommend 10 A); ADT as check-only vs production branch (recommend
check-only); whether/when MoCHII's default band optics switch from D03 to
astrodust (a model choice with gate re-baselining attached).
