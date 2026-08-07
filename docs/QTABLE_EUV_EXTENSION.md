# Extending the SEDust Q table into the EUV — feasibility review

Date: 2026-07-21.  Status when written: review only, no implementation.
**The extension was built upstream and MoCHII now runs on it** — with one of the
review's recommendations left out, recorded below.  Sections 1–8 are the
2026-07-21 review, left as written; this preamble records what was actually
built and what it did to MoCHII, in the order it happened: the table itself
(2026-08-04), then the data layout and host API it is delivered through
(2026-08-07).

## What was executed (upstream SEDust, adopted 2026-08-04)

The Q table was regenerated on a wavelength grid extended past the review's
recommended 10 A floor, all the way to the DH21 dielectric function's own end:

| | before | after |
|---|---|---|
| `q_astrodust_P0.20_Fe0.00_1.400.dat` | 169 x 1129, 0.0912–3.981e4 um | base SEDust table |
| `q_astrodust_P0.20_Fe0.00_1.400_euv.dat` | — | 169 x 1762, **1.0e-4**–3.981e4 um |
| file size | 21 MB | 33 MB (EUV companion) |

MoCHII took the `_euv.dat` companion as its default Q table, because MoCHII is a
photoionization code and the unsuffixed 1129-point file is not an adequate
default for one.  That choice is no longer expressed as a file name: since
2026-08-07 the two are one array in the model's HDF5 product and `include_euv`
picks the view (see below).

The short end is therefore 12.4 keV, not the 1.24 keV of section 4's
recommended floor.  The wavelength axis is a separate grid file
(`DH21_wave_to_12keV`) and the existing columns were appended to, not
recomputed: the flag census of section 2 for lambda >= 0.0912 um is reproduced
row for row by the new table (47,574 T-matrix / 140,320 Rayleigh / 2,907
geometric-optics), and the kext table those columns feed is unchanged to
rounding there (worst 1.4e-10 on C_ext).

**One recommendation was not adopted.**  Section 4 proposed a fourth branch —
Mie for the equal-volume sphere at x > 50 in the EUV — because the
geometric-optics limit assumes an opaque grain and returns g = 0 (section 3).
That branch does not exist: `spheroid_optics.f90` still dispatches on
`X_RAYLEIGH_MAX = 0.1` / `X_GEOMETRIC_MIN = 50` alone, no flag was added, and
the table's flag column below 0.0912 um reads 46,008 T-matrix, 59,834
geometric-optics, 1,135 Rayleigh (measured here on the table in the tree).  The
upstream table header states this openly — "where the size parameter leaves the
T-matrix range the Rayleigh dipole limit (x < 0.1) and the geometric optics
limit (x > 50) are used, exactly as the Q table itself does" — so it is a known
approximation, not an oversight, but section 3's objection to g = 0 for a
translucent grain stands and is **open**.

## What it did to MoCHII

1. **The extension is gone.**  `grain_model_mod` still asked for the grid down
   to `1.23984/par%eion_max`, but SEDust prepends points only when that is
   shorter than the table's first wavelength.  At the default
   `par%eion_max` = 100 eV that is 0.0124 um against the table's 1.0e-4 um, so
   nothing was prepended; astrodust and DL07 came back with NLAM = 1762 and
   lambda(1) = 1.0e-4 um at 100 eV, at 150 eV, and with no `lam_min` at all —
   the same grid all three times.  An extension appeared only above
   `par%eion_max` = 12398.4 eV, and astrodust then *refused* it (build status 6:
   the band's optics are the same spheroid and this tree carries the Q table
   without the T-matrix), which `build_grain_model` turned into a message naming
   `par%eion_max`.  SEDust does not silently substitute an equal-volume sphere.
   Since 2026-08-07 `grain_model_mod` does not ask at all: no `lam_min` is
   passed, so neither the prepend nor the refusal is reachable from MoCHII, and
   the grid is the stored axis in every run.  The reasons are below.
2. **The ionizing-band grains changed identity.**  The old kext table's header
   recorded its own EUV as "Mie for the VOLUME-EQUIVALENT SPHERE on the same
   Draine & Hensley (2021) dielectric function"; the new one is the same
   random-orientation T-matrix on the same b/a = 1.4 oblate spheroid as the rest
   of the table.  That, not the disappearance of a code path, is the physical
   content of this update: a shape approximation was replaced by the shape the
   model is defined on.
3. **How much it moved.**  Above 0.0912 um, nothing beyond rounding.  Below it
   the curve was recomputed (590 -> 633 grid points).  Over the band MoCHII
   transports, 13.598–100 eV (137 points), the astrodust table differs from the
   old one by at most 1.010e-2 (C_ext), 3.798e-2 (C_abs), 4.645e-2 (C_sca), with
   medians 5.44e-3 / 1.11e-2 / 9.48e-3, and the band integrals move **+0.459%**
   (C_ext) and **+0.910%** (C_abs).  Much larger differences, up to 1300x, sit
   in 1.0e-4–1.1e-4 um, i.e. at the 12.4 keV end — **outside** anything MoCHII
   transports.  The two must not be quoted together.
4. **The gate tightened.**  With no extension there is no interpolation:
   every model wavelength is a table wavelength, so `tests/grain_kext` measures
   rounding on both sides of 0.0912 um and its EUV tolerance went 3e-2 -> 1e-9.
   (Held for astrodust; DL07 and Zubko moved to 1e-6 on 2026-08-07, for a
   reason that is a property of the products rather than of the grid — see
   below and `docs/SEDUST_API_AUDIT.md`.)
5. **Cost.**  Building the 1.56x longer grid took 156.4 s (astrodust), about
   73 s (DL07) and 0.55 s (Zubko) per rank at `par%sed_NT` = 200, paid at
   startup.  The corresponding figure on the 1129-point grid was not measured
   before the table was replaced.  Re-measured on the HDF5 products below.

## How the table is delivered now (2026-08-07)

SEDust reorganized what it ships and how a host asks for it, and MoCHII follows.
The table of section 3 did not change; the way MoCHII reaches it did.

**One file per model, and it is HDF5.**  `data/` is now one directory per dust
model — `astrodust/`, `dl07/`, `zubko/` — beside the two the models share,
`dielectric/` (optical constants, PAH cross sections) and `release/` (published
tables and the size distribution).  A model's optics are one product,
`data/<model>/sedust_<model>.h5`, holding its wavelength axis (`/grid/lambda`),
its `(lambda, a_eff)` cross-section tables (`/qtable`) and its size-integrated
extinction curve (`/kext`) together, with the index `i_lyman` at which the axis
crosses the Lyman limit.  The EUV and non-EUV forms are therefore no longer two
files but one array and a cut, and `include_euv` picks the view.

**One host call, and one directory.**  `build_dust(m, model, data_dir, ...)`
resolves everything the model is made of from `data_dir` — the optics product,
the size distribution, the ZDA config, and the dielectric functions the optics
were computed on — so `par%sed_qtable`, `par%sed_sizedist`,
`par%sed_zubko_config` and `par%sed_zubko_dir` are replaced by a single
`par%sed_data_dir`, blank by default and resolved to `<executable
dir>/SEDust/data` at run time.  `par%sed_workdir` is gone with them: SEDust used
to open its dielectric functions on paths hard-coded relative to `sed/` with no
`iostat`, so the host had to `chdir` into a SEDust `sed/` tree or the process
died before a status existed; those opens now go through the same data root and
carry `iostat`, and a bad directory comes back as a status
(`docs/SEDUST_API_AUDIT.md`, Finding 1).  `par%sed_kext` stays as the override;
blank now means `/kext` of the model's own product.  `grain_model_mod` went from
eight build calls to two.

**`include_euv = .true.`, and no `lam_min`.**  Three reasons, all measured or
read out of SEDust: (a) with the whole axis the astrodust grid already starts at
1.0e-4 um = 12398.4 eV, so `lam_min = 1.23984/par%eion_max` does nothing below
that; (b) `lam_min` is a coverage requirement rather than a grid — a model that
cannot reach it refuses to build, and Zubko's grid *is* its own optics table, so
a floor passed uniformly would take that model out rather than widen it;
(c) whether the grid spans the transported band is checked with the numbers in
`gas_opacity_mod` at setup, for every model alike.  This is what retires item 1
above.

**Measured grids**, with the arguments `grain_model_mod` passes:

| model | NLAM | lambda(1) [um] | lambda(N) [um] |
|---|---:|---|---|
| astrodust | 1762 | 1.0000e-04 | 3.9810e+04 |
| dl07 | 1823 | 6.2054e-05 | 3.9810e+04 |
| zubko | 1201 | 1.0000e-03 | 1.0000e+04 |

DL07's is its **own** axis now.  It used to be handed the astrodust Q table to
borrow a grid from, and came back on astrodust's 1762 points from 1.0e-4 um; its
product carries its own 1823 points from 6.2054e-05 um, where the D03 optical
constants its optics are Mie on actually end.

**The optics did not move.**  astrodust built the old way (`build_astrodust` on
the text `_euv` Q table) against the new way (`build_dust` on the product):
the wavelength axis is identical (relative difference 0.0) and the served
extinction agrees to 4.91e-13 (C_ext), 4.85e-13 (C_abs), 4.72e-13 (C_sca) and
5.0e-14 on `<cos>` — the precision the text file was written with, so the two
are the same numbers.

**The text products are gone from this tree.**  SEDust writes the same optics as
text (`q_*.dat` from `make_qtable.x`, `kext_*.dat` from `calc_kext.x`); nothing
in a MoCHII run opens one, so the 257 MB they occupied was removed and
`SEDust/populate_data.sh` no longer copies them.  Everything else in a model
directory is an input rather than a product — Zubko's ZDA config, size
distributions, DustEM optics and calorimetry; astrodust's DH21 tabulations — and
stays.

**Cost, re-measured** on the products, single rank at `par%sed_NT` = 200 (the
build is not OpenMP-parallel; forcing one thread against 72 gives the same
numbers to 0.05 s):

| model | 2026-08-04 (text tables) | now (HDF5 product) |
|---|---:|---:|
| astrodust | 156.4 s | 152.9 s |
| dl07 | ~73 s | 2.4 s |
| zubko | 0.55 s | 0.18 s |

astrodust is unchanged, which says the cost is the temperature-grid work
(enthalpy and the Planck-averaged opacity over `par%sed_NT` points), not reading
the table.  DL07's 30x is real work removed: given the astrodust Q table it had
to solve Mie on the D03 dielectric functions for each of its four populations,
and it now reads its own stored `/qtable`.

**One gate tolerance.**  `tests/grain_kext` holds all three models to
`TOL_KEXT_NODE = 1e-12`, worst measured 3.3e-16 (astrodust), 2.7e-15 (DL07),
7.2e-15 (Zubko) — rounding on a double.  It was two tolerances, 1e-9 for
astrodust and 1e-6 for the other two, because `calc_kext.x` wrote their `/kext`
from the 7-digit text optics while `/qtable` of the same file carried full
double precision; both halves of every product now come from one set of
numbers.  `docs/SEDUST_API_AUDIT.md` records that and four other properties of
the v1.00 API a host used to have to work around, all since resolved.

## Zubko optics (2026-08-07)

Zubko's product now stores two sets of optics and the caller picks:
`zubko_optics = 'zda'` (default) is the tables the Camps et al. (2015) benchmark
distributes, which is what the seven codes it compares against read;
`'mie_d03'` is SEDust's own recomputation, Mie on the Draine (2003) optical
constants, under `/qtable/*_mie_d03` with `/kext_mie_d03` beside it.  MoCHII
takes the default, the same rule a blank `par%sed_kext` follows, and adds no
input for it.

That default changes which grains `par%dust_model = 'zubko'` means.  Against the
curve this tree served before, over the transported band 13.598–100 eV (148
points): worst single wavelength 1.96e-3 (C_ext), 1.64e-3 (C_abs), 3.80e-3
(C_sca), band integrals +0.117% / +0.082% / +0.244%, albedo 5.1e-4 and
&lt;cos&gt; 2.1e-3 absolute, `C_ext`(0.55 um) +0.027%.  Longward of 1 um it
reaches 7.7% at 6958 um, which is where the distributed graphite table and the
D03 Mie recomputation differ by 8% in `Q_abs`.  No MoCHII test or example input
names `zubko` — `par%dust_model` defaults to `astrodust`, whose product did not
change — so nothing recorded moves; `tests/grain_kext` is the only case that
builds it.

Part of the difference is the PAH population.  Measured on the distributed
tables: their PAH `Q_sca` and &lt;cos&gt; are identical to graphite's at the same
radius, every one of 28 radii x 1201 wavelengths, and their `Q_abs` is identical
to graphite at the 303 wavelengths shortward of 0.0578 um (21.2 eV, the DL07
PAH-to-graphite transition).  The recomputation used to store no scattering for
the PAHs at all, so it now carries theirs as well — `C_sca` up to 3.7e-3 higher,
`C_abs` unchanged to the last bit.

## Earlier adoption note (2026-08-02)

One thing this memo assumes changed before the table did: `dust_extinction` no
longer performs the size integral at call time.  It serves the model's
*precomputed* size-integrated curve — then the text tables
`SEDust/data/kext_astrodust_MW_euv.dat`, `kext_dl07_MW_euv.dat`,
`kext_zubko_BARE_GR_S.dat`, since 2026-08-07 `/kext` of the model's own product —
interpolated onto the model's own wavelength grid `m%lam`; the size integral is
now the separate `size_integrated_extinction`, which is what wrote those curves.
Which curve a run serves is `par%sed_kext` (MoCafe's convention: blank takes
SEDust's own standard curve for the named `par%dust_model`, a named file is
mandatory).  The
gate `tests/grain_kext` holds the served curve against the size integral over
the model as built, so a table that does not belong to that model is caught.
The section 7 "decide the default" step (S6) resolved to the grain model: with
`par%ion_dust_kext` blank the transport optics come from `par%dust_model`, and
the dusty smoke tests (`sedust_smoke`, `dustemis_smoke`, `pahlive_smoke`) now
run that way — astrodust, so extinction and emission are one set of grains —
rather than through the D03 `ion_dust_kext` file they carried while this was
under review.

---

Everything below is the 2026-07-21 review as written.

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

Outcome added 2026-08-04; see the preamble for the numbers.

- S1. Regime-map script: (a, lambda) occupancy of Rayleigh/T-matrix/Mie
      zones on the extension grid; pick the Mie/GO lambda threshold.
      — *not verifiable from this tree*; the MoCHII copy of SEDust carries
      `tmatrix/output/` only.
- S2. Mie branch in `run_tmatrix` + gate (ii) overlap comparison; ADT
      cross-check at |m-1| < 0.03.
      — **not done.**  `spheroid_optics.f90` still dispatches on
      `X_RAYLEIGH_MAX = 0.1` / `X_GEOMETRIC_MIN = 50` with no fourth branch and
      no new flag, so 59,834 of the 106,977 EUV rows are the geometric-optics
      limit with g = 0.  Section 3's objection is therefore still open.
- S3. Extend `DH21_wave` to 10 A; regenerate the new columns; gates (i),
      (iii), (iv); append and re-verify the untouched columns.
      — **done, and past the recommended floor**: the new axis
      `DH21_wave_to_12keV` runs to 1.0e-4 um (12.4 keV), 1762 points.  It is not
      the proposed uniform 200/decade continuation: points 0–631 follow the DH21
      dielectric function's own energy axis so its absorption edges are resolved
      as tabulated, point 632 is an interpolated 0.0912*(1-1e-4) placed to
      resolve the Lyman-limit step of the *radiation field*, and 633–1761 are the
      original `DH21_wave`.  The append-and-verify discipline held: the old
      columns reproduce row for row.
- S4. Regenerate qlib graphite from `callqcomp` on the extended grid;
      PAH branch activates automatically.
      — **superseded, not done as written.**  `qlib_gra_D16MGemt_1.400` is
      still 169 x 1009 with the 0.0912 um floor.  The carbonaceous absorption
      above 21.4 eV, where the DL07 PAH cross section is zero by construction,
      instead comes from D03 graphite by Mie on the Draine 2003 dielectric
      functions (1/3 parallel + 2/3 perpendicular), as the kext table header
      records.
- S5. Pipeline smoke (gate v) + regenerate `kext_astrodust_MW.dat`.
      — **done**: `kext_astrodust_MW.dat`, `kext_astrodust_MW_euv.dat` (both
      1762 points, 1.0e-4–3.981e4 um) and `kext_dl07_MW_euv.dat` (1823 points,
      6.2054e-5–3.981e4 um) are regenerated.  Those file names are upstream's;
      MoCHII now takes the same curve from `/kext` of the model's product and
      does not carry the text form.
- S6. MoCHII host adoption: band optics from `dust_extinction`
      (option-gated), d_dusty re-gate, then decide the default.
      — **done** 2026-08-02, resolved to the grain model (see the adoption note
      in the preamble).  The option gating is gone with it: `par%dust_model`
      names the grains and `build_dust` resolves everything else from
      `par%sed_data_dir` (2026-08-07).

Of the three open decisions listed for the author, two are settled: the
extension floor went to 1 A rather than either candidate, and MoCHII's default
band optics are the grain model named by `par%dust_model`.  ADT as a check on
the x > 50 EUV rows is still open, and is now the natural way to bound the
error S2 left in place.
