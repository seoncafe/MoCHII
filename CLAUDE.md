# CLAUDE.md

Guidance to Claude Code when working in this repository, and the persistence
layer for the project: a fresh session must be able to resume from this file.
The full chronological development log lives in **`CLAUDE.md_2026.07.27`**
(2700+ lines) — consult it for the audit trail behind any item below.

## What MoCHII is

**MoCHII — MOnte Carlo for H II regions.** A Monte Carlo radiative-transfer
code for dusty photoionized nebulae on adaptive (octree AMR) grids:
H/He photoionization + thermal balance, metal lines added one ion at a time,
dust + PAH emission, recombination lines, and the nebular continuum — all from
one self-consistent radiation field. Author: Kwang-il Seon (KASI). Sibling of
the author's MoCafe (dust MC RT), LaRT (Lyman-alpha RT), and EXHALE (exoplanet
atmospheric escape) codes.

Built by **copying validated modules from MoCafe v2.00** (octree AMR engine,
event-driven traversal, peel-off, Lucy-style iteration, SED transport,
HDF5/FITS I/O) and adding a gas-physics layer with modern analytic atomic
data. It does not modify the MoCafe tree.

## Read first

- **`docs/PLAN.md`** — founding design: what MoCafe provides, what is new, the
  staged plan G0–G5 with validation gates, dust/PAH coupling, the metal-line
  registry ("adding an ion is a data operation, not a code operation"), and the
  CHIANTI fitting pipeline (Tier-1 cooling fits for the hot loop, Tier-2
  collision-strength fits + n-level solve for diagnostics).
- **`src/PORTING.md`** — module-by-module port map (what comes from which tree,
  what changes on the way in).
- **`docs/MoCHII_physics.pdf`** — atomic data with sources and fit functions,
  algorithms, all gate results. Companions: `MoCHII_UserGuide.pdf` (build/run,
  full `par%` reference, output formats, adding-an-ion runbook),
  `MoCHII_fitting.pdf`, `MoCHII_cooling_analysis.pdf`.

## Status

**All founding stages G0–G5 are built and gated.** Capabilities by area:

- **Ionization + thermal (G0–G2):** ionization equilibrium with opacity
  feedback; thermal balance by bisection on Te (ionization re-solved each
  trial); cooling = recombination + free-free (van Hoof 2014 Gaunt option) +
  collisional ionization + Tier-1 metal line cooling + H-impact [C II] 158 /
  [O I] 63.
- **Metal registry — 11 elements:** C, N, O, Ne, S, Ar, Mg, Fe, Si, Cl, Ca.
  Generic element cascade with charge exchange; Badnell 2023 RR+DR; VFKY96 /
  VY95 photoionization cross sections; adding an ion is a data operation.
- **Diffuse field (G3):** explicit ground-recombination packets (case A/B),
  He I excited channels, He II — replaces case-B on-the-spot when enabled.
  The three ground continua are sampled from the Milne free-bound spectrum
  `E^2 sigma_pi(E) exp[-(E-Eth)/kT]` with the transport cross sections
  (`par%diffuse_energy_model='milne'`, default; `'exponential'` and
  `'threshold'` retained for reproduction and sensitivity).
- **Diagnostics (Tier-2):** n-level statistical-equilibrium solver + line
  emissivities; SH95 H I / He II / Porter He I recombination lines; nebular
  continuum (Ercolano & Storey 2006 free-bound + free-free); leaf emissivity
  output (`<base>_emis`).
- **Dust + PAH:** EUV/FUV absorption and scattering, x_HII-dependent PAH
  survival, equilibrium T_dust + IR, self-contained SEDust stochastic emission
  with x_HII-weighted PAH share.
- **Grid / transport:** octree AMR with solution-driven re-refinement (G4);
  Cartesian (`car`) backend with DDA walk as default; shared-window recycling;
  RQMC (Owen-scrambled Sobol); 3D density-cube input; AMR grid builder +
  slice/cutout.
- **Sources:** composable N point sources + external field; file spectra;
  physical-unit spectra + ISRF presets.
- **Imaging:** peel-off (direct + dust-scattered, bin-resolved cubes,
  unattenuated direct image); Python readers + flux-conserving map maker
  (`tools/python/mochii_output.py`).
- **PDR switches (off by default):** metal electrons in n_e, metal photoheating,
  BT94 grain photoelectric heating; `species_ne` rate cache.
- **Convergence:** volume-integrated criterion (`par%conv_crit='vol'`,
  recommended for production) alongside the historical cell-max.
- **He I 2^3S / 10830** metastable diagnostic (Stages 1+2), incl. MoCHII → LaRT
  10830 coupling.

### Benchmark status (Lexington HII40 / HII20)

Passing after the 2026-07-25 Tier-1 cooling-fit correction. HII40
`<T[NpNe]>` = **7894 K** (inside the published 7720–8199); HII20 = **6376 K**
(0.4% below the published 6402–6749 floor); L(Hβ) ≈ −3% of published; ionization
fronts and R_out captured at both.

**Physics lesson from that fix (do not regress):** the Tier-1 cooling fits had
two defects, found by a fixed-state cooling-budget comparison against Cloudy
c25.00. **D1** — level populations were spread by Boltzmann factors; the correct
n_e ≪ n_crit limit puts the population in the lowest level (fits are now
regenerated from the n_e→0 limit of the n-level solve). **D2** — fits used
Burgess-Tully *scaling* energies; the observed `.elvlc` level differences are
the physical ones. Remaining known limitation: **Λ(T) has no density
suppression** — a Λ(T, n_e) Tier-1 (the `local_ne` cooling path) is the open
item; at nH=100 it biases Te slightly low.

## Next steps

- **Hard-spectrum (PN/AGN) extension** — the main open track. Plan in
  `docs/HARD_SPECTRUM_PLAN.md` (band headroom to ~500 eV, registry up to
  O/Ne VI, inner-shell/Auger, MOCASSIN PN + Cloudy gates). Its density-dependent
  Λ(T, n_e) prerequisite is done (`par%cooling_model='local_ne'`).
- **He I 10833 below the Porter n_e floor** — over-estimated at the ten-percent
  level for n_e < 10; refuse below the floor or feed from the explicit 2^3S
  population. Affects diffuse-medium runs only.
- **SH95 line-table interpolation** — bilinear-in-log is systematically low by
  0.06–0.13% on L(Hβ); replace with a cubic spline / Storey-Hummer form.
- **Cold-gas (T < 1e3 K) cooling boundary** — Tier-1 fits collapse below 1e3 K;
  extend the fits (or call the n-level solve) and add a [C I] H-impact channel.
- **He I 2^3S metastable — remaining gates** (H2b Cloudy anchor, H4 EXHALE
  triplet profile; plan `docs/HEI_10833_PLAN.md`).
- **Plane-parallel slab leftovers** (composing slab illumination with internal
  sources; Chandrasekhar H-function gate).
- Broad MPI collective/IO status wrapper; further ions as science needs them
  (P/K would reuse the VY95 path).

## Reference trees (add as working directories when needed)

| Tree | Path | What MoCHII takes |
|---|---|---|
| MoCafe v2.00 | `~/MoCafe/new/MoCafe_v2.00/` | AMR engine + transport + I/O (see `src/PORTING.md`); build/run conventions; `README_HOWTO.md` for the `par%` reference |
| EXHALE | `~/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/` | `cooling_data/` fitting pipeline; atomic-rate docs (Verner/Badnell/Voronov/charge-exchange); 1D reference profiles |
| MOCASSIN | `~/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/` | `gaunt.f90` (verbatim free-free); `fb_ff` nebular-continuum block + gamma tables; the 11-case regression suite as 3D validation reference |
| CHIANTI v11.0.2 | `~/RT_Codes/CHIANTI/dbase` | raw atomic data for the fitting pipeline (never parsed at runtime) |
| Other references | `~/RT_Codes/{Verner_Fortran_sub, Storey_Hummer, MAPPINGS, CMacIonize}` | cross-checks for rates, emissivities, code comparison |

## Conventions

- **MoCafe conventions, not MOCASSIN.** Namelist input into a single `par`
  (`params_type`); `real(kind=wp)`; MPI-only
  (`mpirun -np N ./MoCHII.x input.in`); HDF5 default with FITS fallback;
  procedure pointers select algorithms at startup. When a ported and a new
  module disagree in style, MoCafe style wins.
- **Documentation in English**, American spelling (all `.md`/`.tex`).
- **Validation by reference comparison, staged gates.** Every stage in
  `docs/PLAN.md` has a quantitative gate (analytic Stromgren, MOCASSIN 3D,
  EXHALE 1D, PyNeb ratios). AMR results must reproduce the equivalent
  uniform-grid run.
- **Atomic data are fitted offline, evaluated online.** The Python pipeline
  under `tools/fitting/` is the auditable source of every coefficient in
  `data/atomic/`; Fortran only evaluates formulas. Every coefficient file
  carries a provenance header.
- Git: the author manages commits (see `NFS_GIT_가이드.md` in the MoCafe tree).
  Do not commit or push unasked.

## Key design decisions already made (do not relitigate)

1. **Why not MOCASSIN + octree:** MOCASSIN's `pathSegment` resists structural
   traversal change (two attempts there failed and were reverted). MoCafe's AMR
   engine already is the event-driven design that route would have required.
2. **Trace-species metals, one ion at a time** — registry + generic element
   cascade + one generic n-level solver; metal cooling enters at G2 because
   H II Te (~8 kK) requires it.
3. **Free-bound is current** (Ercolano & Storey 2006); **free-free** ported
   verbatim (Hummer 1988) with the van Hoof et al. 2014 table swap as an
   option at the `getGauntFF` interface.
4. **PAH/dust coupling is a first-class goal** (EUV grain heating, dusty
   Stromgren competition, x_HII-dependent PAH survival).

## Layout

```
CLAUDE.md               this file (summary)
CLAUDE.md_2026.07.27    full chronological development log
README.md               short overview
docs/PLAN.md            founding design + staged plan (read first)
src/                    Fortran sources (PORTING.md = port map)
data/atomic/            fitted coefficient files (provenance headers required)
tools/fitting/          Python fitting pipeline (extends EXHALE cooling_data)
tools/python/           output readers + map maker (mochii_output.py)
examples/               run templates (MoCafe examples/ pattern)
tests/                  regression cases + baselines (MOCASSIN-suite style)
```
