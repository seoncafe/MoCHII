# The 'chianti' excitation model: data generation (October 11, 2026)

Phase 2 (data part) of
[ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md](../../ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md).
Product: `data/atomic/excitation/chianti/` (model name `'chianti'`), 70 atoms,
70 cooling fits, 3 atomic-hydrogen impact files and `manifest.json`
(SHA-256 of all 143 files). The runtime does not yet accept the model name:
`src/nlevel_mod.f90::metal_excitation_data_file` still allows only
`mocassin_v2` and `extended`. No `src/` file was changed here.

## 1. Coverage of the 77 stages

| Source | Stages |
|---|---|
| CHIANTI 11.0.2, every bound collision-connected level, every radiative and collision row (67) | C I-VI, N I-VII, O I-VII, Ne I-VII, S I-VII, Ar III-VII, Mg I-VII, Fe II-VII, Si I-VII, Cl II-IV, Cl VI, Ca II, Ca V-VII |
| Cloudy c25.00 Stout, CHIANTI has no valence collision data (3) | Cl VII (32 bound outer-shell levels of the 161 in `cl_7.nrg`; Upsilon nodes 9800-9.8e7 K, held below 9800 K), Ar II and Ca IV (ground fine-structure pair, Pelan & Berrington 1995 Table 3; nodes from 1580 K and 1e4 K, held below the first node) |
| none: no line cooling (7) | C VII (bare nucleus); Ar I, Cl I, Ca I, Fe I (no electron-impact data in CHIANTI 11.0.2 or Stout); Cl V, Ca III (CHIANTI has only `.easplups`, inner-shell excitation-autoionization) |

Totals: 6629 levels, 129842 radiative rows, 574482 collision rows,
305 MB (atoms; the fits are 280 kB). Largest: Si VII 460 levels / 105570
collision rows / 59.8 MB, Mg VII 450 / 101025 / 57.1 MB, Ne V 304 / 46056 /
25.8 MB, Fe III 322 / 51681 / 24.2 MB. The table of all stages with the check
results is [coverage_table.md](coverage_table.md); all numbers are in
[coverage.json](coverage.json) and [checks/](checks/).

"Levels" counts the levels written. Levels at or above the ionization limit
(autoionizing; no autoionization rates in the n-level solve) and levels that
no collision or radiative decay connects to the ground are not written; each
atom header states written versus source counts (for example C IV 24 of 923,
N III 86 of 204, Ar VI 155 of 202). Fe II is the full 142-level atom with all
1268 `.wgfa` rows and 4339 `.scups` transitions; Fe III the full 322 levels,
10480 rows and 51681 transitions; Ca VII 27 levels, 86 rows, 159 transitions.

**Fe VII changed source.** The generator took Fe VII from Stout (Witthoeft &
Badnell 2008). CHIANTI 11.0.2 has `fe_7` (Tayal & Zatsarinny 2014 collisions
and A values, Li et al. 2018 / Tayal & Zatsarinny theoretical energies), so
under the plan's rule Fe VII is now CHIANTI; Stout is used only for the three
ions CHIANTI lacks (`chianti_cooling.STOUT_EXCITATION_IONS`). Low-density
cooling Stout / CHIANTI: 1.14 (5e3 K), 0.90 (1e4 K), 0.79 (2e4 K),
0.89 (5e4 K), 0.96 (1e5 K).

**Existing choices for single ions kept** (all earlier checked code paths of the
generator, now applied to full atoms): O II/O III E1 A values at the selected
photon frequencies and the Storey & Zeippen (2000) [O III] optical A values;
Ne I two added Stout A values; Si I NIST MCHF strengths; Mg III level 167
theoretical energy; N IV Stout energies; O V and Cl VI observed/ICFT energy
splices; O IV, Ne V, S VII, Ar V, Mg IV valence-complex energies; S VI
centroid energies; BT lower-temperature holds for C IV, Mg III, Si V, O IV,
O V, Ne I, Si I (written as `BT_TMIN_K`).

## 2. Generator changes (tools/fitting)

1. **Radiative rows whose selected energies reverse the source order**
   (`electric_dipole_probabilities.reversed_electric_dipole_probabilities`,
   applied in `chianti_cooling.bound_excitation_atom`). CHIANTI applies a
   `.wgfa` A as a decay from its second index; with the selected energies
   (observed where measured, otherwise theoretical, CHIANTI's own rule) 51
   rows decay upward: N III 9, Ne III 4, Ne VI 1, Ne VII 2, Si VII 14,
   Fe VII 21. All are E1 (opposite parity). They are written from the higher
   state with A scaled by (nu_selected/nu_source)^3 and the emitting
   statistical weight at fixed line strength; nu_source is the frequency at
   which the source A and gf describe the same strength. The writer used to
   stop on them. A reversed non-E1 row is still refused (none occur).
   For the 21 Fe VII rows (levels 117-142, about 550000 cm^-1) nu_source
   matches neither the `.wgfa` wavelength nor the CHIANTI theoretical gap
   (A from Tayal & Zatsarinny, energies mostly Li et al.): their state
   identification is unverified and listed in the atom header.
2. **Negative scaled collision strengths** (`chianti_cooling.read_excitation_atom`).
   CHIANTI 11.0.2 has negative last nodes (x = 1, the infinite-energy limit)
   in type-2 rows: Mg I (2,12) -0.0095 and (13,15) -0.0007, Fe II (1,9)
   -0.0013, Fe VII 25 rows between -1e-7 and -1.1e-4. They are set to zero
   and the source values are listed in the atom header; a negative node
   anywhere else is refused. The Fortran loader refuses negative nodes, so
   these atoms could not be written before.
3. **Fe VII routed to CHIANTI**; the Fe VII Stout special cases in
   `fit_cooling.py`, `fit_nlevel.py`, `validate_stage_cooling.py`,
   `assess_excitation_temperature_range.py` and
   `validate_bound_level_convergence.py` now follow the routing.
4. **Atom headers** state the source, written versus source level counts, the
   reoriented rows and the zeroed nodes; Cl VII states the 129 core-excited
   Stout levels not read.
5. **`fit_cooling.py`**: in the logarithmic branch (ions whose cooling
   underflows below 1e5 K) the removed excitation factor uses the physical
   ground threshold (`ground_excitation_threshold`) when the padded ladder lies
   below it. O VII failed without it (ladder 3.26e6 K against 6.51e6 K,
   infinite least-squares matrix). C V and N VI were refitted with the change;
   their maximum errors are unchanged (0.86%, 0.79%).
6. **`make_hcoll_fine_structure.py`**: `--output-dir`, and each addressed level
   is checked in the target atom (statistical weight and energy within
   0.5 cm^-1 of NIST) before writing. C I 1-3 = 3P0,1,2 (0, 16.417,
   43.414 cm^-1), C II 1-2 = 2P1/2, 3/2 (0, 63.397), O I 1-3 = 3P2,1,0
   (0, 158.265, 226.977): the full atoms keep the ground term at indices
   1-3, and the rows are identical to the `extended` files.

## 3. Checks (every atom)

Programs: `tools/fitting/validate_chianti_excitation_model.py` (per atom,
results in [checks/](checks/)), the Fortran driver
`tools/fitting/excitation_atom_cooling.f90`, which loads each file with
`nlevel_load` from copies of `src/define.f90` (SHA-256 704cefd6...) and
`src/nlevel_mod.f90` (55e39ab0..., unchanged in `src/` when the checks
finished) compiled with ifx -O2 outside the tree, and
`tools/fitting/chianti_excitation_model_coverage.py` (table, manifest).

| Check | Result (worst atom) |
|---|---|
| Fortran `nlevel_load` accepts the file | 70 of 70 (hcoll files read with C I, C II, O I) |
| Radiative photon energy E_u - E_l > 0 and A > 0 | 0 failures in 129842 rows |
| Detailed balance: Fortran collision matrix with radiation off gives Boltzmann populations (3e3, 1e4, 1e5 K) | relative 7.7e-12 (Ne I, 1e5 K), absolute 5.5e-13 |
| Upsilon from file nodes against the source evaluation, 1e3-1e5 K | 5.4e-16 (Cl VII) |
| Low-density cooling: Fortran at n_e = 1e-6 against `cooling_low_density_limit` on the source, 81 T | <= 4.1e-7 for 67 atoms; Fe IV 3.9e-6; S V 2.6e-3 and Ar VII 4.9e-3 where Lambda > 1e-3 of its maximum (2.7% and 12% at 1e3 K) |
| Cooling fit against the source, 801 T over 1e3-1e5 K | <= 0.98% for 69 atoms; O IV 1.13% |
| Fit terms | 2-21 (469 in all); none negative on 801 T or on 1-1e9 K sampling at fit time |

S V and Ar VII (Mg-like) carry A(3s3p 3P0 -> 3s2 1S0) = 1e-30 s^-1, a
CHIANTI placeholder for a J = 0 -> 0 decay. The strict n_e -> 0 limit
(the fits) lets that level radiate every excitation; any finite density
(1e-6 cm^-3 already) de-excites it collisionally first. The finite-density
solve is the physical one; the fits of these two ions overstate the
low-density cooling by 0.4% (S V) and 1.0% (Ar VII) at 1e4 K.

The cooling fits exceed 16 terms for Ca IV (20), Mg I (20), Ne I (20),
S I (21), O V (19) and Mg III (18). `cooling_mod.f90::cooling_load`
allocates by the file's count; `tests/cooling_ne/test_cool_table.f90`
declares `fitA(16)` and cannot read them.

Not measured: the plan's step-1 costs (table build time, memory on each rank and
line-output time in a MoCHII run). The standalone driver needed up to 23 s
and 109 MB per atom for loading plus 81 temperatures x 4 densities
(Mg VII, Si VII).

Python gates on the changed generators (`tests/fitting`, the 16 modules that
import them, 52 tests): 49 pass. The fixture of
`test_excitation_source_storage.py` gave a radiative row upper-first, which
the writer refuses; it now gives both rows lower-upper (summing is still
tested). Three failures do not involve the changed code and were not
changed: `test_niv_radiative_source` (`niv_level_energies` has no
`selected_energies`), `test_selective_generation` (`make_element_data`, N VII
charge exchange), and `test_icft_excitation_energies` (O V at 2000 K, finite
n_e = 1e-8 against the strict limit, 4.6% on a cooling of 6e-44; O V's
source cooling is identical with and without the reorientation of item 2.1).

## 4. Fe II levels 14-35 (audit E3)

The CHIANTI atom has 141 radiative rows with an upper level in 14-35
(13474-31368 cm^-1); every one of these levels decays (smallest total A
0.023 s^-1). The `mocassin_v2`/`extended` file has none. Low-density
cooling of the `mocassin_v2` fit relative to the new one: 0.99 (5e3 K),
0.88 (1e4 K), 0.91 (2e4 K).

## 5. Audit E5 ions

N I, Ca II, C I, S III, Si III are CHIANTI 11.0.2 (the atom headers carry the
SHA-256 of the CHIANTI `.scups` read). Low-density cooling of the
`mocassin_v2` fits relative to the new fits at 5e3 / 1e4 / 2e4 K: N I 0.018 /
0.093 / 0.32, Ca II 0.37 / 0.24 / 0.12, Si III 0.87 / 0.87 / 0.86; C I and
S III have no `mocassin_v2` fit. The `extended` fits agree with the new ones
to 1.5e-3 (they were fitted to the same full CHIANTI atoms).

## 6. Ca VII

Script `tools/fitting/compare_cavii_radiative_probabilities.py`, output
[cavii/cavii_radiative_probabilities.json](cavii/cavii_radiative_probabilities.json),
against the 55 A values of Landi & Bhatia (2003, ApJ 589, 1075) Tables 3-5
transcribed in `../validation_2026-10-11/cavii_published_review_checked/`.

- Energies: observed where measured, otherwise the 25-configuration
  theoretical values (CHIANTI's rule and the writer's). All 86 rows have a
  positive photon energy. For the 69 rows with gf > 0 the frequency at which
  A and gf agree is 0.99946-0.99995 of the listed wavelength's, so A and gf
  are on the listed (observed where both levels are observed) energies, as
  the article's Section 2.1 states.
- Collision matrix: Boltzmann to 2.4e-15; Fortran low-density cooling against
  the source 2.0e-13; fit 0.10%.
- **Configuration mixing is the same as in the v2 file.** Of the 54 printed
  pairs present in CHIANTI (14-16 is not), the allowed rows of Tables 4
  (17) and 5 (15) match the 25-configuration column, but 21 of the 22
  forbidden rows of Table 3 match the 5-configuration column (2-3 matches
  both within 0.1%). The `.wgfa` header says all A values are
  25-configuration. Ratios CHIANTI / 25-configuration reach 0.25 (1-4),
  2.34 (6-16), 3.34 (13-16), 0.61 (2-5), 1.16 (2-4, 3-4).
- **Consequence, measured** with the installed atom (Python n-level solve):
  replacing the 22 Table-3 A values by the 25-configuration column changes
  the Ca VII cooling by <= 3e-4 at 1e4-2e4 K and n_e = 1e2-1e6 cm^-3; the
  four strongest lines ([Ca VII] 6.15 and 4.09 um, 5620 and 4941 A, 99% of
  the cooling) by <= 0.16%; [Ca VII] 2112 A (2-5, 0.5% of the cooling at
  2e4 K) by +3.5 to +4.3%.
- Verdict: the mixing is an inconsistency of the radiative data set, not a
  physical defect of the atom (energy conservation, detailed balance and the
  A-gf-energy relation hold). The atom is installed as CHIANTI distributes
  it. Replacing the forbidden rows by the article's preferred
  25-configuration values changes CHIANTI data and is left to the author.
- The collision strengths are the distorted-wave 5-configuration calculation
  (scattering energies 8-60 Ryd, five BT nodes); the article (Section 2.4)
  warns that its low-temperature values are an extrapolation. Not quantified.

## 7. Not verified

- The atoms were loaded by a copy of `nlevel_mod.f90` (hash above); later
  edits of the loader are not covered.
- The runtime's use of `DENSITY_RATIO 1` with these fits was not run; no
  MoCHII run (MPI) was made.
- Atomic accuracy of the CHIANTI data themselves (beyond Ca VII rows and the
  checks above) is not verified; the 21 Fe VII rows of item 2.1 have
  unverified identifications.
- The O IV fit (1.13%) exceeds the 1% the fit generator aims at, which it
  tests on 201 temperatures only.
- `data/atomic/README.md`, `CLAUDE.md` and the physics documents do not yet
  describe this model.
