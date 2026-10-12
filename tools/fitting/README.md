# tools/fitting

Current selection (October 11, 2026): `metal_excitation_model` defaults to
`chianti` (full CHIANTI 11.0.2 atoms, 70 atoms; Stout for Cl VII, Ar II,
Ca IV); `mocassin_v2` (the MOCASSIN v2 atoms plus the CHIANTI Fe II forbidden
rows) is kept for parity runs. `src/nlevel_mod.f90` resolves both level and
cooling data in `data/atomic/excitation/<model>/`. The `extended` selection
and the root metal atom files were removed on 2026-10-11. The rate roster is
still 43 metal stages and 32 adjacent transitions until the 77-stage files
are installed (md/ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md, Phase 3).

`data/atomic/excitation/chianti/` (generated October 11, 2026) holds full CHIANTI 11.0.2 atoms for 67 stages, Stout
atoms for Cl VII, Ar II and Ca IV, their cooling fits and the C I, C II, O I
hydrogen-impact files: 70 atoms, 6629 levels, 469 fit terms. See
[its report](../../md/atomic_stage_expansion/chianti_model_2026-10-11/README.md).

The generator descriptions and fit comparisons below record the earlier
CHIANTI fitting workflow. They do not describe the full set now selected at
runtime.

Python fitting scripts: CHIANTI v11.0.2 (`~/RT_Codes/CHIANTI/dbase`) -> rate/cooling
fits. Seeded from EXHALE `cooling_data` (`chianti_cooling.py`). The auditable
source of every coefficient in `data/atomic/` — Fortran never parses CHIANTI
at runtime.

## Scripts

| Script | Product | Role |
|---|---|---|
| `chianti_cooling.py` | (library) | direct CHIANTI ASCII reader (.elvlc/.wgfa/.scups), Burgess-Tully descaling, cooling builders; copied from EXHALE |
| `write_excitation_atom.py` | `nlevel_<el>_<stage>.txt` (`--output-dir`) | NLEVEL_VERSION 2 atom with the source Burgess-Tully nodes or Stout temperature tables; every bound collision-connected level unless `--levels` truncates |
| `electric_dipole_probabilities.py` | (library) | E1 A values at selected photon frequencies at fixed line strength; reorients rows whose selected energies reverse the source order |
| `validate_chianti_excitation_model.py` | `md/atomic_stage_expansion/chianti_model_2026-10-11/checks/` | checks of each atom: photon energies, collision-only Boltzmann balance of the Fortran loader, Upsilon serialization, low-density cooling against the source, fit error |
| `excitation_atom_cooling.f90` | (driver, built outside the tree) | loads an atom with `nlevel_load` and writes its cooling at four densities and the radiation-off populations, for `validate_chianti_excitation_model.py` |
| `excitation_level_populations.f90` | (driver, built outside the tree) | loads an atom with `nlevel_load` and prints its cooling and level populations at one (T, n_e) |
| `feii_spontaneous_probabilities.py` | `excitation/mocassin_v2/nlevel_fe_2.txt` | adds the 141 CHIANTI 11.0.2 fe_2 forbidden A values from levels 14-35 that the MOCASSIN v2 Fe II atom lacks; maps the 142 levels by index, g and energy (<= 0.005 cm^-1) and keeps every v2 row |
| `feii_low_density_cooling.py` | `excitation/mocassin_v2/cooling_fe_2.txt`, `manifest.json` | n_e -> 0 cooling fit of that atom read from the atom file (the `mocassin_cooling_fits` linear program with nonnegativity on 1-1e9 K); updates the two manifest SHA-256 entries |
| `chianti_excitation_model_coverage.py` | `manifest.json`, coverage table | 77-stage source coverage, SHA-256 manifest, cooling compared with the earlier models |
| `compare_cavii_radiative_probabilities.py` | `.../chianti_model_2026-10-11/cavii/` | CHIANTI Ca VII A values against Landi & Bhatia (2003) 5- and 25-configuration columns |
| `fit_cooling.py` | `data/atomic/cooling_<ion>.txt` (`--ion-spec el:stage --output-dir`) | Tier-1 cooling fits Lambda(T) = T^-1/2 sum A_i exp(-T_i/T), 1e3-1e5 K, in the n_e -> 0 limit of statistical equilibrium; for the G2 thermal loop |
| `fit_nlevel.py` | `data/atomic/nlevel_<ion>.txt` | Tier-2 n-level data: level energies/weights, A-values, Upsilon(T) as Chebyshev fits in Burgess-Tully scaled space; for output-time diagnostics |
| `make_cooling_table.py` | `data/atomic/cooling_fit_parameters.txt` | collects root coefficient files into a historical listing; model selection uses the directories above |
| `fit_alpha1_milne.py` | `src/recomb_mod.f90` `alpha1_*` coefficients | ground-level alpha_1 from the Milne relation applied to the transport photoionization cross sections (the same integral the diffuse ground continua are sampled from), fitted to the Mao & Kaastra (2016) form over 1e3-1e5 K; agrees with Cloudy c25.00's level-resolved values to <= 0.4%. Gate `tests/milne/check_alpha1.f90` |
| `make_hydrogenic_recombination.py` | `src/hydrogenic_recombination_data.inc` | H II -> H I and He III -> He II from one zero-density nl-resolved radiative cascade on Storey & Hummer's (1991) CPC routines ACOFZ1/PIXZ1 (compiled from `~/RT_Codes/Storey_Hummer/CPC_code/abzs_v1_0`, unmodified): alpha_B, the 2s share in case A and case B, and the He II case-A Lyman-line shares, log T = 1.5-7.0; validated against the SH95 alpha-tot, O&F Table 4.11 and the SH95 H-beta emissivities (printed). Replaces `make_heii_cascade.py` (2026-09-26). Gate `tests/hydrogenic_recombination/` |
| `ground_term_photoionization_table.py` | the `PHTABLE` record of `data/atomic/element_ne.txt` (through `make_element_data.py`) | Ne III ground-term photoionization cross section from Nahar's (2019) Breit-Pauli R-matrix calculation (`data/atomic/norad/ne3.px.gd.txt`): the 2J+1 average of 3P2,1,0 on the photoelectron-energy scale, every source node kept, cut where it meets the VFKY96 fit above the last target threshold (178.42 eV). Gate `tests/ne3_ground_term_photoionization/` |
| `make_element_data.py` | `data/atomic/element_<el>.txt` | photoionization (VFKY96 outer shell, VY95 subshells as `phfit2`; the Ne III `PHTABLE` below 178.42 eV), Badnell RR+DR (CHIANTI forms outside Badnell's coverage), collisional ionization, charge exchange (with the product-term form for sulfur) and its heat, for the 11 registry elements |
| `make_hcoll_fine_structure.py` | `data/atomic/hcoll_<el>_<stage>.txt` (`--output-dir`) | atomic-hydrogen impact rates on the C I, C II and O I ground terms, added to the n-level statistical equilibrium; checks the addressed level indices in the target atoms |
| `make_heii_h_charge_exchange.py` | `src/heii_h_charge_exchange_data.inc` | non-radiative He+ + H -> He + H+ from the reverse calculation of Loreau, Ryabchenko & Vaeck (2014, J. Phys. B 47, 135204; eq. 3, Table 2 process 1) by microscopic reversibility (weights 1/4), with sigma_r ~ (x - x_th)^(1/2) below 15 eV/u where the quantal data stop, Maxwell-averaged, log T = 3.0-6.0 at 0.02 dex; prints the near-threshold bracket, the frame sensitivity and the reciprocity check. Replaces the KF96 He+ row (2026-09-26). Gate `tests/charge_exchange/check_helium_charge_exchange.f90` |
| `make_h2s_collisional_mixing.py` | `data/atomic/h_2s_2p_mixing.txt` | collisional H(2s) -> H(2p) rate coefficients by collider (e, p, He+, He2+) and A_2s, for the H I two-photon quenching |
| `make_hei_lines.py` | `data/atomic/hei_porter_caseB.txt` | He I case-B collisional-radiative line emissivities of Porter et al. (2012, 2013) |
| `make_hei_metastable.py` | `data/atomic/hei_metastable.txt` | He I 2^3S coefficients (every term of its balance except the recombination source) |
| `make_sh95_lines.py` | `data/atomic/sh95_{hi,heii}_case{A,B}.txt` | Storey & Hummer (1995) H I and He II recombination-line emissivities from the CDS files |
| `reference_cno_coefficients_exhale.txt` | (reference) | EXHALE's ground-term-Boltzmann coronal C/N/O cooling fits, kept for comparison; superseded by the n_e -> 0 fits of `fit_cooling.py` |
| `verify_nlevel_pyneb.py` | (report) | end-to-end check: reads the fitted files, solves statistical equilibrium as the Fortran will, compares diagnostic ratios against PyNeb |

## Historical generator ion set (counted 2026-09-27)

28 Tier-1 cooling ions and 28 Tier-2 n-level models (plus the two
`nlevel_fe_{2,3}_full.txt` variants) across hydrogen and the 11-element
metal registry (C, N, O, Ne, S, Ar, Mg, Fe, Si, Cl, Ca) — every cooling ion has
a model, so the runtime density suppression covers all of them.
This is a dated generator inventory, not the current production set. New
source models require a physically checked atom and matching cooling data in
the selected model directory; editing the generator `IONS` dictionaries alone
does not activate a stage.

## Historical Tier-1 population model

The fitted quantity is the n_e -> 0 limit of the n-level statistical
equilibrium: at n_e << n_crit every ion sits in its lowest level, collisional
excitation happens out of that level only, and each excitation radiates its
full excitation energy through the cascade, so Lambda_ion is the excitation
power out of the ground level.  The T_i of the fit are therefore excitation
temperatures, and `ground_excitation_ladder` derives the starting ladder from
the level energies rather than from hand-tuned values.  Transition energies
are the observed `.elvlc` level differences in both the Boltzmann factor and
the emitted photon; the `.scups` dE is the Burgess-Tully scaling energy and is
used only to build the scaled temperature.

The limit carries no density suppression, so the low-n_crit lines
([C II] 158 um, [O III] 52/88 um, the [O II] and [S II] optical doublets)
are overestimated in dense gas; the density-dependent emissivities used for
the diagnostics come from the n-level solve.

## Historical fit quality (recorded per file in the provenance headers)

- Tier-1 cooling: max fit error <= 0.7% over 1e3-1e5 K for the H/C/N/O/Ne/S
  set, except C III (5.5% full range, 0.50% where the ion cools) and S II
  (2.2%); the heavily depleted Si, Cl, Ca ions reach 4-28%, concentrated at
  the ends of the fitted range.
- Tier-2 Upsilon: worst transition <= 1.7%; the residual outliers are weak
  intercombination transitions (O III 1D/1S - 5S2) and the coarse 10-point
  S II .scups tables.
- End-to-end vs raw CHIANTI (same solver, fitted file vs direct spline):
  [O II] 3726/3729 agrees to 0.05%, [S II] 6717/6731 to 0.5% — the fitted
  files are faithful to CHIANTI v11.
- Cross-check vs PyNeb 1.1.31 (independent atomic data, not CHIANTI):
  [O III] 5007/4363 within 1.2%, [N II] 6584/5755 within 0.3%,
  [O II] 3726/3729 within 6%, [S II] 6717/6731 within 10% (worst at
  n_e = 1e4 cm^-3, near the critical density, where the ratios become
  A-value-limited; PyNeb uses Kal09/TZ10 collision data and Z82-WFD96/RGJ19
  A-values there, so this is database spread, not fit error),
  [S IV] 10.51um/1406 within 9% (PyNeb carries DHKD82 collision strengths,
  CHIANTI v11 carries Tayal 2000; the [S IV] 10.51um emissivity itself is
  1.34x PyNeb's, and the Tier-1 [S IV] cooling matches the Cloudy c25.00
  channel of the HII40 benchmark to 3%).

## Requirements

numpy, scipy; PyNeb only for `verify_nlevel_pyneb.py`.
