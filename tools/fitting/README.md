# tools/fitting

Python fitting scripts: CHIANTI v11.0.2 (`~/RT_Codes/CHIANTI/dbase`) -> rate/cooling
fits. Seeded from EXHALE `cooling_data` (`chianti_cooling.py`). The auditable
source of every coefficient in `data/atomic/` — Fortran never parses CHIANTI
at runtime.

## Scripts

| Script | Product | Role |
|---|---|---|
| `chianti_cooling.py` | (library) | direct CHIANTI ASCII reader (.elvlc/.wgfa/.scups), Burgess-Tully descaling, cooling builders; copied from EXHALE |
| `fit_cooling.py` | `data/atomic/cooling_<ion>.txt` | Tier-1 cooling fits Lambda(T) = T^-1/2 sum A_i exp(-T_i/T), 1e3-1e5 K, in the n_e -> 0 limit of statistical equilibrium; for the G2 thermal loop |
| `fit_nlevel.py` | `data/atomic/nlevel_<ion>.txt` | Tier-2 n-level data: level energies/weights, A-values, Upsilon(T) as Chebyshev fits in Burgess-Tully scaled space; for output-time diagnostics |
| `make_cooling_table.py` | `data/atomic/cooling_fit_parameters.txt` | collects the Tier-1 coefficient files into one listing |
| `fit_alpha1_milne.py` | `src/recomb_mod.f90` `alpha1_*` coefficients | ground-level alpha_1 from the Milne relation applied to the transport photoionization cross sections (the same integral the diffuse ground continua are sampled from), fitted to the Mao & Kaastra (2016) form over 1e3-1e5 K; agrees with Cloudy c25.00's level-resolved values to <= 0.4%. Gate `tests/milne/check_alpha1.f90` |
| `make_hydrogenic_recombination.py` | `src/hydrogenic_recombination_data.inc` | H II -> H I and He III -> He II from one zero-density nl-resolved radiative cascade on Storey & Hummer's (1991) CPC routines ACOFZ1/PIXZ1 (compiled from `~/RT_Codes/Storey_Hummer/CPC_code/abzs_v1_0`, unmodified): alpha_B, the 2s share in case A and case B, and the He II case-A Lyman-line shares, log T = 1.5-7.0; validated against the SH95 alpha-tot, O&F Table 4.11 and the SH95 H-beta emissivities (printed). Replaces `make_heii_cascade.py` (2026-09-26). Gate `tests/hydrogenic_recombination/` |
| `make_element_data.py` | `data/atomic/element_<el>.txt` | photoionization (VFKY96 outer shell, VY95 subshells as `phfit2`), Badnell RR+DR (CHIANTI forms outside Badnell's coverage), collisional ionization, charge exchange (with the product-term form for sulfur) and its heat, for the 11 registry elements |
| `make_hcoll_fine_structure.py` | `data/atomic/hcoll_<el>_<stage>.txt` | atomic-hydrogen impact rates on the C I, C II and O I ground terms, added to the n-level statistical equilibrium |
| `make_heii_h_charge_exchange.py` | `src/heii_h_charge_exchange_data.inc` | non-radiative He+ + H -> He + H+ from the reverse calculation of Loreau, Ryabchenko & Vaeck (2014, J. Phys. B 47, 135204; eq. 3, Table 2 process 1) by microscopic reversibility (weights 1/4), with sigma_r ~ (x - x_th)^(1/2) below 15 eV/u where the quantal data stop, Maxwell-averaged, log T = 3.0-6.0 at 0.02 dex; prints the near-threshold bracket, the frame sensitivity and the reciprocity check. Replaces the KF96 He+ row (2026-09-26). Gate `tests/charge_exchange/check_helium_charge_exchange.f90` |
| `make_h2s_collisional_mixing.py` | `data/atomic/h_2s_2p_mixing.txt` | collisional H(2s) -> H(2p) rate coefficients by collider (e, p, He+, He2+) and A_2s, for the H I two-photon quenching |
| `make_hei_lines.py` | `data/atomic/hei_porter_caseB.txt` | He I case-B collisional-radiative line emissivities of Porter et al. (2012, 2013) |
| `make_hei_metastable.py` | `data/atomic/hei_metastable.txt` | He I 2^3S coefficients (every term of its balance except the recombination source) |
| `make_sh95_lines.py` | `data/atomic/sh95_{hi,heii}_case{A,B}.txt` | Storey & Hummer (1995) H I and He II recombination-line emissivities from the CDS files |
| `reference_cno_coefficients_exhale.txt` | (reference) | EXHALE's ground-term-Boltzmann coronal C/N/O cooling fits, kept for comparison; superseded by the n_e -> 0 fits of `fit_cooling.py` |
| `verify_nlevel_pyneb.py` | (report) | end-to-end check: reads the fitted files, solves statistical equilibrium as the Fortran will, compares diagnostic ratios against PyNeb |

## Current ion set (counted 2026-09-27)

28 Tier-1 cooling ions and 28 Tier-2 n-level models (plus the two
`nlevel_fe_{2,3}_full.txt` variants) across hydrogen and the 11-element
metal registry (C, N, O, Ne, S, Ar, Mg, Fe, Si, Cl, Ca) — every cooling ion has
a model, so the runtime density suppression covers all of them.
Adding an ion = add one entry to the `IONS` dict of each fit script and rerun.

## Tier-1 population model

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

## Fit quality (recorded per file in the provenance headers)

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
