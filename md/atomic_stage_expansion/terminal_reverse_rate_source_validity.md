# Source validity of terminal reverse recombination rates

## Scope and conclusion

The [terminal-stage calculation](terminal_stage_next_implementation.md) samples 1,000–100,000 K with a hypothetical additional ion for each of the eleven original terminal stages. Its reverse coefficients are real source records: ten ground-state `M=1` Badnell-group RR and DR fits and one CHIANTI Fe VIII pair. The published *fit-to-calculation* accuracy intervals do not cover the full sampled range for O VI and Ne VI DR. Both Li-like O VI and B-like Ne VI papers set the DR-fit lower temperature at `5,000 q^2 K`, or 125,000 K for residual charge `q=5`. The source-native current tabulation includes lower-temperature points, but those points cannot extend the published article's stated validity. Fe VIII has a fitted curve at 1,000 K, but that point is below the experiment-derived plasma range. The fit accuracy percentages below are not atomic-theory or experiment uncertainties. This audit does not establish a physically valid thermal cutoff or justify suppressing any outgoing terminal ionization.

The diagnostic reads `source_terminal()` in [validate_terminal_ionization.py](../../tools/fitting/validate_terminal_ionization.py), which calls the exact parsers and formulas in [make_element_data.py](../../tools/fitting/make_element_data.py). The CI row comes from the local Voronov `cfit.f`; no newly inferred CI correction is applied. Its CI source key is `(Z, electrons before ionization)`. The reverse RR/DR source key is `(Z, electrons in the recombining ion, M=1)`. The table gives those keys and the published DR-fit lower temperature. All ten Badnell-group RR fits use [Badnell (2006)](https://doi.org/10.1086/508465), except Ca VII, whose Si-like RR comes from [Kaur et al. (2018)](https://doi.org/10.1051/0004-6361/201731243). The Fe VIII coefficients are CHIANTI `fe_8.rrparams` and `fe_8.drparams` records, checked against [Schmidt et al. (2008)](https://doi.org/10.1051/0004-6361:200810834) Tables 3–4 in [the preceding assessment](terminal_ionization_assessment.md).

| Outgoing CI | CI `(Z,N)` | Reverse `(Z,N,M)` | Reverse DR final article | Published DR-fit lower T (K) | DR share at 10,000 / 100,000 K |
|---|---:|---:|---|---:|---:|
| C V → C VI | (6,2) | (6,1,1) | [Badnell 2006, H-like](https://doi.org/10.1051/0004-6361:20053269) | 250, subject to negligible-rate qualification | 3.96e-134 / 3.71e-13 |
| N VI → N VII | (7,2) | (7,1,1) | [Badnell 2006, H-like](https://doi.org/10.1051/0004-6361:20053269) | 360, subject to negligible-rate qualification | 5.25e-181 / 9.49e-18 |
| O V → O VI | (8,4) | (8,3,1) | [Colgan et al. 2004, Li-like](https://doi.org/10.1051/0004-6361:20034197) | **125,000** | 0.331 / 0.945 |
| Ne V → Ne VI | (10,6) | (10,5,1) | [Altun et al. 2004, B-like](https://doi.org/10.1051/0004-6361:20040964) | **125,000** | 0.799 / 0.906 |
| S VII → S VIII | (16,10) | (16,9,1) | [Zatsarinny et al. 2006, F-like](https://doi.org/10.1051/0004-6361:20053737) | 569 | 0.215 / 0.236 |
| Ar VII → Ar VIII | (18,12) | (18,11,1) | [Altun et al. 2006, Na-like](https://doi.org/10.1051/0004-6361:20053278) | 490 | 0.500 / 0.968 |
| Mg IV → Mg V | (12,9) | (12,8,1) | [Zatsarinny et al. 2003, O-like](https://doi.org/10.1051/0004-6361:20031462) | 11.6 | 0.542 / 0.785 |
| Fe VII → Fe VIII | (26,20) | (26,19,1) | [Schmidt et al. 2008](https://doi.org/10.1051/0004-6361:200810834) | 2,321 experimentally derived plasma range; analytic fit extends lower | 0.891 / 0.935 |
| Si V → Si VI | (14,10) | (14,9,1) | [Zatsarinny et al. 2006, F-like](https://doi.org/10.1051/0004-6361:20053737) | 290 | 0.189 / 0.342 |
| Cl VII → Cl VIII | (17,11) | (17,10,1) | [Zatsarinny et al. 2004, Ne-like](https://doi.org/10.1051/0004-6361:20040463) | 11.6 | 1.03e-49 / 1.11e-5 |
| Ca VI → Ca VII | (20,15) | (20,14,1) | [Kaur et al. 2018, Si-like](https://doi.org/10.1051/0004-6361/201731243) | 360 | 0.972 / 0.962 |

## What the published ranges mean

The final journal versions were read as page-formatted PDFs. Badnell (2006) gives RR fit differences from computed data below 1% for multiply charged ions over `(10^1–10^7)q^2 K`, with a correct asymptotic form outside. Kaur et al. give Ca VII Si-like RR and DR fit differences below 5% over `(10–10^7)q^2 K`. Colgan et al. and Altun et al. give Li-like and B-like DR fit differences below 3% and 5%, respectively, only over `(5,000–10^7)q^2 K`. The O-like paper gives 3% over 0.001–100,000 eV. The F-like paper gives 3% over `0.001 < T/q^2 < 1,000 eV`. The Ne-like paper gives 1% over 0.001–100,000 eV. The Na-like paper gives about 5% over `(10–10^7)q^2 K`. The H-like paper gives 5% over `(10–10^7)q^2 K` or down to the temperature where DR becomes negligible. Each paper measures an analytic-fit difference against its own calculation; none establishes a uniform uncertainty in the underlying atomic rates for these specific eleven ions. For example, the F-like paper estimates a roughly ±20% theoretical uncertainty from comparison to Fe XVIII laboratory data, not from a direct measurement of S VIII or Si VI.

The current source database is not necessarily identical to every printed table in an older paper. Its [RR](https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/) and [DR](https://amdpp.phys.strath.ac.uk/tamoc/DATA/DR/) pages provide ground-level `M=1` formulas and `cfout_K.gz` original-versus-fit comparisons. The source-native grid is nominally `(10–10^7)q^2 K`. On its tabulated nodes from 1,000 to 100,000 K, the largest *RR* residual among these ten records is 3.45% for Ca VII. The O VI DR residual reaches 16.6% at 2,500 K, below the Li-like paper's 125,000 K lower limit. Cl VIII DR has an enormous relative residual where its original rate is effectively zero: at 4,900 K the original is about 1.13e-189 and the fit about 1.45e-111 cm3 s-1, while RR is about 1e-11 cm3 s-1. This does not indicate an enormous total-recombination error. The sampled, rounded original-plus-fit RR+DR comparisons and each component's fit residual are stored in `results.json`. They cannot be treated as theory uncertainty or as a continuous error bound. In particular, the archived source coefficient set and older published fit-error statements should be distinguished.

## Reproducible diagnostic and limits

The [isolated diagnostic](validation_2026-10-10/terminal_reverse_provenance/rate_provenance.py) evaluates the same source coefficients and formulas as the staged terminal diagnostic, directly separates RR and DR, and writes [801-temperature profiles](validation_2026-10-10/terminal_reverse_provenance/full_profiles.npz), [numerical records and file hashes](validation_2026-10-10/terminal_reverse_provenance/results.json), and [the visually inspected full-profile PDF](validation_2026-10-10/terminal_reverse_provenance/reverse_rate_profiles.pdf). It compares current local coefficients with author-hosted `cfout_K.gz` files, whose download addresses are `https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/cfout_K.gz` and `https://amdpp.phys.strath.ac.uk/tamoc/DATA/DR/cfout_K.gz`. Those external gzip files were retained in `/tmp` for this calculation, and their hashes are in `results.json`. Run it from the repository root using `MPLCONFIGDIR=/tmp/mochii_terminal_reverse_mpl python3 md/atomic_stage_expansion/validation_2026-10-10/terminal_reverse_provenance/rate_provenance.py`. Its log is [run.log](validation_2026-10-10/terminal_reverse_provenance/run.log).

| Input condition | This diagnostic |
|---|---|
| Temperature | 801 logarithmic points from 1,000 to 100,000 K; displayed checkpoints at 1,000, 10,000, and 100,000 K |
| Radiation and H charge transfer | None; rate coefficients only, without a spatial or radiation calculation |
| CI | Local Voronov `cfit.f` outgoing row; no terminal `dere_hybrid` calibration |
| RR and DR | Ten local Badnell-group `M=1` rows; Fe VIII local CHIANTI rows |
| Electron density | Not specified: these are zero-density rate coefficients, not finite-density suppressed DR |

Thus the source lookup verifies that reverse coefficients exist and are invoked by the diagnostic. It does not verify the physical error in CI, finite-density DR behavior, terminal H charge transfer, or Auger channels. A thermal-domain judgment should use the complete ratio and energy-loss arrays in the existing terminal diagnostic, with O VI and Ne VI reverse-fit extrapolation explicitly identified wherever those reverse reactions affect a proposed boundary.

## Main-session review and reference archive

The main session inspected the complete eleven-ion rate profiles and directly
checked the Li-like and B-like published fit-range statements and the Si-like
RR fit statement against the final journal PDFs. Nine verified final PDFs
used in this assessment have been preserved in `~/MoCHII/references/`;
[the archive record](validation_2026-10-10/terminal_reverse_provenance/reference_archive.json)
gives filenames and hashes. This review changes neither the production
terminal roster nor any reverse coefficient.
