# Measured readiness of the 21 added excitation ions

The installed eleven `ELEMENT` files declare 43 stages and 32 adjacent transitions. The staged `langevin64_exit_heat/products/atomic` roster declares 64 and 53. The 21 new ions are listed below. Running `python audit.py` imports the current `bound_excitation_atom`, reads both rosters, enumerates private artifacts, and writes the source paths and SHA-256 hashes in `result.json`. It does not fit, compile, register, or install anything. The counts are current reader results, not copied from planning documents. At this preserved pre-installation checkpoint, the 21 newly added ions had **zero** installed atoms or fits. The later [directly remeasured inventory](post_o_iv_installation.json) records one installed atom/fit, O IV, and 20 remaining additions. The production roster still has 43 stages.

| Added ion | Bound levels | Radiative rows | Collision pairs | Source nodes | Private atom/fit files found |
|---|---:|---:|---:|---:|---:|
| Ar V | 197 | 2,646 | 3,196 | 38,322 | 4 / 1 |
| Ar VI | 155 | 2,719 | 2,036 | 24,432 | 2 / 1 |
| Ar VII | 196 | 4,787 | 19,110 | 286,650 | 1 / 1 |
| C V | 49 | 206 | 273 | 4,095 | 2 / 2 |
| Ca VI | 5 | 10 | 10 | 90 | 1 / 1 |
| Cl VI | 171 | 3,760 | 14,535 | 218,025 | 1 / 1 |
| Cl VII | 32 | 496 | 496 | 6,448 | 1 / 1 |
| Fe V | 34 | 361 | 561 | 5,001 | 1 / 1 |
| Fe VI | 96 | 1,417 | 4,560 | 68,400 | 2 / 1 |
| Fe VII | 189 | 17,531 | 17,766 | 230,958 | 4 / 1 |
| Mg IV | 157 | 2,921 | 12,246 | 183,690 | 4 / 1 |
| N IV | 136 | 3,272 | 9,180 | 156,060 | 3 / 3 |
| N V | 24 | 74 | 66 | 537 | 1 / 2 |
| N VI | 49 | 212 | 273 | 4,095 | 2 / 2 |
| Ne IV | 22 | 91 | 139 | 735 | 1 / 1 |
| Ne V | 304 | 7,530 | 46,056 | 690,840 | 5 / 1 |
| O IV | 147 | 2,319 | 10,731 | 160,965 | 2 / 3 |
| O V | 166 | 2,995 | 13,695 | 205,425 | 2 / 1 |
| S V | 159 | 3,157 | 12,561 | 188,415 | 0 / 2 |
| S VI | 32 | 181 | 496 | 7,440 | 1 / 2 |
| S VII | 177 | 3,891 | 15,576 | 233,640 | 1 / 2 |

The private-file columns are discovery counts across the October 9 and October 10 validation trees; they include older unsuccessful candidates. They are **not** acceptance counts. `result.json` gives exact candidate paths and source hashes. The source reader uses CHIANTI `.elvlc`, `.wgfa`, and `.scups` for 19 ions, with separate identified Stout energy vectors for N IV and O V and other guarded valence-energy selections. Fe VII and Cl VII instead use Cloudy c25 Stout `.nrg`, `.tp`, and `.coll`. Every bound selector also uses the identified physical ionization-potential snapshot. The current source path is established by `_source_files` and the invoked `bound_excitation_atom`, not a comment alone.

Previously completed **numerical** checks include the low-density source fits for all 21 ions, but the saved products have unequal depth of finite-density and public registration checks. Direct ion-specific `result.json` and `native_result.json` paths are enumerated in the machine-readable audit. Representative actual result records are [Ar V full-bound native and registered checks](../full_bound_native/ArV/registered/result.json), [O V native cooling](../icft_energy_cooling/native/OV/native_result.json), [Cl VII source cooling](../clvii_constant_strength/products/ClVII/ClVII.json), and [S V cooling fit](../../validation_2026-10-09/stationary_cooling_exponentials/products/SV/SV.json). The O V 29-level model was compared with the 166-level source at three densities through 1e6 cm-3; O IV's older 35-level atom differs from the full 147-level model by up to 0.43245% and fails its declared 0.3% cutoff. Full-bound Ar V, Mg IV, Fe VII, and Ne V public registration has been tested with isolated objects. These checks do not establish all 21 atoms at every density or a coupled nebular calculation. No production addition is justified by private-file presence alone.

## Source temperature and radiation qualifications

The current reader explicitly holds Upsilon constant below a source or adopted lower bound for O IV at 3,200 K and O V at 5,000 K. O IV's final source article identifies 3,200–32,000,000 K. O V's 5,000 K is a source-table boundary, not a verified published accuracy limit; its current upper metadata is 50,000,000 K. Fe VII and Cl VII use the Stout 9,800–98,000,000 K nodes with a cold constant-Upsilon approximation. In their raw CHIANTI headers, Ar V states source collision strengths over 5,000–5,000,000 K, Ar VI over 7,200–7,200,000 K, C V over 5,000–50,000,000 K, N V over 5,000–50,000,000 K, and N VI over 7,200–72,000,000 K. Ne IV gives `3.6 < log T < 6.1`, approximately 3,981–1,259,000 K. Fe V's header says its spline fits are accurate over 5,000–50,000 K. The header language must not be silently promoted to a published model-accuracy range. Other ions' source interval is not established by this audit. The private 1,000–100,000 K fit grid is a numerical use range, not a source-accuracy statement. No new cold floor or high-temperature extrapolation has been selected here.

O IV has a specific unresolved radiative identity issue. The [direct diagnostic](../o_iv_main_radiative_review/README.md) found 343 of 2,064 higher positive-gf/opposite-parity rows whose stored A and gf imply a frequency more than 1% from every available observed, theoretical, or wavelength frequency. It also confirmed one original A value independently. Those facts establish inconsistent frequency metadata but do not identify which datum should be changed. The O IV full-atom serialization and numerical checks subsequently completed. Root verified that the actual cooling path uses imported A and identified energies, not gf or the stored wavelength, and adopted the explicitly qualified source-data model; see the [installation review](../o_iv_full_completion/installation/README.md). A frequency-cubed correction chosen from the wavelength alone would be unjustified. O V's energy selection resolves five mixed observed/theoretical radiative inversions through an identified 166-level Stout correspondence. Its source scattering resonances are not recalculated after the energy splice. Other ions' differing A values or source hashes are not classified here as physical defects without a transition-specific proof.

The most economical next source-temperature review is Ca VI and Fe V: they have small bound atoms, actual staged fit/native products, and no newly adopted cold rule, yet their full source temperature limits need explicit point-of-use treatment. Resolve each original collision grid and physical lower limit before adding any floor. The O IV full atom is now adopted and installed as a source-data model; detailed atomic-source accuracy remains uncertain. For the remaining ions, extract and check each actual source grid or final article interval, then validate any adopted low-temperature continuation against the 1,000–100,000 K target. Follow that with the affected compiled source, 401-by-161 registration, coupled stage/thermal checks, and all-elements terminal-band assessment. No broad test suite or production installation was run for this read-only inventory.

The audit now refuses to overwrite its original snapshot. Root verified that this refusal preserves the previous result. `post_o_iv_installation.json` was produced by the actual same reader after installation; no old private numerical product was regenerated.

## Latest oxygen installation checkpoint

The [latest directly remeasured inventory](post_oxygen_installation.json) records both O IV and O V with installed atoms and fits: 2 of the 21 additional ions, with 19 remaining additions. Both full-source controllers completed, root checked frozen hashes and full plots, and current installed files passed actual consumer checks. Imported spontaneous rates are explicitly qualified source-model estimates, not certified matrix elements. Their unused gf/wavelength metadata do not themselves establish a wrong A. See the [O IV](../o_iv_full_completion/installation/README.md) and [O V](../o_v_full_completion/installation/README.md) installation reviews. Earlier inventory snapshots are preserved. Production stage counts remain 43; no fourth/fifth oxygen stage or full 64-stage roster is active yet.
