# Ne I and Mg I electron-excitation source validity

This audit separates numerical reproduction of CHIANTI cooling from the physical
temperature coverage of the scattering calculations. Ne I subsequently
received conditional approval as the explicit CHIANTI numerical-model
approximation described below; production installation still requires its
product checks. Mg I has not received source acceptance. The raw CHIANTI files are in the local
`RT_Codes/CHIANTI/dbase` installation; the published Mg I tables remain in the
CDS archive and are not copied into this repository.

## Ne I

The local `ne_1.scups` identifies Zatsarinny and Bartschat (2012), *Physical
Review A* **86**, 022717, DOI [10.1103/PhysRevA.86.022717](https://doi.org/10.1103/PhysRevA.86.022717), and says R. P. Dufresne prepared the CHIANTI data on January 15, 2024. The header also identifies the last-three-strength extrapolation to the Burgess–Tully infinite-temperature limit. The article's published abstract describes B-spline R-matrix calculations in jK coupling with 457 target states: 87 bound states and 370 continuum pseudostates. Its stated **incident-electron energy** coverage is excitation threshold through 300 eV. That is not a published Maxwellian-temperature accuracy interval. The APS original PDF returned HTTP 401 in this environment; no final publisher PDF was obtained. An indexed copy displaying the final journal pagination was inspected only for context, not substituted for the requested original methods verification.

Direct parsing of the local CHIANTI source gives 31 collision-connected bound levels, 142 positive radiative rows, and 465 collision rows. All collision rows use Burgess–Tully type 2 and have nonnegative ordinates. Inverting the first nonendpoint source abscissae yields 999.41–1000.36 K across rows; the last finite abscissae yield 998,531–1,001,956 K. The numerical 1,000–100,000 K interval is thus inside the CHIANTI tabulation. These transformed nodes do not independently establish the accuracy of the paper's Maxwellian rates or the details of Dufresne's conversion from energy-dependent cross sections. The maximum first finite-node temperature is 1000.355160 K; the adopted common constant-Upsilon floor is 1000.36 K, including a narrow continuation at 1,000 K. This is inferred from rounded database nodes, not a published physical validity limit. All 142 source A rows retain the physical energy orientation, and none of their source wavelengths differs by over 3% from the selected energy difference.

A private complete original-node atom was serialized at `validation_2026-10-10/source_validity_ne_i_mg_i/ne_i_candidate/nlevel_ne_1.txt`: 31 levels, 142 radiative rows, 465 collision rows, and 6,975 source nodes. The finite native source/profile controller in the adjacent validation directory completed 401 temperatures from 1,000 to 100,000 K at four electron densities (`10^-6`, `10^2`, `10^6`, and `10^10 cm^-3`). Its largest relative effective-strength difference was `6.67e-16`, largest relative cooling difference `2.45e-15`, and largest absolute population difference `4.45e-16`. Full cooling and difference plots were visually inspected: curves coincide and residuals stay near roundoff. These original 31/142 profiles are now physically stale because the subsequent source audit identified two missing spontaneous channels. They remain preserved. The `passed` flag denotes numerical comparison of that earlier model, not acceptance of its radiation completeness. The 2012 final article's detailed near-threshold treatment and the 2024 CHIANTI conversion remain the physical provenance checks before production adoption.

## Mg I

The local `mg_1.scups` identifies Barklem et al. (2017), *Astronomy & Astrophysics* **606**, A11, DOI [10.1051/0004-6361/201730864](https://doi.org/10.1051/0004-6361/201730864). Its header says the original close-coupling calculation is in LS coupling and that CHIANTI splits each LS strength using distorted-wave fine-structure ratios. It also says its Burgess–Tully infinite-temperature points were extrapolated from the last three points and warns that the resulting data were prepared for effective ionization/recombination models, not spectroscopic diagnostics.

The [publisher Version of Record held by Griffith Research Online](https://research-repository.griffith.edu.au/server/api/core/bitstreams/44dd1d66-e9e2-4f29-bf2f-156d8c096ce1/content) was inspected with `pdftotext -layout` and archived privately at `~/MoCHII/references/Barklem_etal_2017_AA_606_A11_VoR.pdf`. Its Methods and Results describe independent convergent close-coupling and B-spline R-matrix calculations and Maxwellian integration of excitation cross sections. The published effective collision strengths cover 25 LS states and 300 pairs at 500 K, 1,000–10,000 K in 1,000 K steps, 15,000 K, and 20,000 K. The [CDS J/A+A/606/A11 ReadMe](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A%2BA/606/A11?format=html&tex=true) confirms the physical-temperature matrix files, their units, and symmetric pair indexing. Neither the final paper nor CDS supplies a published 100,000 K effective-strength table. The CHIANTI source's 1,000–100,000 K numerical values cannot by themselves extend that published range.

Exactly two of the 96 bound CHIANTI collision rows have a negative source ordinate, and each negative is the extrapolated infinite-temperature endpoint: source indices 2–12, `-0.009492`, and 13–15, `-0.0007207`. Every earlier ordinate of those rows is positive. These are mathematically invalid effective strengths at the endpoint, not evidence that a negative physical collision strength was calculated. Their physical LS labels differ from the numeric source indices: the CDS matrices must be joined by the state classifications, not the same integer pair. The present strict original-node writer correctly refuses both rows. Clipping, dropping, or silently replacing their endpoints would create unverified source data.

The focused comparison script in `validation_2026-10-10/source_validity_ne_i_mg_i/compare_mg_i_cds.py` sums CHIANTI fine-structure strengths for the nine LS terms represented among the first 15 levels and compares all 36 interterm pairs with the official BSR matrices. On all 13 published temperatures, the output JSON records source SHA-256 values, all pair ratios, and every temperature's maximum and mean absolute fractional difference. The companion PDF plots the complete set of 36 ratio profiles and the largest and mean differences; it was visually inspected. At 20,000 K the mean is 0.190% and maximum 1.56%; at 1,000 K they are 4.31% and 34.0%; at 500 K they are 23.3% and 118%. The worst pair is LS terms 3–4 throughout. These are measured differences between CHIANTI's split/interpolated representation and the published BSR source, not estimates of the physical CCC–BSR uncertainty. The exact cause of the low-temperature discrepancy requires examination of the CHIANTI preparer's conversion and source interpolation; agreement at high temperature does not establish fidelity at the cold end.

The directly published LS matrices are a coherent nonnegative alternative for **500–20,000 K** after a deliberate LS-level radiation/energy model is established. Simply substituting LS matrices into the current 15-level fine-structure atom would mix incompatible level definitions. Nothing in the final paper or CDS confirms 20,000–100,000 K collision strengths or provides a correction for the two negative CHIANTI infinite-temperature endpoints. Mg I therefore remains a numerical-fit candidate, not a physically validated 1,000–100,000 K production atom.

## Verification and limits

The local CHIANTI files and the official CDS BSR matrix bytes were read directly. The latter were kept only under `/tmp/mgi_cds` for this comparison; the repository contains code, hashes, and results, not their raw numerical rows. No new package was installed. The prior cooling-fit residuals test representation of the selected CHIANTI calculation and do not alter the source-validity findings here.

The main session independently checked the completed Ne I controller state,
all frozen source hashes, and the result artifact hash. It inspected the full
four-density cooling profiles and recorded the review in
[`main_review.json`](validation_2026-10-10/source_validity_ne_i_mg_i/main_review.json).
The maximum relative cooling difference measured in that result is
`2.4424906541753444e-15`. That review concerned the original 31/142 model; its profiles cannot validate the subsequently corrected 31/144 model.

Visual review also found that the Mg I comparison figure's original long
vertical label extended outside the page. The label now states the same ratio
with standard collision-strength notation. Only the figure was regenerated
from its preserved measurement JSON, whose bytes were checked unchanged; no
collision strengths or native profiles were recalculated. The corrected PDF
was rendered and visually checked.

Further final-article access attempts and the precise outstanding source
claims are recorded in the [access audit](source_access_final_articles.md).

## Subsequent Ne I source-model judgment and radiation correction

The official CHIANTI collision representation is adopted as an explicit
numerical approximation over 1,000--100,000 K. Its 465 collision rows and
6,975 source nodes are preserved. Below the common inferred finite-node
floor of 1000.36 K, Upsilon is held constant while the actual gas temperature
enters the Maxwellian factors. This does not verify the unavailable final
article Methods or the 2024 conversion; atomic accuracy remains unquantified.
Local Stout has no independent Ne I electron-collision table and local
MOCASSIN v2 has no Ne I collision file.

The subsequent radiation audit confirmed that CHIANTI source level 15,
3p 1S0 at 152970.733 cm^-1 with statistical weight 1, incorrectly had no
downward spontaneous channel. The first 31 Stout states match the observed
CHIANTI energies and weights. The NIST-dated Stout spontaneous table supplies
exactly two missing transitions within those states:

| Lower source level | Upper source level | Spontaneous probability (s^-1) |
| --- | --- | --- |
| 3 (3s, J=1) | 15 (3p, J=0) | 9.00e5 |
| 5 (3s, J=1) | 15 (3p, J=0) | 6.82e7 |

These two rates are selectively added. All 142 existing CHIANTI A values,
levels, and collision nodes remain unchanged; Stout does not replace the
remaining radiation table. Exact source hashes, configuration/J/energy
identity, and the added rates are checked at the actual reader. The corrected
31/144 atom has only the ground state as a radiatively closed state. This
established radiation defect must be corrected regardless of its cooling
effect size. Earlier 31/142 profiles and fits are preserved as physically
stale products.

The corrected private atom passed the root session's focused structural
check. Full native versus independent raw-source profiles, population and
energy balances, cooling coefficients, the compiled coefficient reader,
and the registered density table are the next checks before installation.
The prepared finite controllers execute those authorized successors without
requiring another user message. The corrected Ne I production atom and fit are now installed with exact
validated bytes. The finite source and product controllers completed all
steps; root checked all 43 and 46 frozen hashes, respectively, and visually
inspected complete cooling, population, residual, and fit curves. Maximum
source/native cooling error is 2.22e-15; the 801-temperature fit error is
0.0317686%; the actual GNU coefficient-reader error is 7.78e-12. The
64,561-cell Intel registry passed with maximum sampled application error
0.0041481%. A production debug reader opened the installed atom at 10,000 K
and four densities; all population and cooling records exactly matched the
validated profile. See [the root review](validation_2026-10-10/ne_i_radiative_completion/main_review.json).
