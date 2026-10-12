# Ca VII verification against the final Landi-Bhatia article

## Physical finding

The supplied final article resolves the publication-access dependency and
confirms the separate 25-configuration radiative calculation and
5-configuration collision calculation. It also reveals a source-data
inconsistency that the source file's attribution alone did not establish.
The MOCASSIN Ca VII file combines radiation rates that numerically match
different configuration models in the published tables.

The checked final version is Landi and Bhatia (2003), *The Astrophysical
Journal*, 589, 1075-1084, DOI 10.1086/374686. The local PDF has ten pages,
the journal masthead, received/accepted dates, and the matching title and
authors. Its identity agrees with the previously verified ADS record.
The methods in Sections 2.1 and 2.2, the reliability discussion in Section
2.3, and Tables 2-7 were read directly from the final article using
`pdftotext -layout`. Tables 2-4 were also visually checked in the PDF.

## Measured radiation comparison

| Published table | Checked pairs | Column matching the source A values |
| --- | ---: | --- |
| Table 3: ground forbidden and metastable transitions | 23 | Five configurations |
| Table 4: allowed ground-to-3s3p3 transitions | 17 | Twenty-five configurations |
| Table 5: allowed ground-to-3s2 3p3d transitions | 15 | Twenty-five configurations |

All 55 comparisons match the stated column within 0.1%, allowing printed
rounding and possible summed multipole contributions in source records.
This is a numeric comparison of specific rows, not verification of all 310
positive source radiation pairs. The other 255 pairs' configuration
provenance remains unverified. The article adopts the 25-configuration
radiative calculation, with five-configuration results shown for comparison.

Examples make the distinction explicit:

| Pair | Source A [s^-1] | Published five-configuration A | Published twenty-five-configuration A |
| --- | ---: | ---: | ---: |
| 1-4 | 2.717e-4 | 2.717e-4 | 1.072e-3 |
| 2-5 | 34.74 | 34.74 | 56.85 |
| 6-16 | 0.4051 | 0.4051 | 0.1732 |
| 9-16 | 1.056 | 1.056 | 2.183 |
| 1-7 | 3.092e8 | 3.038e8 | 3.092e8 |
| 4-27 | 3.909e7 | 5.262e9 | 3.909e7 |

Table 3 includes an explicit note that the 9-16 value contains an added M1
contribution. For 2-3 and 3-4 the source values differ slightly from the
printed five-configuration M1 entries; the printed entries alone do not
establish whether this is rounding or additional multipoles. Do not silently
replace a total A by a single multipole without checking that distinction.
The source header attributing the rates generally to the 25-configuration
calculation is insufficient evidence for every row. External scientific
source files are preserved for direct comparison and were not edited.

## Energy basis and collision calculation

Section 2.1 states that radiation rates and oscillator strengths were corrected
for observed-versus-theoretical energy differences when observations exist
for both transition endpoints. An all-theoretical level vector with unchanged
rates cannot therefore be labeled a consistently rescaled theoretical atom.
The published and current CHIANTI theoretical level vectors also differ:
the largest measured difference is 20.0 cm^-1, at level 10. The vectors have
the same listed state identities, but they are not byte-identical numerical
realizations of Table 2.

The four previously measured nonpositive source photon gaps remain a physical
inconsistency. Their disappearance with a theoretical-only energy vector is
a numerical observation, not sufficient acceptance of that replacement.
The observed energies should not be discarded without evaluating a coherent
choice of radiation probabilities, multipole scaling, and collision thresholds.

Section 2.2 confirms distorted-wave collisions calculated with five
configurations at incident electron energies of 8, 10, 15, 20, 30, 40, and
60 ryd. The results were scaled, fitted using Burgess-Tully methods, and
integrated over a Maxwellian distribution. Table 7 belongs to that collision
model. These seven incident energies are not the MOCASSIN file's 23 Maxwellian
temperature nodes. The article's warning in Section 2.4 concerns low-temperature
extrapolation relative to the ground-term close-coupling results of Galavis
et al. (1995); a cold continuation must remain an explicit approximation.

The printed Table 2 note calls level 18 a 1P1 state, while its row and the
surrounding discussion identify 3S1. This published labeling discrepancy was
recorded; the original PDF was preserved unchanged. The numerical source
mapping uses the actual row's state identity, not that note.

## Acceptance and remaining work

[published_source_review.json](published_source_review.json) records all 55
pairs, the 27 printed theoretical energies, and hashes of the article,
source file, CHIANTI energies, and audit program.
[root_review.json](root_review.json) records direct inspection.
The first audit attempt used an incorrect CHIANTI fixed-width slice; its
failed program is preserved separately. The corrected audit uses the documented
72:87 theoretical-energy field and passes all comparison assertions.

No corrected atom was accepted or installed, and no rate equation was changed.
The source review must resolve the radiation model and nonpositive energy
pairs together before a candidate can pass the native reader. The user subsequently approved the public-input and paired-data changes on
October 11, 2026; the physical Ca VII correction remains unfinished. The
article-access dependency itself is resolved; do not request this paper again.
