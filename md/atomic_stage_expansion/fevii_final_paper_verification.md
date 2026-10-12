# Fe VII final journal source verification

The final typeset article by M. C. Witthoeft and N. R. Badnell (2008),
*Astronomy & Astrophysics* 481, 543--551,
[DOI: 10.1051/0004-6361:20079073](https://doi.org/10.1051/0004-6361:20079073),
was obtained from the [author's public server](https://amdpp.phys.strath.ac.uk/tamoc/PUB/2008/aa_481_543.pdf).
Its journal identification, pages, and first-page appearance were checked;
methods and results were read with `pdftotext -layout`. The PDF is stored at
`~/MoCHII/references/Witthoeft_Badnell_2008_AA_481_543.pdf` (602,821 bytes;
SHA-256 `f266f23796567ee0a6da0305cc3b1c3430083396109adf0c2eb81113de90ae7a`).
Prior ADS metadata verification identified bibcode `2008A&A...481..543W`.
Publisher and ADS PDF routes had failed; the author-hosted copy contains the
final journal article rather than a preprint.

## Published method and limits

Section 2 describes ICFT R-matrix scattering with 189 close-coupling levels
from 89 terms, using a larger configuration-interaction structure calculation.
Collision strengths are Maxwellian averaged. Sections 3.1--3.2 discuss
resonance-resolution tests and sensitivity to configuration interaction;
the tabulated temperature interval is 9,800--98,000,000 K. Section 4 notes
additional configurations beyond the calculation's capability. These tests
do not establish a uniform accuracy bound for every transition or validate
constant strengths below the lowest source temperature.

## Direct correspondence with installed source files

The [reproducible table check](validation_2026-10-10/fevii_final_paper/check_published_tables.py)
extracts all 189 Table 3 entries and compares them with the unchanged Stout
`fe_7.nrg`. Levels 1--9 match the printed NIST column within
4.735e-7 Ry. Levels 10--189 match the printed calculation column within
5.100e-5 Ry, consistent with the journal's four-decimal printing and the
Stout energy rounding. Thus the adopted source uses observed ground-complex
energies and calculated higher energies; it does not use calculated energies
for every level.

All 27 Table 4 radiative rates and all 30 Table 5 rates from the
40-configuration calculation match Stout `fe_7.tp` at the printed precision.
The [complete numerical comparison](validation_2026-10-10/fevii_final_paper/level_correspondence.json)
retains every checked row. The Stout collision file attributes its 189-level
ICFT dataset to this article and has the same temperature endpoints. The
article does not print all 17,766 collision-strength curves, so this check
does not independently authenticate every deposited collision number.

The source reader and existing cutoff diagnostics remain applicable.
The user-approved constant Upsilon continuation at 1,000--9,800 K remains
an explicit approximation with unverified resonance accuracy. No production
stage was activated and no spatial plasma model was run during this audit.
