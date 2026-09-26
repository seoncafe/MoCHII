# He+ + H Charge Exchange Review

Review date: 2026-09-25.

This folder collects the assessment of the singlet approach probability in
He+ + H -> He + H+, the published sources, the implementation paths, and a
reproducible comparison of the rate formulas.

**Current decision:** do not apply an additional factor of 1/4 to the KF96
term without establishing its normalization. The radiative correction is
documented. The normalization and temperature dependence of the term that
MoCHII treats as nonradiative remain unresolved.

## Contents

| File or directory | Purpose |
|---|---|
| [IGARASHI_2023_SUPPLEMENT_2026-09-25.md](IGARASHI_2023_SUPPLEMENT_2026-09-25.md) | Igarashi (2023) supplementary tables (1-100 keV/u; article not obtained): no bearing on thermal rates (Maxwell weight <= 3e-41 at 1e5 K); total capture exceeds the ground-state channel from Loreau by 650x at 0.8 keV to 3x at 80 keV, so a total capture is not the nebular channel; the endothermic n = 2 capture channels open at 5e4-1e5 K |
| [LOREAU_2014_RATE_2026-09-25.md](LOREAU_2014_RATE_2026-09-25.md) | The non-radiative rate from Loreau et al. (2014) through microscopic reversibility (`calculations/loreau2014_rate.py`), its uncertainty band, and what it means: ~1e-16 cm^3 s^-1 at 6e3-5e4 K, 3-9% of the radiative channel; MoCHII's total He+ + H rate is 3.8x high at 1e4 K; relevant only where He+ meets mostly neutral H; the same KF96 term is in EXHALE |
| [ANALYSIS_2026-09-25.md](ANALYSIS_2026-09-25.md) | Maxwell averages of Zygelman (1989) Fig. 8 and of Loreau (2014) through microscopic reversibility: the KF96 fit omits the 1/4 at 5e4-1e5 K and exceeds its own source's bound by >= 11x below 2e4 K; both calculations give ~1e-16 cm^3 s^-1 at 1e4 K |
| [REVIEW_2026-09-25.md](REVIEW_2026-09-25.md) | Physical assessment, evidence limits, numerical comparison, and the next verification steps |
| [CODE_PATHS.md](CODE_PATHS.md) | Inspected rate functions, callers, heating, relevant tests, and the Cloudy comparison |
| [references/README.md](references/README.md) | Reading guide with DOI links and the exact sections, tables, or figures used |
| [references/manifest.json](references/manifest.json) | ADS metadata, original locations, page counts, and SHA-256 checksums for ten published PDFs |
| [calculations/README.md](calculations/README.md) | Scope and limitations of the formula evaluation |
| [calculations/rate_comparison.py](calculations/rate_comparison.py) | Standard-library Python script that writes the comparison as CSV to standard output |
| [calculations/rate_comparison.csv](calculations/rate_comparison.csv) | Recorded output of that script |

Start with the review. Use the source guide when checking a particular claim.
The PDFs are preserved without changes; no preprints are used in this folder.

## Organization and validation

The former `md/HEII_H_CHARGE_EXCHANGE_APPROACH_FACTOR_2026-09-25.md` has been
moved into this topic folder and revised as `REVIEW_2026-09-25.md`. Its
unsupported benchmark bound and unconditional statements about the radiative
fit, triplet channels, and Cloudy have been replaced with qualified statements.
The reference in `CLAUDE.md` has been updated.

The organization work checks PDF readability and copy checksums, local links,
and execution of the formula comparison. It does not validate a replacement
physical rate. Production Fortran, atomic data, and existing simulation
products are unchanged. No simulation or Fortran regression test was run for
this documentation and reference organization.
