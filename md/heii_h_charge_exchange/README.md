# He+ + H Charge Exchange Review

Review date: 2026-09-25.

This folder collects the assessment of the singlet approach probability in
He+ + H -> He + H+, the published sources, the implementation paths, and a
reproducible comparison of the rate formulas.

**Decision (2026-09-26, applied):** the non-radiative term is the
Loreau et al. (2014) rate through microscopic reversibility, tabulated by
`tools/fitting/make_heii_h_charge_exchange.py` into
`src/heii_h_charge_exchange_data.inc` over 1e3-1e6 K and used at all
temperatures (`ion_balance_mod::nonradiative_HeII_H`, validity stated
there), with sigma_r ~ (x - x_th)^(1/2) below 15 eV/u where the quantal
data stop (sigma_f v constant at low energy). EXHALE carries the same
rate and form (9.644e-17 at 1e4 K in both codes). The
KF96 He+ row it replaces is 93x larger at 1e4 K. Above ~3e4 K
the Loreau and Zygelman calculations still disagree by up to 1e3; that
conflict is recorded at the function and in `docs/MoCHII_physics.tex`, not
resolved. EXHALE carried the same KF96 row until the same day and now
carries the same rate. Measured effect on a 1e5 K blackbody (20
iterations, not converged, with the fit carried through the threshold, 1.3x
the adopted rate at 1e4 K): x(He+) +1.3% median where He+ and H0 coexist,
lines within 0.34%; recorded in the project CLAUDE.md. The converged
Lexington runs restarted with this rate and the SLD98 radiative channel
(`results/continuous_energy/m11_hecx_solver/`, 2026-09-27) moved
<T[NpNe]> by <= 0.04 K. The notes
below are the review that led to the decision and describe the code as it
was before it.

## Contents

| File or directory | Purpose |
|---|---|
| [RADIATIVE_CHANNEL_2026-09-26.md](RADIATIVE_CHANNEL_2026-09-26.md) | GJ07 R26 (1.25e-15) against SLD98 (19) (1.20e-15) for the radiative channel: the same Zygelman x 1/4 data; against Courtney et al. (2021) Fig. 3 (digitized) SLD98 is within 2% over 1e3-1e4 K, GJ07 3-6% high; both codes now use SLD98 |
| [IGARASHI_2023_SUPPLEMENT_2026-09-25.md](IGARASHI_2023_SUPPLEMENT_2026-09-25.md) | Igarashi (2023) supplementary tables (1-100 keV/u; article not obtained): no bearing on thermal rates (Maxwell weight <= 3e-41 at 1e5 K); total capture exceeds the ground-state channel from Loreau by 650x at 0.8 keV to 3x at 80 keV, so a total capture is not the nebular channel; the endothermic n = 2 capture channels open at 5e4-1e5 K |
| [LOREAU_2014_RATE_2026-09-25.md](LOREAU_2014_RATE_2026-09-25.md) | The non-radiative rate from Loreau et al. (2014) through microscopic reversibility (`calculations/loreau2014_rate.py`), its uncertainty band, and what it means: ~1e-16 cm^3 s^-1 at 6e3-5e4 K, 3-9% of the radiative channel; with the KF96 term MoCHII's total He+ + H rate was 3.8x high at 1e4 K (as of 2026-09-25, before the decision above); relevant only where He+ meets mostly neutral H; EXHALE carried the same KF96 term then |
| [ANALYSIS_2026-09-25.md](ANALYSIS_2026-09-25.md) | Maxwell averages of Zygelman (1989) Fig. 8 and of Loreau (2014) through microscopic reversibility: the KF96 fit omits the 1/4 at 5e4-1e5 K and exceeds its own source's bound by >= 11x below 2e4 K; both calculations give ~1e-16 cm^3 s^-1 at 1e4 K |
| [REVIEW_2026-09-25.md](REVIEW_2026-09-25.md) | Physical assessment, evidence limits, numerical comparison, and the next verification steps |
| [CODE_PATHS.md](CODE_PATHS.md) | Inspected rate functions, callers, heating, relevant tests, and the Cloudy comparison (as of 2026-09-25) |
| [references/README.md](references/README.md) | Reading guide with DOI links and the exact sections, tables, or figures used |
| [references/manifest.json](references/manifest.json) | ADS metadata, original locations, page counts, and SHA-256 checksums for the eleven PDFs (ten articles and the Igarashi 2023 supplement) |
| [calculations/README.md](calculations/README.md) | What each script computes, its output, and its limits |
| [calculations/rate_comparison.py](calculations/rate_comparison.py) | Standard-library Python script that writes the 2026-09-25 formula comparison as CSV to standard output |
| [calculations/rate_comparison.csv](calculations/rate_comparison.csv) | Recorded output of that script |
| `calculations/loreau2014_rate.py`, `loreau2014_reciprocity_rate.py`, `zygelman1989_fig8_rate.py`, `radiative_channel_comparison.py` | the Maxwell averages of the ANALYSIS and LOREAU notes and the radiative-channel figure (listed in `calculations/README.md`) |

Start with the review. Use the source guide when checking a particular claim.
The PDFs are preserved without changes; no preprints are used in this folder.

## Organization and validation

The former `md/HEII_H_CHARGE_EXCHANGE_APPROACH_FACTOR_2026-09-25.md` has been
moved into this topic folder and revised as `REVIEW_2026-09-25.md`. Its
unsupported benchmark bound and unconditional statements about the radiative
fit, triplet channels, and Cloudy have been replaced with qualified statements.
The reference in `CLAUDE.md` has been updated.

The organization work of 2026-09-25 checked PDF readability and copy
checksums, local links, and execution of the formula comparison. It did not
validate a replacement physical rate, and it changed no production Fortran,
atomic data, or simulation product. The replacement came the next day (the
decision above), with its own gate
(`tests/charge_exchange/check_helium_charge_exchange.f90`).
