# Osorio et al. (2019) Ca I model-atom source audit

Date: October 10, 2026. This is a source-access audit, not an adopted Ca I cooling atom. No source rows or publication files were copied into the repository.

## Publication and physical scope

NASA ADS identifies Osorio, Lind, Barklem, Allende Prieto, and Zatsarinny, [*Ca line formation in late-type stellar atmospheres. I. The model atom*, A&A **623**, A103 (2019)](https://doi.org/10.1051/0004-6361/201834680), bibcode `2019A&A...623A.103O`. The ADS record marks an open publisher PDF but contains only `SIMBAD:2` in its `data` field; it exposes no model-atom download. Its `ASSOCIATED` property does not identify a numerical archive.

The indexed text of the final publisher PDF describes a 127-level non-LTE atom: 96 Ca I levels, 30 Ca II levels, and the Ca III ground level. It lists 1,656 Ca I and 152 Ca II bound-bound transitions. Some Ca I fine-structure components remain separate while higher triplets and Rydberg states are merged. Electron excitation for the lower Ca I states through `4s.8s(1S)` uses an extension of the B-spline R-matrix calculation of Zatsarinny et al. (2006). Higher states use an impact-parameter approximation; Ca II electron excitation uses another source. The paper also includes electron ionization, hydrogen excitation, and charge exchange. These mixed processes and merged states cannot be treated as a complete fine-structure electron-strength table simply by copying the reported model size. The indexed text does not establish a tabulated temperature interval for every electron rate.

A note in the paper says its model atom **will be available** from A&A and a future model-atom collection. This is a promise of distribution, not a verified download. The final A&A PDF and HTML currently return HTTP 403 from this environment, so the complete final methods, any later corrections, and the numerical model have not been inspected directly. The publicly indexed excerpts above are provisional evidence, not a substitute for reading the complete final PDF.

## Access checks

| Resource | URL | Observed response and content |
| --- | --- | --- |
| ADS article record | `https://api.adsabs.harvard.edu/v1/search/query?q=bibcode:2019A%26A...623A.103O` | HTTP 200 JSON; official PDF link, `data=["SIMBAD:2"]`, no numerical model link |
| Final publisher PDF | `https://www.aanda.org/10.1051/0004-6361/201834680/pdf` | HTTP 403, 93-byte HTML access-denial page; no PDF saved |
| Final publisher HTML | `https://www.aanda.org/articles/aa/full_html/2019/03/aa34680-18/aa34680-18.html` | HTTP 403, same 93-byte HTML access-denial page |
| Candidate publisher online material | `https://www.aanda.org/articles/aa/olm/2019/03/aa34680-18/aa34680-18-olm.zip` | HTTP 403, same denial page; this does not establish that a ZIP exists |
| Candidate CDS catalogue | `https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A103/ReadMe` | HTTP 404 HTML; no catalogue established |
| Local reference search | `~/MoCHII` and `~/RT_Codes`, filenames matching Osorio/Ca atom/collision terms | No matching local Osorio model file found |

The publisher's 403 response body SHA-256 begins `0d3e98ca727fc120`; the CDS 404 response body SHA-256 begins `f5a8d05b5522789f`. These hashes identify error pages, not scientific data. ADS and the indexed final-article text were checked before wider web searches. Searches for an author or institutional model archive, CDS record, and public model-atom collection produced no numerical file tied to this article. No new PDF was obtained for `~/MoCHII/references`.

## Requirement for a usable Ca I source

Obtain the actual electronic Osorio model with its redistribution conditions, license, state identifiers, energies, statistical weights, radiative probabilities, electron collision-rate tables or cross sections, temperature grid, and distinction between BSR-derived and approximate rows. Then verify pair completeness and radiative closure for a coherent bound Ca I subset. Any merged-term to fine-structure conversion requires a physical recoupling rule, not an index join. The 2019 Zatsarinny full cross-section set remains separately unresolved as recorded in [the prior access audit](../fe_ca_source_completion/README.md). Neither source currently supports a production Ca I cooling table.
