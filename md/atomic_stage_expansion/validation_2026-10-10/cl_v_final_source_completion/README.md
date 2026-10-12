# Cl V final-source access and numerical identity audit

## Physical judgment

The MOCASSIN v2 `clv.dat` collision table is consistent, at ten common
temperatures, with the two-decimal Cl V effective collision strengths exposed
by a search index of the typeset Saraph and Storey (1999) paper. This is strong
evidence that its low-temperature entries derive from the calculation reported
there. It is **not** direct verification of the published final article. The
publisher PDF and full final methods and table remain inaccessible here, and
the article must not yet be accepted as proof for a 1,000–100,000 K production
Cl V cooling model. The MOCASSIN nodes at 50,000, 100,000, and 200,000 K lie
beyond the indexed printed table ending at 40,000 K; their provenance and
accuracy remain unverified. This is a ground-term **two-level** model and omits
cooling through higher Cl V states.

Saraph and Storey, *Astronomy and Astrophysics Supplement Series* **134**,
369–375 (1999), DOI [10.1051/aas:1999437](https://doi.org/10.1051/aas:1999437),
ADS bibcode `1999A&AS..134..369S`. The ADS API directly confirmed the title,
authors, journal, year, DOI, and an abstract identifying the Cl V ground-term
fine-structure transition and resonance-sensitive collision calculation. The
ADS record marks the article open access, although its linked publisher PDF
returned HTTP 403 in this environment.

## Access routes checked on 2026-10-10

| Route | Result |
| --- | --- |
| NASA ADS API bibcode lookup | Successful metadata and abstract; no full article text in response. |
| DOI and `aas.aanda.org/10.1051/aas:1999437/pdf` | HTTP 403. |
| Publisher `aas.aanda.org/articles/aas/pdf/1999/01/ds1586.pdf`, with and without download query | HTTP 403. |
| Publisher HTML under `aas.aanda.org/articles/aas/full/1999/01/ds1586/` | HTTP 403. |
| `www.aanda.org/articles/aas/pdf/1999/01/ds1586.pdf` | HTTP 403. |
| ADS `articles.adsabs.harvard.edu/pdf/1999A%26AS..134..369S` | Timed out. |
| ADS `articles.adsabs.harvard.edu/full/1999A%26AS..134..369S` | Redirected to a scan-manifest web application; no readable article returned. |
| ADS `link_gateway/.../PUB_PDF` | HTTP 403. |
| CDS `J/A+AS/134/369` | CDS responded that the catalog is not found or unavailable; no numerical supplement. |
| CiteSeerX indexed, journal-formatted PDF URL | Search index exposed page text, but direct PDF GET returned HTTP 404. |
| ResearchGate author-uploaded, journal-formatted text | Search index exposed excerpts; the actual final file was not obtained. |

No file was added to `~/MoCHII/references/` because no legitimate final PDF
could be downloaded. The journal-formatted indexed text appears to reproduce
pages 369–375, but it is not a substitute for inspecting the final article.
The exact methods, table transcription, and temperature-validity statements
still require direct final-article verification. The article's 21-term
close-coupling description and a possible discussion of reliability through
30,000 K are **indexed claims only**, not accepted source facts.

## Direct local data and provisional comparison

The external original file is
`/home/kiseon/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/data/clv.dat`,
SHA-256 `819052a750275325438e384bed46649b69627b2fc04d86125161d81ded415a52`.
Direct parsing gives two levels, one collision pair, 15 Upsilon nodes from
2,000 to 200,000 K, lower and upper statistical weights 2 and 4, gap
1490.8 cm^-1, and spontaneous rate 0.0298 s^-1. Its header attributes the
collision strengths to Saraph and Storey (1996, in preparation), the rate to
unpublished Storey data, and the energy to a NIST compilation. The local
MOCASSIN source invokes this file through `makeElements` and `forLines` /
`equilibrium`; these are live data rather than an orphan table. The spontaneous
rate and gap are not supplied by the collision article, so that paper cannot
validate them.

The [indexed Table 5 excerpt](https://citeseerx.ist.psu.edu/document?doi=15e037d5725168fa42aab91856e8d82d4ab83d20&repid=rep1&type=pdf)
exposes 25 complete Cl V rows from 1,000 to 36,000 K. It begins a 40,000 K
row, but the Cl V value was cut off in the available snippet. Ten printed
temperatures coincide exactly with the MOCASSIN nodes. Their absolute
differences in Upsilon are at most 0.003, within rounding to the table's two
decimal places. For example, at 2,000/6,000/20,000 K the local values are
2.718/3.487/4.427 and the indexed values are 2.72/3.49/4.43. The complete
machine-readable local result is [`result.json`](result.json); the parsing and
comparison procedure is [`compare_indexed_table.py`](compare_indexed_table.py).
These are numeric comparisons with an **unverified search-index transcription**,
not published-table validation.

The user-approved cold approximation of constant Upsilon below the local first
node gives 2.718 at 1,000 K. The indexed table gives 2.59, a 4.9421% larger
local Upsilon. The source is therefore *not* numerically identical to the
indexed calculation at 1,000 K even though the adopted cold continuation is
explicit and bounded. The Maxwellian Boltzmann factor must use actual gas
temperature, not the clamped source temperature. If the final printed 1,000 K
value is confirmed, it provides a direct source value for replacing this
particular approximation after a physical review.

The exact remaining checks are the final paper's full Section 3 and Section 4
method, complete Table 5 including its last Cl V entry, any stated reliability
or temperature bounds, and whether the MOCASSIN values above 40,000 K are
published elsewhere. A final PDF from the publisher or a verified copy of the
published pages is needed. No cooling fit, source atom, or production data was
generated in this audit.
