# Update: final article obtained and directly reviewed

The user supplied the final Landi and Bhatia (2003) PDF. The previous access
attempts below are historical. Methods and printed tables are now verified
against the final version. The review also finds mixed configuration-model
radiation data: 23 Table 3 rows match five-configuration values, while 32
Tables 4-5 rows match twenty-five-configuration values. A theoretical-only
energy substitution is not yet an accepted corrected atom. See
[the final-publication review](../cavii_published_review_checked/README.md).

# Ca VII observed and theoretical energy comparison

The local source-table comparison narrows the nonpositive-photon-energy
problem to mixed observed and theoretical energies. This is an audit of data
and definitions, not acceptance of a replacement physical model.

All 27 statistical weights agree between the MOCASSIN file and current
CHIANTI. The original energies differ from CHIANTI's observed-if-available,
otherwise-theoretical selection by at most 0.031 cm^-1. The supplied MOCASSIN
atom has 310 positive radiative pairs. Current CHIANTI retains 86; all 86 are
present in the original and their A values agree exactly at stored precision.
Thus simply replacing 310 rates with current CHIANTI would omit 224 source
pairs, rather than prove full source coverage.

| Pair (lower, upper) | Source A (s^-1) | Mixed-energy gap (cm^-1) | Theoretical gap (cm^-1) |
| --- | ---: | ---: | ---: |
| (15, 17) | 5.723e-09 | -331.200000 | 6759.719000 |
| (22, 23) | 2.202e-05 | -1672.700000 | 5029.531000 |
| (22, 24) | 0.01185 | -1038.700000 | 5724.812000 |
| (22, 25) | 0.007001 | -678.700000 | 6187.969000 |

All four problematic pairs are absent from the current 86-pair table. Every
original positive A has a positive gap when compared with the all-theoretical
energies. This establishes numerical consistency of that energy basis; it
does not establish its absolute atomic accuracy, justify neglecting observed
energies everywhere, or independently verify the direction and normalization
of each transition. No A values or level energies have been installed or
accepted on this evidence alone. Original external source files remain intact.

[The full comparison](source_comparison.json) includes all 27 energy rows and
source hashes. ADS identifies the cited final article as Landi and Bhatia
(2003), *The Astrophysical Journal*, 589, 1075, “Atomic Data and Emission-Line
Intensities for Ca VII,” DOI 10.1086/374686, bibcode 2003ApJ...589.1075L.
This citation metadata is verified by [the ADS response](ads_published_source.json).
Crossref's [publisher metadata](crossref_metadata.json) additionally verifies
the full page range, 1075-1084. The methods, tables, theoretical energy model,
and A-value normalization have not yet been verified against the final article.

## Article access and dependent work

Publisher PDF access through IOP returned HTML rather than PDF. The original
Chicago publisher PDF route returned HTTP 403. The legitimate ADS article
PDF request timed out, and the browser could not read either IOP or ADS.
[Access attempts](article_access.json) record these results. No access
restrictions were bypassed and no preprint or abstract was substituted.

Additional legitimate requests are preserved in
[the second access record](article_access_second_attempt.json) and
[the Crossref-link access record](article_access_crossref_links.json).
The ADS PDF request again timed out. Crossref's HTTP and HTTPS final-version
PDF links redirected to the publisher's HTML access challenge, not a PDF.
Challenge query parameters were omitted from the record. Citation metadata
does not establish the missing atomic-model methods.

A final journal article is still required to verify the use of the
25-configuration radiative model with the 5-configuration collision model,
the theoretical energy ordering, and any required A-value energy correction.
These physical choices block accepting a corrected Ca VII atom, but do not
block independent options or other ions. The measured source comparison is
retained as evidence for the next source review.
