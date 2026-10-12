# New Atomic Expansion Reference Storage

Newly obtained papers and their associated author data are stored centrally
in `/nfs/mocafe/kiseon/MoCHII/references/` (`~/MoCHII/references/`). The local
`md/atomic_stage_expansion/references/` entries are verified symbolic links,
so existing relative document links and historical source hashes remain usable.
Path relocation did not regenerate atomic data, figures, or model outputs.

| File | Source and verification | SHA-256 |
|---|---|---|
| `Kingdon_Ferland_1999_ApJ_516_L107.pdf` | Final publisher PDF; methods and Table 1 read and table visually checked; DOI 10.1086/312008 | `2c681fe32532baddee5d19a67e472cefc51f6a0135de01ab82eb9aa182e78def` |
| `Ferland_Korista_Verner_Dalgarno_1997_ApJ_481_L115.pdf` | Harvard DASH authorized final-journal copy with repository cover; final ApJ L115--L117 and Section 2 verified; DOI 10.1086/310655 | `612bf229d498dd20e4019dd8440282373af410237c2c7a0d4a6554ae68cd8af0` |
| `Zhang_2026_Fe3plus_to_Fe6plus_supplement.zip` | Official author Zenodo dataset, DOI 10.5281/zenodo.20264991; ground-parent ADF09 DR tables read directly. Associated final publisher article was supplied by the user and verified. | `ed44a5f76b824ec462c6807e3566d4bb82ebb0e4064e5856e398882a5f9d0033` |
| `Zhang_2026_ApJS_285_50.pdf` | User-supplied final ApJS 285:50; Tables 6/7 and Eqs. 4--6 verified; DOI 10.3847/1538-4365/ae80c0 | `ba91543363d26b11593c19ade6c731484a7ed8b39da7679493399c35cfa40d50` |
| `Fernandez-Menchero_2017JPhysB_50_065203.pdf` | User-supplied final J. Phys. B 50:065203; methods and Tables 1/2 verified; DOI 10.1088/1361-6455/aa5fc4 | `0bcf693a59608c5236082ecb4e5f6e35bd8d6a8e3f98e530851794886cc8de3e` |

The Fe supplemental archive is data rather than a journal PDF. The user
supplied both final journal PDFs on October 10, 2026; they were read directly
with `pdftotext -layout`. Methods, equations, tables, and local atomic-data
identities were checked. See [Fe verification](zhang2026_final_paper_verification.md)
and [N IV verification](niv_final_paper_verification.md). Earlier publisher
access failures no longer block these two checks. No new literature search
was needed to read the supplied files.

## Additional final journal articles for unresolved cooling energies

The ADS API was attempted first; DNS resolution failed. The A&A publisher PDF
URLs returned HTTP 403. Final published copies were then obtained from the
author-affiliated University of Strathclyde atomic-data archive. Their printed
A&A identifiers, DOI, and EDP Sciences publication headers were checked, and
Sections 2 (Structure) and 3 (Scattering) were read with `pdftotext -layout`.
No atomic indices have been replaced based on a paper's level count alone.

| Central file | Final journal DOI | Author-affiliated final PDF | SHA-256 |
|---|---|---|---|
| `Fernandez-Menchero_2014_AA_566_A104.pdf` | 10.1051/0004-6361/201423864 | https://amdpp.phys.strath.ac.uk/tamoc/PUB/2014/aa_566_A104.pdf | `a893e36ea05605c83b6dff0d30498af483a5c193dc45de44460ea266f573ac79` |
| `Fernandez-Menchero_2014_AA_572_A115.pdf` | 10.1051/0004-6361/201424849 | https://amdpp.phys.strath.ac.uk/tamoc/PUB/2014/aa_572_A115.pdf | `29afcae994beeaf3dbdf47079ff582faf59b4feedee6ddccd827efdece9bbde7` |

The Be-like article calculates a 238-level ICFT model; the Mg-like article
calculates a 283-level ICFT model. The local O V and Cl VI CHIANTI files have
different model sizes. A verified mapping of level identities and transition
indices is required before using either corresponding Stout energy vector.
The methods use calculated target energies for scattering; replacing threshold
energies does not establish low-temperature resonance accuracy.

## Fe VII final-source verification

Fe VII's final journal article was obtained from the author-affiliated server
and saved as `~/MoCHII/references/Witthoeft_Badnell_2008_AA_481_543.pdf`.
Its SHA-256 is
`f266f23796567ee0a6da0305cc3b1c3430083396109adf0c2eb81113de90ae7a`.
See [Fe VII verification](fevii_final_paper_verification.md) for final
methods, all 189 checked energy rows, 57 printed radiative rows, and the
remaining collision-data and cold-temperature limitations.

## Cl VII final-source verification

ADS was attempted first and DNS resolution failed. A&A returned HTTP 403;
IOP returned a non-PDF access page. Final journal PDFs were obtained from the
author-affiliated Strathclyde archive and read using `pdftotext -layout`.
The outer paper Sections 2.1 and 2.3 and inner paper Sections 2.1 and 2.3
were read. Their temperature interval is 200(q+1)^2--2e6(q+1)^2 K, where
q=Z-11. For Cl VII q=6, so the lowest supplied temperature is 9800 K.

| Central final PDF | DOI | SHA-256 |
|---|---|---|
| `Liang_Whiteford_Badnell_2009_AA_500_1263.pdf` | 10.1051/0004-6361/200911866 | `9ba98f031d4cceae3b9f1c6096d9dd228f65f6cd26f23144168c91fe2bfc20df` |
| `Liang_Whiteford_Badnell_2009_JPhysB_42_225002.pdf` | 10.1088/0953-4075/42/22/225002 | `ba6c4af9be8377d04fe96222d80b5b9d33565c6bd8789fd78a749dfa53b38bf4` |

Both PDF download URLs are under `https://amdpp.phys.strath.ac.uk/tamoc/PUB/2009/`.
The outer paper has 32 levels and the inner paper has 134, overlapping in
five low levels. The 161-level Stout configuration set is consistent with
combining these models, but the table header cites only the outer paper;
complete numeric provenance of the core-excited collision rows is not verified.
Do not infer a bound-state atom from the total level count.

## Further excitation and reverse-rate references

Nine final journal PDFs used to check the terminal reverse-rate sources are
now in the central reference directory. Their filenames, measured byte counts,
and hashes are recorded in the
[reference archive](validation_2026-10-10/terminal_reverse_provenance/reference_archive.json).
The [terminal reverse-rate audit](terminal_reverse_rate_source_validity.md)
distinguishes the published fit intervals from the current author coefficient
sets and sampled fit residuals. Availability of these papers does not establish
terminal ionization closure.

The final Storey and Zeippen (2000), Froese Fischer (2006), Wiese and Fuhr
(2007), and Froese Fischer and Tachiev (2004) radiative compilations have also
been preserved centrally. The
[radiative-source audit](oxygen_radiative_source_provenance.md) records their
filenames and hashes, the O III optical correction, and why differences in
other ions' A values do not by themselves establish an incorrect operator.

`Barklem_etal_2017_AA_606_A11_VoR.pdf` is the final published Mg I article
obtained from Griffith Research Online. Its methods and associated official
CDS temperature tables were checked as documented in the
[Mg I source audit](source_validity_ne_i_mg_i.md).
Its published 500--20,000 K matrices do not verify a 100,000 K continuation.

The final Ne I and Cl V articles remain unavailable through the tested
publisher and institutional routes. The
[access log](source_access_final_articles.md) identifies both complete
citations, access results, and the exact claims that still require their
original final PDFs. No substitute preprint or indexed copy has been archived
as either final article.
