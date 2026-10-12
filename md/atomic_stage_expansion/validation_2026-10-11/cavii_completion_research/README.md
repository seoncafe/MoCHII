# Ca VII full-scope source assessment

This assessment uses the final journal article by Landi and Bhatia (2003),
*The Astrophysical Journal*, 589, 1075-1084, DOI 10.1086/374686, the
MOCASSIN v2 `cavii.dat` atom, and the CHIANTI Ca VII files included with
Cloudy c13.05 and c17.03. It does not select or install an atom.

The MOCASSIN file has 27 active levels and 310 positive spontaneous
probabilities. Every transition pair and stored A value matches the older
Cloudy copies of `ca_7.wgfa` exactly; the c13.05 and c17.03 files are also
identical to each other at the parsed-row level. Thus the full 310-row
radiative set is locally available. The current CHIANTI file contains only
86 of these rows. The old file contains 69 rows with positive oscillator
strength, 241 with zero oscillator strength, and 102 with a negative signed
wavelength. The negative sign denotes a theoretical-energy wavelength in
this file; all four transition pairs whose MOCASSIN level-energy differences
are nonpositive have negative signed wavelengths.

The old CHIANTI header attributes its radiative data to the paper's
25-configuration calculation. A direct numerical test contradicts that
attribution for all 23 rows printed in the article's Table 3: each source A
value matches the five-configuration column within 0.1%, while every row
differs by more than 0.1% from the 25-configuration column. The 32 checked
rows in Tables 4 and 5 match the 25-configuration column within 0.1%.
The article adopts 25 configurations for radiation and five for collisions.
The other 255 radiative rows are not printed in these tables, so their
configuration provenance is unverified. The published Table 3 provides
25-configuration values for its 23 listed transitions, but a printed row
may report only one multipole whereas a source A value can be a sum.

The four nonpositive MOCASSIN energy differences are (15,17), (22,23),
(22,24), and (22,25). Their old CHIANTI wavelengths are negative, indicating
theoretical energy differences. The article says it corrected transition
rates for observed-versus-theoretical energies when both endpoints had
observations. Substituting theoretical energies for observed levels 17 and
23-25 would affect 91 radiative rows, including 73 rows currently carrying
positive signed wavelengths and 76 rows with zero oscillator strength.
The oscillator strength is insufficient to undo the energy correction for
those 76 rows: the file does not split their total A values into M1, E2, or
higher multipole components. Replacing all 27 energies with theoretical
values would affect still more corrected rows. A positive energy difference
alone does not establish a coherent unscaled radiative calculation.

The article's Table 2 provides theoretical energies and identifies the
relevant level states, but the actual energy ordering of the unobserved
levels cannot be established from the selected observed/theoretical level
vector. A physically consistent 310-transition atom requires the complete
25-configuration multipole-resolved radiative calculation or an explicitly
defined approximation with a justified energy basis and corrected rates.
The paper's printed tables and the local CHIANTI files do not provide all
the quantities needed to reconstruct such an exact atom. In particular,
replacing only Table 3's 23 A values would leave the other 255 unpublished
rows' configuration provenance unverified and would not resolve the four
energy inversions.

## Public source routes checked

The article states that its data are available as electronic files and would
be included in the next CHIANTI release. Its only specifically identified
machine-readable journal table is Table 7, which contains collision data.
The [official CHIANTI Ca VII page](https://db.chiantidatabase.org/ca/ca_7.html)
lists `elvlc`, `wgfa`, and `scups` files for the 27-level model; it does not
list multipole-resolved radiative components or the separate full
25-configuration calculation. Its current `wgfa` records the 2017 branching
filter. The [official CHIANTI wgfa format report](https://www.chiantidatabase.org/tech_reports/10_wgfa/chianti_report_10.pdf)
defines the stored columns as level indices, wavelength, weighted oscillator
strength, and total radiative decay rate. It recommends theoretical energy
differences when either endpoint lacks an observed energy and explains the
negative wavelength sign. Direct inspection of all 310 old rows agrees:
208 positive wavelengths have both observed endpoints, and 102 negative
wavelengths lack at least one observed endpoint, with no exceptions.

The [NASA NTRS article record](https://ntrs.nasa.gov/citations/20040015301)
is public but lists no downloadable attachment. An ADS API query was attempted
first for bibcode `2003ApJ...589.1075L`; this environment returned a DNS
resolution error. The previously saved
[`ADS` response](../cavii_energy_basis/ads_published_source.json) verifies
the DOI and publisher metadata. The publisher web route was blocked to the
web reader, and no supplementary multipole files were located through the
official source routes checked. This search does not establish that no such
files exist outside those routes.

The [current NIST ASD Ca VII level list](https://physics.nist.gov/cgi-bin/ASD/energy1.pl?biblio=1&conf_out=1&j_out=1&lande_out=1&level_out=1&perc_out=1&spectrum=Ca+VII&term_out=1)
does not fill the relevant observational gaps: its 27 listed levels include
four `3s2 3p 4s` states while omitting the `3s2 3p 3d` `3F` and `1D`
states at issue. Its listed `3s3p3 3P0` energy, 185356.6 cm^-1, also differs
by 100 cm^-1 from the 185256.6 cm^-1 printed in the final article and used
in the CHIANTI source. This discrepancy is recorded without selecting a
replacement value.

The source comparisons are direct local executions. No production path,
cooling table, or regression calculation was modified or regenerated.

## Main-session verification and saved reference

The main session independently rechecked every old CHIANTI transition pair
and A value against the MOCASSIN source: 310 exact matches in each of c13.05
and c17.03, including 102 negative signed wavelengths. It opened the official
CHIANTI ion page and technical report directly, downloaded the original
five-page report, and read Sections 2 and 3 using pdftotext -layout. The PDF
was saved without overwriting another reference at
`~/MoCHII/references/CHIANTI_Technical_Report_10_wgfa_2019.pdf`.
`root_review.json` records the source comparisons and reference hash.
The technical report is an official data-format document, not a journal paper.


## Additional source investigation after the goal resumed

The ADS API is now reachable through the authorized network execution route.
The saved [query result](ads_additional_sources.json) and
[candidate metadata](ads_candidate_radiation_metadata.json) identify the following
published candidates. Their methods and complete radiation coverage have
**not** been verified; titles alone do not establish an acceptable replacement.

| Candidate | Published identification | Verification state |
| --- | --- | --- |
| Tayal and Gupta (2007), *Energy levels, oscillator strengths and lifetimes in Si-like Ca VII* | Physica Scripta 75, 331; [DOI 10.1088/0031-8949/75/3/019](https://doi.org/10.1088/0031-8949/75/3/019); ADS 2007PhyS...75..331T | The publisher PDF route returned an HTML browser-validation page. The ADS original-article route returned HTTP 403. No final PDF was acquired. |
| Tayal, Gupta, and Sharma (2009), *Correlation effects on fine-structure energy levels, oscillator strengths and lifetimes in Ca VII* | Journal of Physics: Conference Series 194, 042028; [DOI 10.1088/1742-6596/194/4/042028](https://doi.org/10.1088/1742-6596/194/4/042028); ADS 2009JPhCS.194d2028T | The publisher PDF route returned an HTML browser-validation page. The ADS PDF request timed out. No final PDF was acquired. ADS classifies this conference article as not refereed. |

Access results are preserved in [the publisher record](additional_publisher_access.json)
and [the ADS PDF record](additional_ads_pdf_access.json). No access restrictions
were bypassed. In particular, the articles have not been used to justify new
energies, rates, multipole decompositions, or production data. The unresolved
checks are whether their methods provide the same 27 states, complete forbidden
radiation components, and a consistent energy basis for those components.

A directly acquired [official NORAD Ca VII data file](https://norad.astronomy.osu.edu/ca7/ca7.f.fs.txt)
attributes itself to Nahar (1993), Physica Scripta 48, 297,
[DOI 10.1088/0031-8949/48/3/008](https://doi.org/10.1088/0031-8949/48/3/008).
Its own format description specifies dipole-allowed fine-structure transitions
obtained by algebraic transformation of LS multiplets. It does not supply the
forbidden multipole components needed by the present correction. This is a
verified description of the data file, not verification of the final article's
methods. The unchanged 8,056,701-byte original is preserved at
`~/MoCHII/references/NORAD_CaVII_Nahar1993_f_fs.txt`; its SHA-256 is
`89f84bc0c39507860ed72b55ed596c9772b08cab692ef4a1cd574df713fdfd95`.

The Ca VII Stout files in local Cloudy c17.03, c22.02, c23.00, c23.01,
and c25.00 were also inspected directly. Each has 26 levels, 18 radiation
component rows for 15 distinct pairs, and no collision records. Their levels
omit five states from the supplied 27-state model, including the metastable
3F states, and add four 4s states. This is insufficient to replace the full
source model. Counts and hashes are saved in
[the additional inventory](additional_source_inventory.json).

No production atomic data or numerical cooling products were changed during
this investigation. The comparison document was corrected to describe the
already installed Mg I model, selected reference and extended atoms, default
outer-shell treatment, and remaining Ca VII/77-stage work accurately.
The next physical step is to verify the candidate final articles or obtain
complete original radiation components. Neither dipole-only NORAD data nor
the reduced Stout atom resolves the full requested model.


## Conditional radiation-frequency correction assessment

A new read-only diagnostic evaluates the published theoretical energy basis
against the original signed wavelengths for all 310 positive radiation pairs.
It does not generate or install a corrected atom. The source four inverted
radiation gaps become positive in this basis, and no other positive radiation
row acquires a nonpositive gap. The positive collision row connecting states
16 and 18 requires its lower/upper energy orientation to be reversed; the
existing converter already determines collision direction from state energies,
so no converter correction is required by this finding.

For a fixed multipole line strength, a radiation component scales as
`A_new = A_source * (nu_new/nu_source)^(2L+1)`.
Electric and magnetic dipoles have exponent 3; quadrupoles have exponent 5.
The diagnostic explicitly assumes only E1, M1, E2, and M2 components. This
assumption is conditional: the source stores total rates, not complete
component identities. The result must not be called a complete uncertainty
bound if omitted higher multipoles or configuration changes are admitted.

Using the published parity assignments and statistical weights, 111 pairs
permit a single component under that assumption and 199 permit two components.
For a sum of nonnegative components, the frequency-correction factor lies
between the smallest and largest allowed component factors. Across these
310 pairs the conditional factors range from 0.00153375758918 to
3.56296319229. This describes changes to individual radiation probabilities;
no total cooling difference, line-intensity difference, or atomic accuracy
has been measured. The values are not a single uncertainty interval for a
common multiplicative correction.

[The executed diagnostic](../../../../tools/fitting/audit_cavii_radiation_energy_scaling.py)
and [the complete 310-pair record](conditional_energy_rescaling.json) preserve
the formulas, conditional scope, original wavelengths, candidate energy gaps,
selection rules, and source hashes. This assessment makes explicit why an
unchanged total A combined with a new energy vector is not an exact frequency
correction, and which rows still require component information. It does not
resolve the previously identified five/25-configuration source differences.
The next accepted model still requires qualified original components or an
explicitly justified new approximation; neither has been adopted here.
