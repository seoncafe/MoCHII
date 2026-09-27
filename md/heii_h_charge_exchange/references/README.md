# Published Sources

The eleven PDFs below (ten articles and one supplement) are external scientific sources preserved without
changes for comparison with the original publications. They are final
journal versions or ADS scans of the published pages. Publication details
were checked through NASA ADS during the review. The
[manifest](manifest.json) records the original local paths, retrieval URLs
when recorded, page counts, and SHA-256 checksums.

## Primary evidence

| Published source and local PDF | DOI | Passages used and what they establish |
|---|---|---|
| [Zygelman, Dalgarno, Kimura & Lane (1989), Physical Review A 40, 2340](Zygelman_1989_PhRvA_40_2340.pdf) | [10.1103/PhysRevA.40.2340](https://doi.org/10.1103/PhysRevA.40.2340) | Sections II.A and II.B separate the collision calculations. Table I contains radiative rates. Figure 8 compares radiative and nonradiative cross sections. Page 2345 discusses the neglected triplet contribution. The inspected text does not settle the spin normalization of the nonradiative points. |
| [Kingdon & Ferland (1996), ApJS 106, 205](Kingdon_Ferland_1996_ApJS_106_205.pdf) | [10.1086/192335](https://doi.org/10.1086/192335) | Section 3.1, equation (7), and Table 1 give the fit convention, the He+ row, the 6000-100000 K interval, and the reference to Zygelman. No explicit He+ normalization correction was found in these passages. |
| [Stancil & Zygelman (1996), ApJ 472, 102](Stancil_Zygelman_1996_ApJ_472_102.pdf) | [10.1086/178044](https://doi.org/10.1086/178044) | Page 107, footnote 1, explicitly corrects the He radiative transfer and association rates by 1/4. The paper's main collision system is Li + H+; it is included for this footnote. The footnote was checked visually in the published scan. |
| [Stancil, Lepp & Dalgarno (1998), ApJ 509, 1](Stancil_Lepp_Dalgarno_1998_ApJ_509_1.pdf) | [10.1086/306473](https://doi.org/10.1086/306473) | Table 1, reaction (19), gives the radiative He+ + H rate. Reference (10) in the table notes explicitly applies the approach probability correction. The title is *The Deuterium Chemistry of the Early Universe*. |
| [Glover & Jappsen (2007), ApJ 666, 1](Glover_Jappsen_2007_ApJ_666_1.pdf) | [10.1086/519445](https://doi.org/10.1086/519445) | Table 1, R26, supplies the `1.25e-15 * (T/300)^0.25` expression MoCHII used until 2026-09-26, when it moved to SLD98 reaction (19), `1.20e-15 * (T/300)^0.25`. R27 is the older fit for the collisional reverse. |
| [Courtney et al. (2021), ApJ 919, 70](Courtney_2021_ApJ_919_70.pdf) | [10.3847/1538-4357/ac0f57](https://doi.org/10.3847/1538-4357/ac0f57) | Section 2.1, equations (4)-(6), define the approach probability. Section 2.2 describes radiative charge transfer. The radiative calculations do not adjudicate the KF96 nonradiative normalization. |

## Independent comparisons and supporting context

| Published source and local PDF | DOI | Passages used and limits |
|---|---|---|
| [Kimura, Lane, Dalgarno & Dixson (1993), ApJ 405, 801](Kimura_1993_ApJ_405_801.pdf) | [10.1086/172410](https://doi.org/10.1086/172410) | The reverse reaction H+ + He -> H + He+; p. 802, Table 1 and Table 3. The tables were inspected visually. Neither the tabulation nor its later fits are accepted here as a replacement without threshold and consistency checks. |
| [Loreau, Ryabchenko & Vaeck (2014), Journal of Physics B 47, 135204](Loreau_2014_JPhB_47_135204.pdf) | [10.1088/0953-4075/47/13/135204](https://doi.org/10.1088/0953-4075/47/13/135204) | Sections 2 and 3.1 describe independent quantum calculations and their disagreement with Kimura (1993). Section 3.4 and Table 2 give broad-range fits. Process (1) is capture into H(1s); process (2) sums over final H states. The paper uses `eV/u`; conversion and threshold behavior must be checked before thermal integration. |
| [Loreau et al. (2018), Journal of Physics B 51, 085205](Loreau_2018_JPhB_51_085205.pdf) | [10.1088/1361-6455/aab425](https://doi.org/10.1088/1361-6455/aab425) | Methods, Table 1, and Section 3 describe singlet and triplet manifolds and excited n = 2, 3 reactants. These excited-state rates are not direct substitutes for the ground-state rate under review. |
| [Del Zanna et al. (2020), ApJ 898, 72](Del_Zanna_2020_ApJ_898_72.pdf) | [10.3847/1538-4357/ab9d84](https://doi.org/10.3847/1538-4357/ab9d84) | Section 2.7 describes charge exchange, including Zygelman ground-state data and a high-energy extension using Loreau (2014). It also discusses He++ and excited-reactant processes, which must be distinguished from He+(1s) + H(1s). |

## Supplementary data without the article

| Source and local file | DOI | What it contains and its limits |
|---|---|---|
| [Igarashi (2023), Physica Scripta 98, 055402, supplementary data](Igarashi_2023_psacc823supp1.pdf) | [10.1088/1402-4896/acc823](https://doi.org/10.1088/1402-4896/acc823) | Tables behind the article's figures 1-5 in He+ + H at 1-100 keV/u (semiclassical atomic-orbital close coupling, two-electron treatment, per the ADS abstract): elastic, total electron capture, excitation and ionization cross sections, and H(2l), H(3l), Balmer-alpha excitation. The article was not obtained, so the energy convention (taken as lab energy per nucleon of the He+ projectile) and the final states summed in "capture" are not confirmed from the text. See `../IGARASHI_2023_SUPPLEMENT_2026-09-25.md`. |

## Reading and provenance notes

- Extract searchable text with `pdftotext -layout` when checking the papers.
  Kimura (1993) and Stancil & Zygelman (1996) have scan-only body text;
  extraction alone returns essentially page identifiers. Inspect the PDF
  pages or use OCR, and verify any quotation visually against the scan.
- PDF viewer page numbers can differ from printed page numbers. The
  Stancil & Zygelman footnote is on printed p. 107, PDF page 6. The Kimura
  tables are on printed p. 802, PDF page 2. The Loreau (2018) file begins
  with an IOP cover page.
- The local reference directory supplied six PDFs. The preceding research
  session supplied the ADS scans for Kimura (1993) and Stancil & Zygelman
  (1996), and the publisher PDF for Loreau (2014) from the ULB repository.
  An existing research copy supplied the published Stancil et al. (1998)
  deuterium paper.
- Temporary research files with a `.pdf` suffix were checked for actual
  PDF content before inclusion. HTML responses and the different Stancil
  (1998) carbon-transfer paper, ApJ 502, 1006, were excluded. The filenames
  used here identify the journal, volume, and page or article number.
- The ADS search results contain unrelated hits. Only the relevant
  bibliographic records are retained in the manifest. API credentials are
  not included.

The read scope above records the passages used in this assessment. It does
not assert that all results in these publications have been independently
reproduced.
