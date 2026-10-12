# O IV original radiative source identity

This is a source audit, not an adopted O IV cooling atom. The [reproducible
comparison](audit_original_radiation.py) produces [aggregate results](result.json)
from the installed CHIANTI files and an independently downloaded original
OPEN-ADAS ADF04 file. The full ADF04 numerical table remains in `/tmp` for
private scientific inspection because its download terms prohibit
redistribution without permission. No production atomic or fitting data were
changed.

## Published calculation and original data

[Liang, Badnell, and Zhao (2012), A&A 547, A87](https://strathprints.strath.ac.uk/56538/1/Liang_etal_AA_2012_excitation_data_for_the_B_like_iso_electronic_sequence.pdf)
is the final published article, DOI
[`10.1051/0004-6361/201220277`](https://doi.org/10.1051/0004-6361/201220277).
The final PDF is archived locally as
`~/MoCHII/references/Liang_Badnell_Zhao_2012_AA_547_A87.pdf` (SHA-256
`843ee799bce4eefd536f1332bbc92bca3c5cc40f78a1d9a164e6c35672072138`).
Section 2 describes the 204-level AUTOSTRUCTURE target and correction of the
lowest 15 theoretical energies through the Hamiltonian. Section 2.2 explicitly
relates electric-dipole line strength, oscillator strength, transition energy,
and spontaneous rate. The article describes effective collision strengths for
O IV over 3,200--32,000,000 K; this interval does not establish accuracy for
every radiative transition.

The original [OPEN-ADAS O IV ADF04 record](https://open.adas.ac.uk/element/detail/adf04/copaw%5D%5Bb/blike_lgy12%5D%5Bo3.dat)
is `blike_lgy12#o3.dat`, obtained from the official download endpoint. Its
local private copy `/tmp/blike_lgy12_o3_adf04.dat` has SHA-256
`f49d7a8b38d0b82ccff9f4f21871f4e8c785d875521c5be00b017cb9b61407fb`.
The header attributes level energies and radiative rates to AUTOSTRUCTURE and
effective collision strengths to ICFT R-matrix scattering; the original
calculation date is June 18, 2012, and a collision-data revision is dated May
28, 2013. The ADF04 file contains 204 indexed states and 20,706 radiative
pair entries. Each of the 204 statistical weights agrees with CHIANTI's raw
`.elvlc` state of the same original index. The maximum absolute difference
between ADF04 energy and CHIANTI raw theoretical energy is 0.531 cm^-1,
consistent with CHIANTI's energy rounding. This establishes direct original
index correspondence without a configuration-name-only assignment; parent
terms were not silently reassigned.

The final journal source for the separate low-lying rates is Corrégé and
Hibbert (2004), *Atomic Data and Nuclear Data Tables* **86**, 19--34, DOI
[`10.1016/j.adt.2003.11.002`](https://doi.org/10.1016/j.adt.2003.11.002),
ADS bibcode `2004ADNDT..86...19C`. Its publisher final Methods and tabulated
rates were not obtained in this audit. The CHIANTI `.wgfa` footer says rates
for transitions from levels 3--20 were replaced from that source. The exact
scope of the replacement cannot be proved from this sentence alone.

[Polito et al. (2016), A&A 594, A64, Appendix A.2](https://api.repository.cam.ac.uk/server/api/core/bitstreams/8ba92b0b-7cfa-4a0a-b752-02d396c2a6b4/content)
is a final published independent account of the O IV CHIANTI v8 construction.
Its archived final PDF is
`~/MoCHII/references/Polito_et_al_2016_AA_594_A64.pdf` (SHA-256
`bfbd6344989cc94451728e2855533b18e87ddf7e0ccaf87c467b9cb250e0015d`).
Appendix A.2 says CHIANTI selected Corrégé--Hibbert rates for lower levels and
Rynkun et al. (2012) GRASP2K/MCDHF rates for higher levels. Table A.1 lists
six specific O IV transitions among the lowest five levels; the installed
CHIANTI rates agree, after published rounding, with the table's
Corrégé--Hibbert column (the ground fine-structure transition also has the
same value in its Tachiev--Froese Fischer column). The Polito statement is
not evidence that every transition above level 20 uses Rynkun rates.

ADS identifies Rynkun et al. (2012) as *Atomic Data and Nuclear Data Tables*
**98**, 481, DOI [`10.1016/j.adt.2011.08.004`](https://doi.org/10.1016/j.adt.2011.08.004).
Its final published title explicitly limits the calculation to the
`2s2 2p`, `2s 2p2`, and `2p3` configurations; it does not cover all 204
Liang states. An accessible Malmö University PDF was inspected, but its first
page identifies it as a manuscript submitted in July 2011, so it is archived
as `~/MoCHII/references/Rynkun_et_al_2011_submitted_manuscript_not_final.pdf`,
not cited here as final Methods. The publisher final Rynkun article was not
obtained. The accessible final Polito comparison supports low-level rate
selection only.

The official [NASA SolarSoft 2017 CHIANTI archive](https://sohoftp.nascom.nasa.gov/solarsoft/packages/chianti.20170124/dbase/o/o_4/o_4.wgfa)
was downloaded privately to `/tmp/o_4_2017.wgfa` (SHA-256
`6d9c291af2d28f553d659b43f24384bbe2b2e2373066efac70561a20a2c99209`).
Direct comparison, recorded in [source_history_result.json](source_history_result.json),
found the same 4,117 original pair IDs and **zero** `gf`
or `A` differences against the installed file. There are 2,115 wavelength
field differences, reflecting later energy/wavelength processing. Thus the
measured `gf/A` inconsistencies and rate deviations from original ADF04
already existed in the 2017 archived CHIANTI model; the 2019--2023
`wgfa_tidy` passes did not create them.

## Measured comparisons

The current bound selector retains 147 states, 2,319 radiative pairs, and
10,731 electron-excitation pairs. Every selected radiative pair exists in the
original ADF04. Among the selected CHIANTI radiative pairs, 303 spontaneous
rates differ from the ADF04 rate by more than 0.1%; 224 of these have a lower
level numbered 20 or less. Only one pair whose lower level is above 20
differs by more than 1%. Some lower-state differences are orders of magnitude,
so they cannot be described as numerical rounding. These measurements compare
the actual CHIANTI `A` field and original ADF04 `A` field; they do not assert
which value is more accurate. The low-level replacement note does not fully
explain changes with upper levels above 20, so their detailed provenance
remains unresolved.

The [prior frequency diagnostic](../o_iv_main_radiative_review/README.md)
found 343 of 2,064 selected higher-state, positive-`gf`, opposite-parity
transitions whose Einstein frequency from CHIANTI `gf` and `A` differs by more
than 1% from every available raw theoretical gap, observed gap, and stored
wavelength. For **256 of these 343 pairs**, CHIANTI `A` equals the original
ADF04 rate within 0.1%. Thus the selected mixed observed/theoretical energy
vector alone cannot explain these 256 inconsistencies. The original ADF04
does not include the corresponding oscillator strengths needed to determine
which `gf` value entered each original `A` calculation. The final paper's
formulas imply a common transition frequency only if the compared `gf` and
`A` represent the same electric-dipole matrix element and calculation.

Two source-identical examples illustrate distinct issues:

| Original pair | CHIANTI `A` (s^-1) | Original ADF04 `A` agrees? | Frequency from CHIANTI `gf/A` (cm^-1) | Original theoretical gap (cm^-1) | Finding |
| --- | ---: | --- | ---: | ---: | --- |
| 77--87 | 805,000 | Yes | 6,902.1 | 6,928.8 | The stored 209,205 Å wavelength implies approximately 478 cm^-1, whereas the selected energy gap is 6,929 cm^-1 and agrees with the original rate's frequency. |
| 50--118 | 48,300 | Yes | 58,101.5 | 96,537.5 | Even the original theoretical gap does not reconcile CHIANTI `gf` with the source-identical `A`. |

For 77--87, multiplying the published source `A` by the cube of the selected
theoretical-gap-to-stored-wavelength frequency ratio would change the source
rate by a factor of about 3,046. That transformation is not justified by the
original data: the published `A` is consistent with the theoretical gap, while
the stored wavelength reflects a different observed assignment. The
calculation's line strength and level identification have not been verified
against that observed assignment. For 50--118, a frequency-only
adjustment cannot make both stored `gf` and `A` originate from the original
theoretical energies. Neither example licenses a guessed replacement rate.

## Physical status and remaining source checks

The ADF04 supplies a coherent indexed set of calculated states and original
spontaneous rates. It does **not** by itself establish a reliable selected
observed-energy/full-atom cooling model. Some CHIANTI rates differ strongly
from ADF04, and a substantial subset of the retained higher transitions has
an unresolved `gf/A` relation even when the rate itself agrees with the
original source. CHIANTI's lower-state rate changes may be deliberate
corrections, but the available footer is insufficient to trace every changed
pair or establish the final transition multipolarity and line strength.

The NIST-hosted final Tachiev and Froese Fischer (2000), *Journal of Physics B*
**33**, 2419--2435, [publisher-layout PDF](https://nlte.nist.gov/MCHF/all/p256.pdf),
DOI [`10.1088/0953-4075/33/13/304`](https://doi.org/10.1088/0953-4075/33/13/304),
describes MCHF/Breit--Pauli low-state line strengths, oscillator strengths,
and some energy-adjusted rates. The available comparative source literature
places that calculation among the lowest 25 O IV states. It cannot identify
or replace the 147-state atom's many higher radiative rates by numerical
index. The paper was read through the public NIST PDF; direct automated
download returned HTTP 403, so no local final copy was archived. The final
Corrégé--Hibbert publisher DOI resolved to the Elsevier landing page, but
direct ScienceDirect access returned HTTP 403, and the Queen's University
Belfast record had no final PDF. These are access limitations, not evidence
that the papers or data do not exist.

Before changing an O IV radiative rate or adopting a full-atom cooling table,
the remaining physical evidence is the original AUTOSTRUCTURE oscillator
strength or line-strength output for the discrepant high transitions, the
specific CHIANTI conversion and rate-replacement records, and the final
Corrégé--Hibbert final rate table and Methods for low states. The ADF04 and CHIANTI
arrays must remain directly comparable to their original sources. Numerical
population agreement against either selected array would validate only the
implementation, not the missing atomic-source identity.

## Subsequent cooling-path adoption

The preceding sections preserve the unresolved source-quality questions. Root independently reproduced this audit and checked the writer and invoked native consumers: cooling uses imported positive A values and identified physical level energies, without gf or the stored wavelength. The inconsistent unused fields do not themselves establish an incorrect A. Root adopted the complete bound excitation model with imported spontaneous rates explicitly described as atomic-model estimates and atomic accuracy unverified. This is separate from confirming the origin of every rate or using gf for photoexcitation. No A was changed or guessed from the mismatched wavelength.

The selected/source gap comparison is a measurement of source correspondence, not an acceptance rule based on the size of the difference. Observed valence energies and rounded theoretical higher energies are an explicitly qualified atomic-data selection; the scattering resonances are not recalculated. Photon energies and power are computed consistently from that selected vector. See the [installed-model review](../o_iv_full_completion/installation/README.md) for the model's scope, approximations, conserved energy, and exact installed reader checks. The original conditional records remain unchanged; they are superseded by a separate qualified adoption record rather than relabeled as absolute source verification.
