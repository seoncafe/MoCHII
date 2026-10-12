# N IV atomic energy source decision

## Candidate model and acceptance status

The user-supplied final publisher article has now been read and checked;
see [final-paper verification](niv_final_paper_verification.md). A shared
energy selector implements the controlled model below in both low-density
cooling and n-level generation, with source hash, statistical weights,
complete pair counts, and positive-gap checks. No new production N IV atom
has been activated. Numerical validation and selected-level convergence
remain separate from the paper's physical resonance uncertainties.

## Physical decision

Use the 238 level energies in Cloudy c25.00 Stout `data/stout/n/n_4/n_4.nrg`
for the CHIANTI N IV 136 level cooling model. Retain the CHIANTI `.wgfa`
Einstein A values and `.scups` effective collision strengths with their original
level indices. Determine the emitted photon energy, Boltzmann factor, and
detailed balance from the same selected energy differences. Do not use the
`.wgfa` wavelength as a substitute transition energy, or the `.scups` `de`
field as a physical gap; the latter is the Burgess--Tully scaling parameter.

This is an ion-specific source selection, not a general switch to theoretical
CHIANTI energies. The Stout file preserves more precision for high levels and
therefore resolves two artificial zero gaps introduced by rounding the
CHIANTI theoretical energies to whole cm^-1. It also uses measured energies
for its lowest levels. The Stout header states 20 updated levels, but direct inspection shows
that levels 1--29 use observed energies and the remainder agrees with BSR
theoretical energies. This header defect is reported without editing the
external source. The file derives from the ADF04 data provided by
Luis Fernandez-Menchero. The CHIANTI files cite Fernandez-Menchero,
Zatsarinny, and Bartschat (2017), *Journal of Physics B* **50**, 065203,
DOI [10.1088/1361-6455/aa5fc4](https://doi.org/10.1088/1361-6455/aa5fc4),
for their theoretical energies, A values, and effective collision strengths.
The [official CHIANTI N IV record](https://db.chiantidatabase.org/n/n_4.html)
confirms these file roles and says negative `.wgfa` wavelengths flag
approximate theoretical wavelengths.

This selection makes the energies consistent with the theoretical level
identities of the A and collision data; it does not assert that every
theoretical energy is closer to the observed physical level. The
[NIST N IV level table](https://physics.nist.gov/cgi-bin/ASD/energy1.pl?biblio=1&conf_out=1&j_out=1&lande_out=1&level_out=1&perc_out=1&spectrum=N+IV&term_out=1)
independently lists the observed energies of levels labeled `2p 4p 3S1`,
`2p 4p 1S0`, and `2p 4d 1F3` at 595,296.3, 594,182.6, and 591,047.8
cm^-1. These are reported experimental energies, not transcription errors.
Their conflict with the theoretical ordering means that inserting them into
the theoretical rate matrix is not a physically consistent way to improve
its accuracy. The selected model remains an approximation to the observed
spectrum.

## Direct source checks

The local sources inspected were `~/RT_Codes/CHIANTI/dbase/n/n_4/n_4.elvlc`,
`.wgfa`, and `.scups`; CHIANTI `pro/read/read_wgfa2.pro`; and Cloudy c25.00
`data/stout/n/n_4/n_4.nrg`, `.tp`, and `.coll`. CHIANTI's reader describes the
first `.wgfa` index as the lower level and the second as the upper level. Thus
reversing an Einstein A solely because an adopted level energy is inverted
would assign the decay to a different state. Cloudy's CHIANTI reader sorts a
selected energy set and remaps indices, but its source contains a note that
the relation between file index order and physical order may depend on the
energy selection (`source/atmdat_chianti.cpp`, near line 1450).

Direct parsing found 238 energy rows, 3,272 positive-A radiative rows, and
9,180 `.scups` collision rows in the CHIANTI files. The radiative and
collision rows refer only to indices 1--136. All 238 CHIANTI and Stout indices
have the same statistical weight. All 3,272 CHIANTI positive-A pairs exist in
Stout `.tp`: 3,251 A values agree at the printed precision, and 21 differ only
by the final decimal rounding. All 9,180 CHIANTI collision pairs exist in
Stout `.coll`. At 3,200 K, descaled CHIANTI effective collision strengths
agree with the Stout values to their displayed precision for sampled pairs
1--2, 1--5, 29--30, 72--73, 99--100, 107--116, 107--135, 115--135, and
134--135. This supports common index identities and common primary collision
data; the Stout `.coll` temperature grid and interpolation have not been
proven numerically identical to the CHIANTI descaling at every temperature.
The NIST table gives the N IV ionization limit as 624,866 cm^-1; Stout level
137 starts at 625,844.3 cm^-1. The 136-level cut therefore separates
the supplied bound-level model from levels above the ionization threshold.

Level labels are not identical in every row. In particular, CHIANTI level 29
is labeled `2p 3p 3S1`, whereas Stout level 29 is `2s 4s 3S1`. The NIST
table lists `2s 4s 3S1` at 487,607.4 cm^-1, matching Stout's energy and
label; CHIANTI's printed configuration is incorrect for that observed level.
NIST also lists `2p 3p 3S1` at 498,045.5 cm^-1. CHIANTI level 37 carries
that energy but is labeled `2s 4s 3S1`, and Stout level 37 at 497,406.9
cm^-1 is likewise labeled `2S1 4S1 3S1`. Thus neither file's printed
configuration labels identify both mixed `3S1` states consistently.
The Stout high-level configuration notation also differs from CHIANTI
notation. The complete index, statistical-weight, radiative-rate, and
collision-pair checks above establish the operational mapping; they do not
prove that every printed high-level spectroscopic label is correct.

With CHIANTI's observed-when-available energies, 14 radiative rows have a
negative upper-minus-lower gap, 81 collision rows have a negative gap, and
17 collision rows have an exact zero gap. Replacing only observed energies for
levels 107, 134, and 135 with Stout energies removes all 14 radiative
inversions, but leaves 30 collision inversions and 17 exact collision ties.
With the Stout energy set, every one of the 3,272 radiative and 9,180
collision pairs has a positive gap. The minimum adjacent Stout gap among all
238 levels is 0.1 cm^-1, between levels 176 and 177.

| CHIANTI lower--upper | Mixed gap (cm^-1) | Stout gap (cm^-1) | `.wgfa` wavelength (Angstrom) | Stout-gap wavelength (Angstrom) |
|---|---:|---:|---:|---:|
| 107--116 | -1,606 | 5,088.4 | 62,266.50 | 19,652.5 |
| 107--117 | -1,583 | 5,117.0 | 63,171.19 | 19,542.7 |
| 107--123 | -807 | 5,816.2 | 123,915.74 | 17,193.4 |
| 107--124 | -749 | 5,864.6 | 133,511.35 | 17,051.5 |
| 107--125 | -723 | 5,887.5 | 138,312.59 | 16,985.1 |
| 124--134 | -365 | 3,245.2 | 273,972.60 | 30,814.7 |
| 115--135 | -13 | 7,098.6 | 7,692,307.69 | 14,087.3 |
| 119--135 | -2,758 | 4,088.9 | -24,455.857 | 24,456.5 |
| 121--135 | -2,840 | 4,007.2 | -24,956.326 | 24,955.1 |
| 122--135 | -2,861 | 3,986.2 | -25,087.807 | 25,086.6 |
| 127--135 | -3,363 | 3,483.6 | -28,702.641 | 28,706.0 |
| 129--135 | -3,577 | 3,269.6 | -30,581.040 | 30,584.8 |
| 131--135 | -5,163 | 1,683.6 | -59,382.420 | 59,396.5 |
| 133--135 | -5,469 | 1,377.5 | -72,568.940 | 72,595.3 |

The first seven positive `.wgfa` wavelengths above reproduce the magnitude
of the *inverted observed* gap rather than a physically directed emitted
photon. The seven negative wavelengths are based on approximate theoretical
energies, as documented by CHIANTI, and are broadly consistent with the
Stout gaps. These 14 records implicate the CHIANTI observed entries of
107 (`2p 4p 3S1`, 595,296 versus Stout 588,134.6 cm^-1), 134
(`2p 4p 1S0`, 594,182 versus 597,244.4 cm^-1), and 135
(`2p 4d 1F3`, 591,047 versus 597,894.4 cm^-1). Other observed high-level
entries also disrupt the collision ordering, which is why changing only
those three is insufficient for the complete rate matrix.

CHIANTI rounds theoretical levels 72 and 73 to the same 553,257 cm^-1 and
levels 99 and 100 to the same 586,262 cm^-1. Stout gives levels 72--73
553,257.0 and 553,257.5 cm^-1, and 99--100 586,262.2 and 586,262.8
cm^-1. Their `2s 5g 1G4`/`2s 5g 3G5` and `2s 7p 3P0`/`2s 7p 3P1`
identifications are distinct. The small splittings are physically resolved
by the ADF04 precision; neither collision row should be discarded or
mistaken for a radiative zero-energy line. Their `.scups` `de=1 Ry` values
control spline scaling, not the thermal threshold. At 1,000 K the two Stout
gaps are less than 0.1% of kT, but they remain positive throughout the model.

## Implementation and acceptance scope

Read the N IV Stout energies as an explicit ion-specific source with a fixed
1--238 index mapping and select indices 1--136 for the present CHIANTI
collision and radiative model. Validate the 238-row count, contiguous indices,
statistical weights against `.elvlc`, positive gaps for every included
`.wgfa` and `.scups` pair, and the demonstrated A/pair mapping. Fail with a
source-specific diagnostic if those checks change. Keep Einstein A attached to
its original transition and use selected energy differences for all photon
energies and forward/reverse collision rates. A `.wgfa` weighted oscillator
strength derived from its old wavelength must not be treated as independent
of that wavelength if downstream code needs it.

The next validation tested the N IV rate-matrix solve and low-density cooling over the intended
temperature range, checking finite nonnegative populations, positive total
cooling, and convergence as electron density tends to zero. Compare line
contributions involving 107, 134, 135 and the two formerly tied collision
pairs before accepting any regenerated N IV coefficients. The initial source audit did not regenerate cooling products or run a full model.
Subsequent scoped staging tests have now passed; see
[N IV staging validation](niv_staging_validation.md). Production activation and
a transported model comparison remain separate acceptance steps.

The earlier ADS DNS and publisher access failures were resolved for this
article by the user-supplied final publisher PDF. The final methods and
Table 1 were read directly. Printed original BSR indices beyond level 87
cannot be equated to the selected 238-state indices. The complete author
supplement was not independently retrieved; therefore exact high-level
label assignments and convergence beyond 87 remain uncertified. The
implemented approximation uses the checked local pair mapping and does
not claim an improved observed spectrum or established low-temperature
resonance accuracy.
