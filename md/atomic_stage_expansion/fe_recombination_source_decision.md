# Fe V--VII recombination source decision

## Physical decision

A current, stage-specific source **does exist** for all three requested transitions. Zhang et al. (2026), *ApJS* **285**, 50, [doi:10.3847/1538-4365/ae80c0](https://doi.org/10.3847/1538-4365/ae80c0), calculated dielectronic and radiative recombination (DR and RR) for capture by Fe4+, Fe5+, and Fe6+ (Fe V, Fe VI, and Fe VII). This reverses the earlier conclusion based only on installed, older tables. The paper presents separate ground-level DR and RR fits, so its DR can replace the inappropriate CHIANTI/Mazzotta-era low-temperature DR and its RR can be used consistently. The authors recommend their `S1` DR calculation. These lower iron charge states lack a direct experimental benchmark, and the lowest temperatures remain sensitive to near-threshold resonance energies.

The final journal PDF has now been supplied and verified. Tables 6/7 and
Eqs. (4)--(6) establish the selected ground-level, zero-density S1 DR and RR
for recombining charges q=3--6. The fit domain is q^2*(20--1e7) K;
analytic extrapolation outside it has no certified accuracy. The compact
DR fits exceed the stated 2% difference at some author ADF09 nodes, which
is an explicit source limitation rather than a transcription error.

The generator selects these fits and permits Fe V with the reviewed Fe4
charge-transfer estimate. Fe VI/VII still fail preflight because their
higher-charge electron-transfer rates remain unselected. (Corrected
2026-10-11, audit Section 5 item 7: stale. The Langevin rows for Fe5+ and
Fe6+ + H exist (`data/atomic/element_fe.txt`, used with
`par%langevin_charge_exchange`, default `.false.`), and Fe VI and Fe VII are
production stages since the 77-stage roster was installed on 2026-10-11.) The existing
Fe IV recombination data were also corrected: the previous high-temperature
DR fit missed low-temperature resonances. See
[final verification](zhang2026_final_paper_verification.md) and
[implementation and measured validation](final_fe_recombination_validation.md).

The remainder records the historical source audit. Its former publisher
access and provisional-fit statements are superseded by the verification
above; old audit products are preserved as historical evidence.

## New published source and coverage

The 2026 paper's title is *Dielectronic Recombination of Iron M-shell Ions: Fe6+--Fe3+*, by C. Y. Zhang, N. R. Badnell, P. J. Storey, G. Del Zanna, M. G. O'Mullane, and X. Ma. The accessible indexed typeset text describes multiconfiguration Breit--Pauli AUTOSTRUCTURE calculations with interactions between `n=3` and `n=4` configurations and semiempirical corrections to near-threshold resonance energies. Its Table 6 gives `S1` ground-level DR fits for Fe6+ through Fe3+; Table 7 gives separate RR fits. This includes an additional Fe IV to Fe III transition outside the present three-stage request. The companion [author data archive](https://doi.org/10.5281/zenodo.20264991) contains `adf09_Fe{4,5,6}+_S1` and associated atomic calculation files, which independently establishes ion coverage and the selected calculation. Its ADF09 files also contain total ground-level DR grids under `ALFT(1)`, in addition to level-resolved data. The [official ADAS `xxdata_09` documentation](https://www.adas.ac.uk/python_lib/xxdata_09.html) identifies these as total DR coefficients in `cm^3 s^-1` at electron temperatures in K. This is not a direct copy of the journal's compact Table 6/7 fits.

The final-publisher fit expressions are, with electron temperature `T` in K, `alpha_DR = T^(-3/2) sum_i c_i exp(-E_i/T)` and `alpha_RR = A/[sqrt(T/T0) (1+sqrt(T/T0))^(1-B_T) (1+sqrt(T/T1))^(1+B_T)]`, where `B_T=B+C exp(-T2/T)`. The Table 6 `c_i` units are `cm^3 s^-1 K^(3/2)` and `E_i` units are K. The authors report that the fits reproduce their calculated rates within 2% over `Z^2 (20--10^7) K`, where `Z` is the residual ionic charge. Thus the reported intervals are `320--1.6e8 K` for Fe4+, `500--2.5e8 K` for Fe5+, and `720--3.6e8 K` for Fe6+. That claimed fit accuracy is **not reproduced** by the archived total grids, as quantified below; it is distinct from the physical uncertainty of the atomic calculation.

The final-paper Table 6/7 values give the following comparison at `10,000 K`. These numbers were computed directly from the coefficients. The subsequent publisher check confirmed the transcription; linked rate validation is recorded above.

The verified `S1` Table 6 terms are ordered pairs `(c_i, E_i)` with units stated above:

| Recombining ion | Ordered `(c_i, E_i)` pairs |
| --- | --- |
| Fe4+ | `(2.012e-6,188.5)`, `(8.398e-6,1624)`, `(4.158e-5,9845)`, `(8.589e-5,45530)`, `(1.216e-3,262200)`, `(3.703e-2,743800)` |
| Fe5+ | `(7.126e-6,501)`, `(5.411e-5,3001)`, `(1.810e-4,16270)`, `(2.245e-4,62810)`, `(4.365e-3,379200)`, `(3.994e-2,784700)`, `(5.979e-5,14360000)` |
| Fe6+ | `(2.325e-5,1017)`, `(9.656e-5,9864)`, `(4.079e-4,47770)`, `(7.939e-3,279400)`, `(1.336e-1,776800)`, `(4.337e-3,1638000)` |

The verified Table 7 RR parameters are ordered `(A, B, T0, T1, C, T2)`:

| Recombining ion | RR parameters |
| --- | --- |
| Fe4+ | `(4.242e-10,0.5864,26.78,3.415e7,0.2466,3.415e5)` |
| Fe5+ | `(2.053e-9,0.6747,4.608,3.438e7,0.1254,3.840e5)` |
| Fe6+ | `(2.598e-9,0.7035,6.016,3.394e7,0.0792,3.937e5)` |

| Capture by | Recombination | Zhang `S1` DR | Zhang RR | Existing CHIANTI DR | Existing CHIANTI RR |
| --- | --- | ---: | ---: | ---: | ---: |
| Fe4+ | Fe V to Fe IV | 2.5554e-11 | 6.1489e-12 | 6.0480e-27 | 7.8000e-12 |
| Fe5+ | Fe VI to Fe V | 8.2849e-11 | 1.2195e-11 | 3.8722e-35 | 1.5100e-11 |
| Fe6+ | Fe VII to Fe VI | 6.0446e-11 | 2.0470e-11 | 1.0797e-30 | 2.6200e-11 |

All columns have units `cm^3 s^-1`. The previous CHIANTI values were remeasured from local version 11.0.2 files with `tools/fitting/fe_recombination_expansion.py`. The Zhang `S1` values were evaluated independently from the table transcription. At `1,000 K`, the transcribed `S1` DR rates are `1.0511e-10`, `2.2165e-10`, and `2.6608e-10 cm^3 s^-1` for Fe4+, Fe5+, and Fe6+, respectively. These are within the quoted fit ranges, while the authors caution that uncertainty at very low temperature may be larger. The paper reports `S1`/`S2` differences at photoionized-plasma temperatures of about 10%, 26%, and 30% for these three ions, respectively; these indicate uncertainty, but are not an experimental error bound.

Direct comparison against the authors' ADF09 `ALFT(1)` total DR grids, using `tools/fitting/zhang_fe_recombination.py`, finds 19 temperature nodes per ion. The grids cover `160--1.6e8 K` for Fe4+, `250--2.5e8 K` for Fe5+, and `360--3.6e8 K` for Fe6+; their first nodes fall below the provisionally transcribed fit-validity range. Among nodes *within* that range, the largest fit-minus-ADF09 relative differences are **+4.713% at 320 K** (Fe4+), **+5.119% at 500 K** (Fe5+), and **+3.068% at 1,800 K** (Fe6+). These exceed the paper's reported 2% fit accuracy. At selected nodes `800 K`, `1,250 K`, and `720 K` respectively, agreement is within 0.03%. The discrepancy is concentrated toward the cool end, rather than a systematic unit or charge-index error. The subsequent final-publisher check confirmed those coefficients and methods. The measured disagreement is therefore not a transcription error. Its origin remains unestablished; the implementation retains the published fit with this limitation explicitly recorded. The [complete comparison JSON](fe_recombination_expansion/zhang2026_adf09_comparison.json) records every node. The author archive is stored at `/nfs/mocafe/kiseon/MoCHII/references/Zhang_2026_Fe3plus_to_Fe6plus_supplement.zip`, has SHA-256 `ed44a5f76b824ec462c6807e3566d4bb82ebb0e4064e5856e398882a5f9d0033`, and is 1,730,357 bytes.

Potential older sources are Nahar, Bautista, and Pradhan (1998),
*Phys. Rev. A* 58, 4593, DOI 10.1103/PhysRevA.58.4593, and Nahar and
Bautista (1999), *ApJS* 120, 327. The indexed 2026 discussion describes
unified RR+DR results and differences from its new target and resonance
model. Their final methods have not been verified in this investigation,
and no numerical coefficients from those papers are adopted. Before any
use, determine whether a coefficient already includes RR; adding a unified
total to a separate RR contribution would double count that contribution.

## Cloudy GuPC and density limits

Cloudy c25 `ion_recomb_Badnell.cpp` labels its `Fe_Gu_c[0]` and `Fe_Gu_E[0]` row `Fe+6`, stores it for recombination into Fe+5, and labels its provenance `GuPC: Gu, M. private communication`. That row is only Fe VII to Fe VI. The expression uses eV for temperature and energies, unlike the 2026 journal fit's Kelvin coefficients. At `10,000 K`, the local Cloudy GuPC expression gives `5.8929e-11 cm^3 s^-1`, numerically close to the verified Zhang Fe6+ fit value. Numerical agreement does **not** demonstrate common provenance or a published validity range for GuPC. No later published derivation of this particular GuPC row was established.

Cloudy fills other missing Fe DR with a configurable charge-dependent mean plus an Arnaud--Raymond high-temperature term; that is an estimate, not a separate Fe V/VI low-temperature calculation. Cloudy may subsequently apply collisional DR suppression using density, charge, and temperature. The Zhang fits are for the recombining ion's initial ground level in the zero-density atomic calculation. The authors' archive supplies rate data from other ground-term levels, useful when those levels become populated. A density suppression factor cannot supply a missing baseline coefficient, and it should be assessed separately from replacing that coefficient. A model that populates metastable levels needs level populations and the corresponding initial-level recombination data; simply applying the ground-level fit to all populations is an additional approximation.

## Source-access and verification boundary

The NASA ADS API was attempted first using the local token without printing it, but DNS resolution failed. The publisher's [DOI landing/PDF route](https://iopscience.iop.org/article/10.3847/1538-4365/ae80c0/pdf) was then tried through the web tool and elevated `curl`: the web tool could not access the site and `curl` received a Radware Bot Manager HTML challenge despite an HTTP 200 response. A candidate `content.cld.iop.org` route returned HTTP 403. The [author's Zenodo archive](https://doi.org/10.5281/zenodo.20264991) was accessible and its `supplement.zip` was obtained directly. Search indexing exposes text and Tables 6/7 from an apparently final, typeset, CC BY 4.0 paper at a [ResearchGate page](https://www.researchgate.net/publication/410735149_Dielectronic_Recombination_of_Iron_M-shell_Ions_Fe_6_-Fe_3), but this is a repost and the publisher original PDF was not retrieved. The exact journal-table transcription, equation typography, full figure axes, and any publisher corrections therefore remain unverified against the publisher's original file. That historical access blocker was resolved by the user-supplied final PDF; the measured checks are linked above. No access restriction was bypassed.

The historical generator rejected Fe models with more than four stages before writing an atomic file. Without that guard, stages `i=4,5,6` would read `fe_5`, `fe_6`, and `fe_7` older CHIANTI RR/DR fits. The local Badnell tables have no corresponding ground-level entries, which is only a coverage observation. No production files or reference products were changed during this source audit. Full plasma-profile comparison remains downstream of the complete Fe stage expansion and its matched physical inputs.
