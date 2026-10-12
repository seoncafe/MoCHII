# Zhang et al. (2026) final-paper verification for Fe recombination

## Source identity and scope

The supplied IOPP PDF identifies itself as C. Y. Zhang et al., *The Astrophysical Journal Supplement Series* **285**, 50 (2026 August), DOI [10.3847/1538-4365/ae80c0](https://doi.org/10.3847/1538-4365/ae80c0). Its metadata lists IOPP as creator. It has ten pages, was published on 2026 July 23, and has SHA-256 `ba91543363d26b11593c19ade6c731484a7ed8b39da7679493399c35cfa40d50`. The examined file is `/nfs/mocafe/kiseon/MoCHII/references/Zhang_2026_ApJS_285_50.pdf`. I read its methods, results, and conclusions using `pdftotext -layout` and visually inspected page 9, including Tables 6 and 7 and Equations (4)–(6). The paper has no appendix. Its author data archive, `/nfs/mocafe/kiseon/MoCHII/references/Zhang_2026_Fe3plus_to_Fe6plus_supplement.zip`, has SHA-256 `ed44a5f76b824ec462c6807e3566d4bb82ebb0e4064e5856e398882a5f9d0033` and is identified in the paper by DOI [10.5281/zenodo.20264991](https://doi.org/10.5281/zenodo.20264991).

The authors calculated zero-density dielectronic recombination (DR) and radiative recombination (RR) separately for capture by Fe4+, Fe5+, and Fe6+ as well as Fe3+. They used multiconfiguration Breit–Pauli AUTOSTRUCTURE in the independent-processes distorted-wave approximation. Their target configurations include both `n=3` and `n=4`; low-lying resonance energies received semiempirical corrections. They calculated DR from the `n=3,4` configurations separately from higher Rydberg states, then combined them. The compact fits assume the recombining ion occupies its initial ground level. The authors recommend structure calculation `S1`; `S2` is a sensitivity comparison using different orbital optimization. The paper reports no direct experimental benchmark for these lower charge states. It describes its low-temperature uncertainty through `S1`/`S2` differences and near-threshold resonance sensitivity, rather than a measured experimental error bound.

## Exact published equations and coefficients

The final-paper Table 6 exactly matches `DR_FITS` in `tools/fitting/zhang_fe_recombination.py` for the three requested ions. Equation (4) is `alpha_DR(T) = T^(-3/2) sum_i c_i exp(-E_i/T)` with `T` and `E_i` in K and `c_i` in `cm^3 s^-1 K^(3/2)`. Table 6 gives the following ordered `(c_i, E_i)` terms for the initial ground level of the recombining ion:

| Capture by | Published `S1` DR terms `(c_i, E_i)` |
| --- | --- |
| Fe4+ | `(2.012e-6,188.5)`, `(8.398e-6,1624)`, `(4.158e-5,9845)`, `(8.589e-5,45530)`, `(1.216e-3,262200)`, `(3.703e-2,743800)` |
| Fe5+ | `(7.126e-6,501)`, `(5.411e-5,3001)`, `(1.810e-4,16270)`, `(2.245e-4,62810)`, `(4.365e-3,379200)`, `(3.994e-2,784700)`, `(5.979e-5,14360000)` |
| Fe6+ | `(2.325e-5,1017)`, `(9.656e-5,9864)`, `(4.079e-4,47770)`, `(7.939e-3,279400)`, `(1.336e-1,776800)`, `(4.337e-3,1638000)` |

Equation (5) gives `alpha_RR(T) = A / {sqrt(T/T0) [1+sqrt(T/T0)]^(1-B_T) [1+sqrt(T/T1)]^(1+B_T)}`, with Equation (6) `B_T = B + C exp(-T2/T)`. Table 7 exactly matches the RR parameters transcribed in `fe_recombination_source_decision.md`:

| Capture by | Published RR parameters `(A, B, T0, T1, C, T2)` |
| --- | --- |
| Fe4+ | `(4.242e-10,0.5864,26.78,3.415e7,0.2466,3.415e5)` |
| Fe5+ | `(2.053e-9,0.6747,4.608,3.438e7,0.1254,3.840e5)` |
| Fe6+ | `(2.598e-9,0.7035,6.016,3.394e7,0.0792,3.937e5)` |

Here `A` has units `cm^3 s^-1`, `T0`, `T1`, and `T2` have units K, and `B` and `C` are dimensionless. Direct evaluation at `10,000 K` gives DR/RR in `cm^3 s^-1` of `2.555384524e-11 / 6.148856084e-12` for Fe4+, `8.284940911e-11 / 1.219453689e-11` for Fe5+, and `6.044559831e-11 / 2.046981496e-11` for Fe6+. These independently reproduce the earlier transcribed table calculations. DR and RR must be added only when a total recombination rate is needed; the DR fit is not a unified DR+RR coefficient.

## Temperature validity and author-grid discrepancy

The final paper states that its fits are accurate to better than 2% over `Z^2 (2 x 10^1–1 x 10^7) K`, explicitly calling `Z` the residual charge of the target ion. For Fe4+, Fe5+, and Fe6+, this gives `320–1.6e8 K`, `500–2.5e8 K`, and `720–3.6e8 K`, respectively. These are fit domains, not empirical bounds on the physical atomic rates. The compact formulas have no stated density correction and are ground-level coefficients. At sufficient density, other levels in the ground term can be populated; the paper directs users to its level-resolved archive and population-weighted rates for that case.

I reran `tools/fitting/zhang_fe_recombination.py` against the author archive and checked the ADF09 headers and parent indices. Each `adf09_Fe{4,5,6}+_S1` file has `NUCCHG=26`, `ITYPE=3`, and parent index 1 at zero excitation energy for the expected `3d^4`, `3d^3`, or `3d^2` ground level. The final `ALFT(1)` section supplies 19 total DR temperature nodes. This is a comparison of the published DR fit with the author-provided DR totals for the same `S1` ground level; it does not test RR.

| Capture by | Largest absolute fit-minus-ADF09 difference inside published domain | Temperature |
| --- | ---: | ---: |
| Fe4+ | `+4.712592%` | `320 K` |
| Fe5+ | `+5.119215%` | `500 K` |
| Fe6+ | `+3.068374%` | `1,800 K` |

Thus the publisher-original coefficients confirm the discrepancy with the stated 2% fit accuracy. It cannot be attributed to a copied coefficient, the DR/RR distinction, charge mapping, temperature units, or an `S1`/`S2` substitution. Several grid nodes agree almost exactly while intervening nodes differ more, which is consistent with fitting to a sparse subset of temperatures; the paper and archive do not document the fitting grid, so this is an inference, not a verified cause. At `360 K`, below the Fe6+ fit domain, the fit differs from ADF09 by `-47.700%`; no accuracy is claimed there. The archive's DR totals are printed to three significant digits, far too fine a rounding effect to explain the 3%–5% discrepancies. The final paper and archive do not establish whether the archived totals were revised independently of the fit or whether the stated accuracy was evaluated on another grid.

## Physical acceptance decision

The `S1` ground-level Table 6 DR and Table 7 RR fits are suitable published, ion-specific sources for Fe V to Fe IV, Fe VI to Fe V, and Fe VII to Fe VI within their stated temperature domains. Their adoption should carry an explicit ground-level, zero-density assumption and should not claim verified 2% agreement with the supplied author ADF09 totals. The `S1`/`S2` differences at photoionized-plasma temperatures, quoted in the final paper as up to 10% for Fe4+, 26% for Fe5+, and 30% for Fe6+, indicate model sensitivity; they are not experimental uncertainty limits. The absence of direct benchmarks and the authors' caution about very low temperatures remain material limitations. The journal coefficient tables are now directly verified, so publisher-table access no longer blocks coefficient import. The fit/archive discrepancy is a documented limitation rather than an unresolved coefficient identity; using the ADF09 totals instead would require an explicit choice of interpolation and separate validation. This audit changed no production code or atomic input files, ran no plasma model, and did not validate density-dependent or metastable-level recombination.
