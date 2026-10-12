# O V full-bound cooling installation

The 166-level source-data excitation atom and 19-term cooling fit are installed with every validated numerical record preserved. The current source generator reproduced all atom records exactly. The actual installed native reader, populations, cooling, and emissivity consumers at 10,000 K and four densities exactly matched saved native profiles. The actual installed GNU coefficient reader passed 801 temperatures, with maximum resolved relative difference 6.12732087291e-12. See the [manifest](manifest.json) and [installed consumer result](installed_native_result.json).

The adopted model uses imported CHIANTI spontaneous rates as atomic-model estimates and computes photons from identified physical level-energy differences. It uses neither source gf nor the wavelength column in cooling. Their inconsistency alone does not establish a wrong imported A. This adoption is based on the actual input quantities and conserved rates/power, not on treating numerical agreement as proof of atomic accuracy. Detailed source-rate accuracy and the observed/theoretical scattering-energy selection remain unverified. No arbitrary wavelength-based rescaling was made.

The distributed O V collision grid starts at 5,000 K, unlike the general temperature interval printed in the final 2014 paper. The documented constant-Upsilon approximation applies below the actual table boundary; actual gas temperature enters Maxwellian factors and detailed balance. This distributed boundary is not a published accuracy limit. Numerical use is 1,000--100,000 K. Scattering resonances are not recalculated after the identified energy selection. Photoexcitation, autoionization, radiation trapping, and excited-state ionization are excluded.

| Verification | Measured result |
| --- | ---: |
| Full-source/native cooling relative difference | 2.02060590482e-14 |
| Gross-power-normalized energy residual | 6.63342718617e-14 |
| Maximum 19-term fit error, 801 temperatures | 0.139092942% |
| Registered application error | 0.099054140% |

All seven full-profile PDF pages and the original conditional results are preserved. The model's qualified adoption is a separate record; absolute source accuracy remains unverified. The combined coefficient listing now contains 39 species including H I and 279 terms. No existing lower-stage fits were recalculated. The production oxygen ELEMENT record still has three stages, and the full metal roster remains 43 stages/32 transitions. O IV/O V installed data do not activate additional stages. No coupled spatial plasma calculation was run. The remaining 19 added-ion data and thermal/charge closure requirements are separate.
