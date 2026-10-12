# O IV full-bound cooling installation

The 147-level electron-excitation/spontaneous-radiation atom and 20-term low-density cooling fit are installed in `data/atomic`. The [manifest](manifest.json) distinguishes installation from stage activation. Every numerical record equals the full-bound candidate validated at 401 temperatures and four densities; the current generator also reproduced all records exactly. The actual installed native reader, populations, cooling, and emissivity consumers passed at 10,000 K and four densities, exactly matching saved native profiles. The actual installed GNU coefficient reader passed all 801 temperatures.

The adopted model imports CHIANTI spontaneous probabilities as atomic-model estimates. It calculates photons from identified physical level-energy differences and conserves collisional input and radiative power. `gf` and the wavelength column are not input quantities to this cooling model. A discrepancy between those unused fields and A is not sufficient evidence that the imported A is incorrect. Direct correspondence with original indexed target states and source A values was checked; some deliberate CHIANTI rate replacements remain incompletely traced. No rate was rescaled from an unidentified stored wavelength. This is acceptance of an explicitly qualified source-data cooling model, not certification of atomic accuracy or of every underlying matrix element. The residual source uncertainty remains in the original audit.

The imported A rates are estimates, and the observed/theoretical level selection does not recalculate scattering resonances. Numerical use is 1,000--100,000 K. The collision calculation's published interval begins at 3,200 K; the user's approved constant-Upsilon approximation applies below it, while actual gas temperature enters Maxwellian rates and detailed balance. Photoexcitation, autoionization, radiation trapping, and excited-state ionization are outside this model. The complete bound atom is required for the density ratio: the previous 35-level selection failed the adopted full-source cutoff comparison.

| Verification | Measured result |
| --- | ---: |
| Full-source/native cooling relative difference | 1.39888101103e-14 |
| Gross-power-normalized energy residual | 4.37526002907e-14 |
| Maximum fit error, 801 temperatures | 0.287515276% |
| Registered application error | 0.277913875% |
| Installed coefficient-reader relative difference | 7.95230548079e-12 |

The source/native input conditions and all seven full-profile plot pages are retained in the parent validation directory. No coupled spatial plasma calculation was run. The production ELEMENT roster remains 43 metal stages and 32 adjacent transitions. O IV has prepared installed data, but production oxygen still has three stages. Existing lower-stage numerical products were preserved; the coefficient listing was recollected from existing fits, now 38 species including H I and 260 terms. This is progress toward the approved 64-metal-stage roster, not completion or activation of that roster.
