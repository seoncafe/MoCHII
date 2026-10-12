# Preliminary terminal electron-impact ionization assessment

This is a direct evaluation of local source fits, not a compiled terminal-rate
implementation or a verified thermal closure. The approved 64-stage roster
was chosen for photons below 100 eV. That condition does not make the next
collisional ionization rate identically zero at arbitrary temperature.

The diagnostic compares raw Voronov C_i(T) with next-stage ground-parent
RR+DR alpha_(i+1)(T), at 801 temperatures from 1e3 to 1e5 K. It assumes no
photons above the terminal edge, no neutral-hydrogen charge transfer, and only
adjacent electron processes. The ratio is the two-stage stationary population
ratio in that restricted case, not a predicted nebular fraction. It uses raw
Voronov rather than the runtime HYBRID correction; source-fit extrapolation
accuracy and next-stage DR final-paper provenance are not yet accepted.

| Element | Terminal stage | Physical terminal edge (eV) | Raw C_i/(RR+DR) at 1e5 K |
|---|---:|---:|---:|
| C | 5 | 392.090518 | 2.630576816e-19 |
| N | 6 | 552.067329 | 5.249053251e-28 |
| O | 5 | 113.898952 | 3.636577331e-05 |
| NE | 5 | 126.246910 | 2.003955627e-05 |
| S | 7 | 280.954393 | 2.262123708e-13 |
| AR | 7 | 124.411944 | 4.388555312e-06 |
| MG | 4 | 109.265414 | 1.296216954e-03 |
| FE | 7 | 124.967145 | 2.485619215e-05 |
| SI | 5 | 166.767426 | 1.212742461e-06 |
| CL | 7 | 114.201349 | 1.755617540e-04 |
| CA | 6 | 108.783736 | 6.675500845e-05 |

For Fe, absence of local Badnell M=1 records was not treated as absence of
a recombination model. The actual CHIANTI fe_8 RR/DR APIs and files were
inspected. The RR file provides a six-parameter Badnell form, and the DR file
provides nine terms attributed to Schmidt et al. 2008, A&A 492, 265. These
were directly evaluated for the diagnostic; their final published source is
not yet verified. No missing coefficient was replaced by zero.

Mg IV is the most sensitive terminal boundary in this restricted assessment:
the ratio reaches 0.0012962169536965206 at 1e5 K. This does not justify
omission merely because the effect is small. It establishes that a 100 eV
photon criterion alone cannot certify a thermally unrestricted ion roster.
The implementation must either include physically required higher stages or
state and enforce the accepted truncation approximation in its temperature
and ionization domain. No such new runtime policy has been implemented here.

Required continuation: verify next-stage source validity, apply the actual
HYBRID factors and permitted charge-transfer terms, compare complete
ion-fraction profiles under matched inputs, and assess collisional energy
loss and nonadjacent Auger channels. This source diagnostic is not marked
complete for production acceptance.

Measured source output:
[JSON](validation_2026-10-10/terminal_ci_source_assessment.json).
