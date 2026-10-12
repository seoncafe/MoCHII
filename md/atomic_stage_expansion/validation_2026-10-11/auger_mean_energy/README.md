# Auger mean heating and the published KM energy correction

## Physical judgment

The current metal heating consumer fully thermalizes metal photoelectron
energy. For this linear calculation, a physically consistent unconditional
mean total Auger-electron kinetic energy is sufficient. At fixed photon
energy E and initial vacancy s, the mean electron heat is

`Q_mean(E,s) = E - I_s + <EA_s>`.

This follows from linearity of expectation; it does not assert that every
charge outcome has the same energy. Charge populations separately require
the probabilities of the final charges. Conditional outputs need their
joint distribution. Nonlinear secondary-electron deposition requires an
electron spectrum, since in general `<D(K)> != D(<K>)`. Correlation with final
charge is additionally required if the requested observable is conditional
on that charge. The earlier blanket demand for a joint charge/energy cascade
before *any* thermal calculation was too restrictive and has been corrected.

Actual current consumers are `species_rate_path_score`,
`species_gamma_compute`, and `metal_heating` in `src/species_mod.f90`.
`src/ion_score_mod.f90` scores secondary-electron energies only for H I,
He I, and He II. The proposed mean formula is not installed in those paths.
Existing inner-shell charge and heat defects remain unresolved.

## Final journal source and direct arithmetic

Kaastra and Mewe (1993), *Astronomy and Astrophysics Supplement Series* 97,
443-482, bibcode `1993A&AS...97..443K`, is available as the final journal scan
at `~/MoCHII/references/Kaastra_Mewe_1993_AAS_97_443.pdf`. Its SHA-256 is
recorded in [source_means.json](source_means.json). `pdftotext -layout` cannot
extract the raster methods; journal pages 445, 446, and 447 were rendered
and visually read directly in this review. Section 2.2 rejects Auger
transitions with nonpositive emitted-electron energy. Section 3 identifies
the mean-energy discrepancy, omitted subsequent excitation/decay processes,
and subshell correlations. Equation (7) proposes the approximate correction

`EA' = EA - epsilon I`,

using the source's own vacancy energy I and discrepancy epsilon. It is a
correction to an unconditional mean, not a PH1 threshold conversion, a
published electron spectrum, or a prescription for arbitrary final-state
excitation. The published probabilities omit entries below 0.00005, and
printed precision is not atomic accuracy.

The new read-only script `tools/fitting/audit_auger_mean_energy.py` evaluated
all 14 KM subshell rows mapping to the nine physical PH1 channels of the
historical 100 eV audit. The current target remains 77 metal stages; these
counts are not a new exhaustive 77-stage reachability measurement.

| Initial ion | KM vacancy | I (eV) | EA (eV) | epsilon | EA' (eV) |
| --- | --- | ---: | ---: | ---: | ---: |
| Mg I | L1 | 92.0 | 27.7 | 0.087 | 19.696 |
| Mg I | L2 | 54.3 | 35.9 | 0.080 | 31.556 |
| Mg I | L3 | 54.0 | 35.6 | 0.081 | 31.226 |
| Mg II | L1 | 104.0 | 20.1 | 0.108 | 8.868 |
| Fe I | M2 | 60.9 | 43.0 | 0.101 | 36.8491 |
| Fe I | M3 | 59.0 | 41.1 | 0.104 | 34.964 |
| Fe II | M2 | 82.9 | 48.3 | 0.148 | 36.0308 |
| Fe II | M3 | 81.0 | 46.4 | 0.151 | 34.169 |
| Fe III | M2 | 104.9 | 40.8 | 0.203 | 19.5053 |
| Fe III | M3 | 103.0 | 38.8 | 0.207 | 17.479 |
| Ca I | M1 | 46.0 | 11.7 | 0.206 | 2.224 |
| Ca I | M2 | 28.4 | 15.7 | 0.187 | 10.3892 |
| Ca I | M3 | 28.0 | 15.3 | 0.189 | 10.008 |
| Ca II | M1 | 58.0 | 7.1 | 0.209 | **-5.022** |

Ca II M1's negative corrected mean exceeds the conservative storage-rounding
bound of 0.089475 eV in magnitude. A negative kinetic-energy mean is
physically inadmissible. Do not clip it to zero or interpret the rounding
bound as physical uncertainty. This demonstrates that the paper's suggested
mean correction cannot be adopted unconditionally for this complete set.
The original positive EA is not itself a negative energy; it remains a mean
from a different approximate energy convention.

## Remaining implementation requirements

1. Establish consistent vacancy, final-binding, fluorescence, and excitation
   energies. KM, PH1, and current NIST potentials cannot be mixed unchanged.
2. Establish valid charge yields near discrepant onsets. Ca II's PH1 3s
   threshold is below its two-electron final ground-state threshold; Fe III
   and Mg II have KM vacancy energies above 100 eV although PH1 is active
   below 100 eV. A mean correction cannot repair these charge thresholds.
3. For full thermalization, validate the mean ledger before integrating
   shell heat and charge rates. A state-resolved model is a possible route,
   but is not the only mathematically sufficient representation.
4. Require electron spectra before adding nonlinear metal secondary
   deposition. Do not infer this new physics from H/He's existing estimator.

Validation was limited to published definitions, source parsing, Eq. (7)
arithmetic, the sign and rounding checks, and source/artifact hashes.
No rate tables, production source, reference products, spatial profiles,
or plasma outputs were regenerated. No build or plasma regression is
appropriate for this read-only source audit.
