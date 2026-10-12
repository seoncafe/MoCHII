# Existing compact atoms and high-density cooling

## Physical judgment

A low-density fit residual cannot be assumed to emit without collisional
de-excitation. For a converged excitation atom, numerical fit or table errors
are not evidence of omitted unsuppressed emission. The existing partial-atom
combination creates such emission when the fitted low-density coefficient is
slightly above its registered reference. A separate problem is that the
installed five-level C I, C II, and N I atoms omit much of the bound-state
excitation represented by those fits. Changing their coverage flag without
checking the source atom would retain the wrong density dependence.

## Matched source comparison

The following values were read from executed 801-temperature comparisons
from 1,000 to 100,000 K. The maximum across densities uses 1e-6, 100, 1e6,
and 1e10 cm^-3. The full source and compact source use identical selected
physical energies, original collision interpolation, Maxwellian electrons,
and spontaneous optically thin radiation. The compact source retains only
its first five states. The installed atom also has fitted collision strengths;
that representation difference was measured separately. No photon field,
spatial domain, gravity, flow, trapping, or ionization balance was calculated.

| Ion | Installed levels | Full bound levels | Maximum five-level cooling deficit at ne=1e-6 (%) | Maximum deficit across densities (%) | Low-density fit/full-source difference (%) |
| --- | ---: | ---: | ---: | ---: | ---: |
| CI | 5 | 42 | 92.678236 | 99.996769 | 0.948500 |
| CII | 5 | 53 | 93.174537 | 97.820077 | 0.295932 |
| NI | 5 | 50 | 83.925896 | 99.997399 | 0.589366 |

The main session inspected all three complete cooling profiles. The full
atom with density-ratio cooling reproduces the independently computed source
density dependence; its remaining fit difference is the low-density fit
error. The legacy partial combination instead produces large artificial
cold-temperature residual peaks at high density. These are ion-coefficient
comparisons, not measured fractions of a nebular cooling budget.

The source-node writer now supports exactly degenerate electron transitions,
which occur in full C I and C II. No level energies were split artificially.
The BT scaling energy remains separate from the physical level gap.

## Validation and installation sequence

Full C I, C II, and N I native source comparisons and actual cooling-table
registration are prerequisites for installation. Source completeness is
limited to the identified bound excitation model and does not independently
certify scattering accuracy, missing radiative channels, or other colliders.
The source-convergence diagnostic does not establish accuracy outside the
1,000--100,000 K comparison interval.

The main session separately ran the existing hydrogen-impact checks with
private full C I and C II atoms and unchanged H(1s) rate files. Published
low-density rates agreed within 0.302%, ground-term detailed balance within
1.14e-9, and the thermal increment/emissivity difference within 4.44e-16.
This verifies the changed downstream H-impact consumer and preserves the
source identity of its ground-term transitions. O I was an unchanged control.

Measured arrays, source/input hashes, complete plots, controller state, and
native results are preserved in [source comparison](validation_2026-10-10/existing_density_convergence/result.json),
[fit and density comparison](validation_2026-10-10/existing_density_convergence/fit_ratio_result.json),
and [independent H-impact validation](validation_2026-10-10/existing_density_convergence/hcoll_main_review/result.json).
Full C I (42 levels), C II (53), and N I (50) are now installed with
`DENSITY_RATIO 1` after full-source native comparisons and current-module
registration passed. Maximum checked application differences were 0.151758%,
0.055870%, and 0.084574%, respectively. The installed file bytes equal the
accepted candidate hashes. Previous atoms and the installation record are
preserved in [installation](validation_2026-10-10/existing_density_convergence/installation/result.json).
Low-density fit files and hydrogen-impact rate files were preserved.
No expanded ionization stage was activated. The earlier nine-ion gate failure
remains a historical result; that entire gate has not been rerun after these
three corrections. Other failing ions require separate source review.

## Additional installed corrections

N II (58 levels, 386 radiative rows, 327 electron pairs) and S II
(70 levels, 619 radiative rows, 2,415 electron pairs) also passed full-source
native population/cooling and current-module registration checks and are
installed. Their maximum checked application differences are 0.138828% and
0.136925%, respectively. The main session inspected both complete four-density
profiles. Their maximum low-density fit/source differences are 0.6625% and
0.3793%; these are fit representation errors, not scattering accuracy.
Original five-level atoms are preserved in the installation record.

O II and O III have subsequently been replaced with the complete 196- and
177-level source atoms. The initial audit identified 31 and 48 reversed
radiative rows under mixed observed/theoretical energies. The correction
preserves E1 line strength with physical frequency cubed and the emitting
statistical weight; the single O II ground fine-structure M1 row retains its
observed-gap source A. Unidentified same-parity rates remain explicitly
approximate effective rates. O III optical A totals use final Storey and
Zeippen (2000) values. Independent native populations, collision/radiation
energy closure, and full registered tables passed. Maximum registered
application errors are 0.197768% and 0.143563%; the low-density fits were
preserved because the coronal sum was directly verified unchanged. These
checks establish numerical representation, not scattering or wavefunction
accuracy. Previous installed atoms and failed diagnostics remain archived.
See [oxygen source correction and full validation](oxygen_valence_convergence.md).
