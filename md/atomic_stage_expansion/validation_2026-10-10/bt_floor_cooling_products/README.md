# Cold-source collision cooling products

These staged products use original CHIANTI scaled collision nodes and the
adopted constant effective collision strength below the identified source
minimum. Both BT scaling and descaling use the source minimum in that region.
The Maxwellian rate coefficient and detailed balance use the actual gas
temperature. The generation controller changed no production file. Following main-session review, the full atoms and annotated coefficients were installed in `data/atomic`; installation records are in `installation/installation.json`.

| Input or model condition | C IV | Si V |
| --- | --- | --- |
| Bound collision-connected levels | 24 | 169 |
| Source collision temperature interval | 3,200--32,000,000 K | 5,000--50,000,000 K |
| Fit and dense check | 801 temperatures, 1,000--100,000 K | 801 temperatures, 1,000--100,000 K |
| Population model | Stationary electron collisions and spontaneous radiation as electron density tends to zero | Same |
| Excluded processes | Photoexcitation, excited-state ionization, finite-density suppression | Same |
| Maximum rounded fit relative error | 0.344146% | 0.929616% |
| Direct source underflow points | 0 | 94; logarithmic source calculation used |
| Native coefficient-reader relative error, resolved points | 4.44e-16 | 4.44e-16 |
| Registered grid / interior / application relative errors | 1.22e-14 / 3.29e-5 / 8.57e-8 | 1.41e-14 / 3.54e-5 / 2.63e-5 |

The rounded fits and source curves appear in `CIV/CIV.pdf` and `SiV/SiV.pdf`.
Both complete profiles were visually inspected. The source and fitted curves
overlap at plot scale, including Si V's representable cold tail. Registered
cooling used the copied, hash-identical source atoms in each `atom` directory,
with a 401-temperature by 161-density table and 40 grid/interior probes.
Those probes extend to 100 K and 1,000,000 K to check the software table; they
do not establish source accuracy outside the fitted interval.

The `candidate` coefficient files add the source temperature interval, cold
approximation, and raw source hashes at the point of use. Their numeric rows
are identical to the original staged fit files; the actual compiled reader
opened each annotated copy and produced byte-identical 801-point results.
The candidate copies are preserved; their exact bytes are now installed as production data. Main-session checks opened both installed atoms and verified actual native coefficient reading. Numerical products were not regenerated after these location changes.

The constant-strength continuation below the source minimum is an adopted
physical approximation with unverified accuracy there. The reported fit error
measures numerical representation of that adopted source, not its physical
accuracy. The original scaled nodes and source files were not modified.
`controller/state.json` records complete commands, frozen source hashes,
finite resource limits, and step outcomes. `build_result.json` records the
isolated current-source native build. Mg III remains excluded because its
identified 27-state atom omits higher source states whose radiative A values
require separate physical review.
