# Mg I restricted cold excitation model

The actual MOCASSIN v2 `mgi.dat` supplies five levels, ten spontaneous
probabilities, and ten effective collision strengths at a single temperature,
10,000 K. `emission_mod.equilibrium` explicitly holds this single strength
constant while using the actual temperature in Maxwellian excitation and
de-excitation rates. Its high-temperature constancy is an approximation,
not additional collision data.

The candidate [native atom](products/nlevel_mg_1.txt) retains those levels,
probabilities, and strengths. Following the user's adopted cold rule, it
prepares only 1,000--10,000 K. A synthetic 1,000 K representation node carries
the same strengths; it is explicitly identified in the file and does not
claim a second measured source temperature. Evaluation above 10,000 K is
rejected. The existing native cold rule extends below 1,000 K, but this
calculation does not validate that colder interval.

| Input condition | Native calculation | Independent matrix calculation |
| --- | --- | --- |
| Levels, energies, statistical weights | All five source states | Identical states |
| Radiation | All ten positive source rates | Identical rates |
| Collision strengths | Constant supplied 10,000 K strengths | Identical constant strengths |
| Maxwellian temperature | Actual temperature | Same actual temperature |
| Temperature grid | 101 points, 1,000--10,000 K | Identical points |
| Electron density (cm^-3) | 1e-6, 100, 1e6, 1e10 | Identical values |
| Geometry, spectrum, ionization equilibrium | Not calculated | Not calculated |

The compiled native reader, level-population solution, emission, and cooling
passed 404 samples. Maximum relative cooling difference was
2.1860291354869332e-13; maximum absolute population difference was
2.220446049250313e-16. An explicit 10,000.1 K evaluation was rejected with the
source-table interval error. This tests numerical reproduction of the adopted
approximation, not atomic accuracy or a full nebular thermal calculation.

All temperature and population curves in [the figure](products/cold_excitation_profiles.pdf)
were visually inspected. [Root review](root_review.json) verifies source and
artifact hashes. The finite controller completed; no calculation remains live.

The atom is not installed. It adds one restricted candidate to the 54 earlier
full supplied-table candidates; it does not complete the 1,000--100,000 K Mg I
requirement. Ca VII remains unresolved, and published Mg I high-temperature
physics remains an independent requirement. No collision temperature was
invented to support high-temperature cooling.
