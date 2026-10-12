# Si I electric dipole source audit

## Result and scope

This audit preserves the official NIST Si I MCHF original transition file and its formatted display table, then matches their state-resolved records against every positive-`gf`, positive-A row of the local CHIANTI `si_1.wgfa`. All 36 CHIANTI state identities are checked by configuration, term, parity, and J. NIST's original calculation has records for states 1–35; CHIANTI state 36, above the bound-state threshold, is absent. No transition probability, source file, native atom, fit, or production code was changed.

The 35-state bound CHIANTI atom has 229 positive-`gf`, positive-A rows. Exactly 185 match a unique original NIST E1 record. The other 44 are not E1-allowed: 43 join states of equal parity and one joins opposite-parity states with `Delta J = 2` (a NIST M2 record, pair 1–6). Ten of these 44 have a matching original NIST E2, M1, or M2 record; the other 34 have no record in this particular NIST calculation, and their multipoles remain unclassified by this audit. All 44 are oriented toward lower observed energy. The single reversed positive-A row among the NIST E1 matches is CHIANTI 15–14. These checks establish **complete NIST E1 identity coverage of the supplied bound radiation rows**; they do not establish complete spontaneous radiation or accuracy of the 34 unmatched same-parity rates.

## Original calculation and duplicated CHIANTI values

The [NIST collection page](https://nlte.nist.gov/MCHF/Elements/Si/Si_14.16.html) describes the official 2002 MCHF *ab initio* calculation. The exact downloaded [original `.lsj` file](Si_14.16.cff-mchf.lsj.db) has SHA-256 `438150a064b12c202433f590122da7eaf819287e4656d0a6c5ba88586744d904`; it contains 235 original state-resolved transitions: 196 E1, 13 E2, 10 M1, and 16 M2. The exact downloaded [NIST formatted transition table](Si_14.16.cff-mchf-lin.dat.mp) has SHA-256 `049b438f9513a045f9cc2c36dba052f88734a2710d2015da88174b6713b592ee`. These are unchanged external scientific data files retained for direct comparison.

For 170 of the 185 matched E1 rows, CHIANTI's `gf` agrees with the one original NIST length-form `GF` to within `6.23e-4` fraction, consistent with printed rounding. For the other **15 rows**, CHIANTI's `gf` is approximately twice the one original `GF`, and CHIANTI A is likewise approximately twice the original A. Each of these 15 original transitions is displayed twice in the NIST formatted table but occurs once in the original `.lsj`; all 170 singly displayed transitions have the near-unity `gf` ratio. This one-to-one correspondence strongly suggests duplicate display rows were counted twice during source preparation. The preparation code was not available, so the mechanism is an inference, while the numerical duplication is directly established. The current CHIANTI `.wgfa` header names multiple calculations and a later formatting step; it does not identify which step caused the duplication.

| CHIANTI indices | One original NIST A (s^-1) | CHIANTI A (s^-1) | CHIANTI `gf` / original `GF` |
|---|---:|---:|---:|
| 4–10 | 2.3744e8 | 4.748e8 | 2.000394 |
| 5–10 | 1.30767e7 | 2.616e7 | 1.999969 |
| 10–14 | 8.66695e6 | 1.733e7 | 1.999887 |
| **15–14** | **66.4201** | **132.8** | **1.999899** |
| 4–15 | 4.06728e7 | 8.134e7 | 2.000644 |
| 10–26 | 2.15864e7 | 4.318e7 | 2.000228 |
| 15–26 | 164167 | 328400 | 2.000790 |
| 10–30 | 2.28858e7 | 4.578e7 | 2.000677 |
| 4–31 | 3.30655e8 | 6.614e8 | 2.000114 |
| 26–31 | 2.1533e6 | 4.306e6 | 1.999836 |
| 4–32 | 3.73573e7 | 7.472e7 | 2.000437 |
| 5–32 | 1.22833e8 | 2.456e8 | 2.000089 |
| 14–32 | 4.42937e6 | 8.858e6 | 1.998881 |
| 26–32 | 171968 | 344000 | 2.000115 |
| 30–32 | 272560 | 545200 | 1.999973 |

The first and second indices in the table follow CHIANTI's row order; physical lower and upper state identities are stored separately in `result.json`. The remaining 170 source `gf` ratios span `0.999531965–1.000622956`. The 15 doubled A ratios span `1.999394762–2.000504714`. Comparing CHIANTI A with original A tests source-row duplication; it is **not** a comparison of probabilities transformed to observed wavelengths. The full pointwise match, source line numbers, E1 strengths in both gauges, energies, emitting weights, and selected observed-energy transformations are recorded in [result.json](result.json).

## Physical energy order and gauge sensitivity

The only NIST E1 ordering inversion against both CHIANTI's observed and its own theoretical level energies is pair 15–14. NIST's calculated state 15 (`3s2 3p 3d 1D`, odd, `J=2`, `g=5`) lies at 47188.06 cm^-1; state 14 (`3s2 3p 4p 1P`, even, `J=1`, `g=3`) lies at 47317.64 cm^-1. The selected observed levels reverse that order: state 15 is 47351.555 cm^-1 and state 14 is 47284.062 cm^-1. Their physical emission is 15 to 14 with a 67.493 cm^-1 photon. The NIST original record has one length-form E1 `S=45.1996` atomic units and A `66.4201 s^-1` at its calculated 129.58 cm^-1 gap and emitting weight 3. The original velocity form has `S=13.3406` atomic units, a factor 3.388 lower. Applying the standard E1 frequency-cubed and emitting-weight relation to the observed physical gap gives length-form `A(15 -> 14)=5.631292575 s^-1`; the conditional velocity-form result is `1.662068287 s^-1`. Neither is certified here as an accurate absolute rate. NIST's collection page cautions that its Breit-Pauli velocity operator lacks some relativistic terms when term mixing is significant; gauge spread is a sensitivity diagnostic, not a confidence interval.

The NIST original has 11 weak E1 records between bound states that CHIANTI does not include. Their state pairs and rates are in `result.json`; the largest original A among them is `3.03704 s^-1` for pair 11–26. These absences do not explain the reversed 15–14 row and should be assessed before claiming a complete line list. The original calculation does not include CHIANTI state 36, so the 13 CHIANTI rows to that above-threshold state have no counterpart here.

The final [Froese Fischer (2005), *Physical Review A* 71, 042506](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.71.042506), DOI `10.1103/PhysRevA.71.042506`, was identified through ADS. Its publisher PDF returned HTTP 401, and its individual 3p3d-to-3p4p transition data and methods were not checked. That final article's accessible abstract discusses good length/velocity agreement for **ground-configuration-to-3p3d** allowed transitions; it cannot validate this near-degenerate excited-to-excited pair or be silently substituted for the NIST 2002 calculation. The 2002 NIST calculation is nevertheless a directly preserved physical source. Adopting its length-form strengths as a defined approximation requires explicit physical judgment and full finite-density source validation; it does not require inventing an A value.

## Reproduction

Run `python md/atomic_stage_expansion/validation_2026-10-10/si_i_electric_dipole_source_audit/audit_si_i_electric_dipole.py` from the project root. The [script](audit_si_i_electric_dipole.py) reads the two preserved NIST files and the unchanged local CHIANTI `.elvlc` and `.wgfa`, verifies all four input SHA-256 values, and writes [result.json](result.json). It uses only Python's standard library. CHIANTI input hashes are `1538e98411f8c2f03414aebd36f901f7b88bc79109b1f0aea7254af070bd4f57` (`.elvlc`) and `be54f54d2e91695b766c0cb545e744d0d540165add1231775b7361f82ab21c42` (`.wgfa`). The script checks that the 15 doubled `gf` pairs are exactly the 15 twice-displayed original pairs and rejects unexpected state or E1 selection mappings. It is a read-only source audit, not an atomic-data generator.

## Private completed radiation candidate

The [candidate module](si_i_spontaneous_probabilities.py) exposes `source_atom()`, returning the 35 selected observed-energy states, corrected radiation, and the 595 original bound collision rows. It checks the hashes of all original inputs, all 35 source state identities, all 235 original transition multipoles against parity and angular-momentum selection, statistical weights, and the NIST printed strength-to-A relation. The original collision abscissae, ordinates, and Burgess–Tully scaling parameters remain unchanged. Its sole collision continuation is a constant effective strength below an adopted 1003.02 K floor. The first positive nodes of the rounded source fits imply 997.095–1003.020 K; this floor is an inferred numerical choice, not a published physical validity or accuracy bound. Maxwellian rates still use the actual gas temperature.

The private [source-node atom](candidate_v2/nlevel_si_1.txt) has 35 levels, 250 unique radiative pairs, 595 collision rows, and 8,925 original collision nodes. Its [result record](candidate_v2/candidate_result.json) contains the atom hash, module hash, original input hashes, full paths, counts, and floor. The 250 radiative pairs consist of the original 229 CHIANTI pairs plus 21 previously absent NIST pairs. The latter combine 26 original records: 11 E1, 6 E2, 5 M1, and 4 M2. Eleven already represented E1 pairs also have a small identified M2 contribution, which is summed into the same pair. Original E1 strengths are converted with the selected physical photon frequency, emitting statistical weight, and the E1 frequency-cubed law. E2 and M2 use the fifth power; M1 uses the third power. Existing non-E1 CHIANTI rates on 44 pairs remain unchanged approximate effective rates, including 10 pairs for which this NIST calculation lists a forbidden multipole. These are not counted twice. Original positive `gf` on those rows is not treated as a physical E1 oscillator strength.

The previous [provisional candidate](candidate/nlevel_si_1.txt) is retained as an audit artifact from before the final state-selection guards. `candidate_v2` is the source-hash-frozen input for private numerical validation. Neither candidate is installed in the project atomic-data tree. The original 2002 NIST source supplies a defined radiative approximation, but its absolute atomic accuracy is not established here, and the final 2005 article's individual transition methods remain unverified. Finite-density populations, emitted power, and cooling still require independent numerical checks before physical acceptance.
