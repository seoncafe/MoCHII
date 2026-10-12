# Si I installation and physical review

The adopted Si I electron-excitation model is installed in `data/atomic/nlevel_si_1.txt` and `data/atomic/cooling_si_1.txt`. It has 35 bound levels, 250 spontaneous pairs, 595 collision rows, and 8,925 original scaled collision nodes. Production remains at 43 metal stages and 32 adjacent transitions; the task retains the staged 64/53 design. Direct file inspection now finds 36 ordinary metal cooling fits, with seven missing: Ar I, Ca I, Ca III, Cl I, Cl V, Fe I, and Mg I.

## Physical judgment

The official original NIST 2002 MCHF line strengths identify the radiative calculation independently of the inaccessible final 2005 article. Identified length-form E1/E2/M1/M2 strengths use the selected observed photon frequencies and the actual emitting state's statistical weight. This corrects 15 doubled E1 rates and the reversed 14/15 transition, adds 21 missing pairs, and includes identified M2 contributions on 11 existing E1 pairs. The physical 15-to-14 rate is 5.631292575070227 s^-1. The 44 existing non-E1 CHIANTI rates remain explicitly approximate; their individual multipole provenance is incomplete. Absolute atomic accuracy and the final 2005 Methods are unverified.

The adopted constant-Upsilon floor is 1003.02 K, inferred from rounded first finite BT nodes. Maxwellian factors use the actual gas temperature. The numerical application interval is 1,000--100,000 K; numerical success does not establish the atomic source's accuracy throughout that interval. The conditional length/velocity calculation gives a maximum cooling difference of 18.8905642% at density 1e10 cm^-3 and temperature 3845.92 K. This operator diagnostic is not an atomic uncertainty interval.

## Matched numerical inputs and measured checks

| Condition | Compiled calculation | Independent balance calculation |
| --- | --- | --- |
| Level energies and weights | Identical selected 35-state model | Identical selected 35-state model |
| Spontaneous probabilities | Identical corrected 250 pairs | Identical corrected 250 pairs |
| Electron strengths | 595 original BT rows, 8,925 nodes | Same rows and nodes |
| Temperature | 401 points, 1,000--100,000 K | Same points |
| Electron densities | 1e-6, 100, 1e6, 1e10 cm^-3 | Same densities |
| Radiation and other colliders | None; optically thin spontaneous emission | Same assumptions |
| Spatial domain, gravity, boundary conditions | No spatial plasma model | No spatial plasma model |

Maximum relative cooling disagreement is 1.08358e-13. Maximum normalized energy-balance error is 8.64013e-14. All expected controller outputs and 63 frozen input hashes were checked before generator edits. Root visually inspected the complete [cooling profiles](../si_i_radiative_completion/cooling_profiles.pdf), [numerical differences](../si_i_radiative_completion/native_differences.pdf), all six pages of [gauge and population profiles](../si_i_radiative_gauge_sensitivity/full_profiles.pdf), and the [fit curve](candidate/fit/SiI_fit.pdf).

The seven-term serialized fit has maximum relative error 0.000673568012 on 801 independent temperatures, or 0.0673568%. The actual GNU `cooling_load/cooling_eval` reader differs by at most 4.44089e-16. The Intel registered table has 64,561 cells and 40 diagnostic samples; its maximum application error is 0.00136932578, or 0.136933%. The installed atom was opened by the current native module driver at 10,000 K and all four densities; cooling and 35 populations exactly match the saved validated profiles. The installed cooling fit was also opened through the actual coefficient reader at all 801 temperatures and reproduced the validated coefficients exactly. All 35 states have positive-energy spontaneous paths to the unique ground state. This checks actual file reading and downstream population, cooling, and emissivity routines, not a full spatial calculation or 64-stage activation.

The [installation manifest](installation_manifest.json), [main review](main_source_review.json), and [installed reader result](installed_open_check/result.json) distinguish file installation, validated numerical records, provenance annotations, and actual native invocation. Only installation-status and fit-metadata comments differ from the validated products. Unrelated fits were not recomputed. The aggregate coefficient listing now contains 37 species including H I and 240 terms.

## Existing defect corrected separately

The coefficient collector previously opened its output before formatting every input. The private Si I fit used a different metadata header, which caused a TypeError and left a partial aggregate listing. The production header now records the canonical range and the measured 801-point error. The collector validates metadata before writing and replaces its output atomically. A deliberately incomplete input was rejected while preserving an existing completed listing. The full listing was then collected from existing primary coefficients; no numerical fit was regenerated.

## MOCASSIN v2 comparison boundary

The inspected local `mocassin.2.02.73.2` lists `data/si1.dat` in `fileNames.dat`, but that file is absent. Its `hydro_mod.f90::makeElements` opens the listed atomic file when reading the selected element. A file-list entry does not demonstrate an available Si I excitation model or a successful calculation. The installed MoCHII Si I atom therefore adds an available numerical excitation source relative to that local data set. Electron-impact excitation and spontaneous cooling themselves are within the process class MOCASSIN v2 treats; this installation does not add a new physical process beyond that class. Auger joint charge/energy treatment and justified terminal-stage truncation remain separate unfinished requirements.

## Shared generation path

`chianti_cooling.bound_excitation_atom` invokes the source-identified Si I probability transformation after bound-state selection. The copied original NIST scientific file is byte-identical to the inspected source. The parser and transformation reside in `tools/fitting`, with no dependency on the private audit directory. `write_excitation_atom` records both scientific and selection-source hashes; `fit_cooling` now includes Si I. Root reran the focused source identity, photon direction, doubled-strength, shared-multipole, unchanged-collision-node, and serialization checks. Generated numerical atom records are byte-identical to the validated candidate, and the current shared low-density calculation exactly matches all 801 saved source values. Existing production numerical products were preserved. See the [shared generator review](shared_generator_review.json). Exact shared source bytes used by the completed controller are retained in the [executed-source archive](../executed_source_archive/record.json), because integration necessarily changes the current generator hashes.
