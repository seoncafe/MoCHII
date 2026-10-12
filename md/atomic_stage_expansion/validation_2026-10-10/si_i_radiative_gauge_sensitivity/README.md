# Si I electric dipole gauge sensitivity

This is a conditional comparison of two forms of the original NIST 2002 Si I MCHF electric dipole strengths. It does not establish an uncertainty interval or the accuracy of either operator. NIST cautions that the velocity-form Breit-Pauli operator omits relativistic terms that can matter when terms mix strongly. The final 2005 journal article's individual transition methods remain unverified in this audit.

The [runner](run_gauge_sensitivity.py) takes the frozen 35-state [length-form candidate](../si_i_electric_dipole_source_audit/candidate_v2/nlevel_si_1.txt) and the independently calculated length-form [profiles](../si_i_radiative_completion/direct_profiles.npz). It changes only the NIST length-form E1 contribution on 196 identified pairs to the corresponding original velocity-form line strength, using each selected positive photon gap and physical emitting statistical weight. Identified M2 contributions on 11 of those pairs remain unchanged. All 44 original CHIANTI non-E1 effective rates, 10 new non-E1 source pairs, selected energies, 595 collision rows, and the adopted 1003.02 K numerical effective-strength floor remain unchanged. The source file hashes and candidate hash are checked before calculation.

The planned calculation covers 401 logarithmically spaced temperatures from 1,000 to 100,000 K and four electron densities, `1e-6`, `100`, `1e6`, and `1e10 cm^-3`. It solves the same 35-state stationary balance for the velocity-form option. The output contract is:

| Output | Contents and completion check |
|---|---|
| `full_profiles.npz` | Complete temperature and density axes, both cooling profiles, and both 35-state population arrays at each density; loaded arrays must exactly match values passed to the writer. |
| `full_profiles.pdf` | Six pages: cooling profiles, relative cooling differences, and one complete population/difference page for each density; `pdfinfo` must report six pages. |
| `result.json` | `passed: true` only after both products pass checks, with source hashes, changed E1 count, maximum cooling and population differences at each density, elapsed time, product paths, and product hashes. |

The controller successor command, after the Si I direct-profile and energy-balance manifest succeeds, is:

```bash
MPLCONFIGDIR=/tmp/mochii_sii_gauge_mpl OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python md/atomic_stage_expansion/validation_2026-10-10/si_i_radiative_gauge_sensitivity/run_gauge_sensitivity.py
```

Freeze the runner, this README, the candidate module/atom/result, both original NIST files, the CHIANTI `.elvlc`, `.wgfa`, and `.scups`, `tools/fitting/validate_excitation_density_profiles.py`, and the prior `direct_profiles.npz`/`direct_step.json` after those prior files complete. Set a finite controller timeout and memory limit. The successor must preserve partial outputs and stop on any hash, array, population, or PDF failure; do not restart an uncertain partial run automatically. The result is a sensitivity diagnostic and cannot by itself accept or reject the physical source based on the size of its numerical difference.

Preparation checks passed for all 196 velocity-form E1 strengths and a small 3-temperature, 4-density finite-population probe. The full 401-temperature calculation and PDF inspection have **not** yet run.
