# Ar II and Ca IV effective collision strength interpolation audit

The actual MOCASSIN v2.02.73.2 interpolation and the adopted Stout log-PCHIP interpolation differ by at most **0.2673% for Ar II** and **0.5063% for Ca IV** on the sampled 1,000–100,000 K range. These are differences in effective collision strength, not measured differences in total cooling, ion populations, or nebular observables. The larger Ca IV difference occurs near 34,475 K; Ar II peaks near 31,261 K. Both lie within the respective tabulated temperature intervals.

## Invoked methods and input conditions

The local MOCASSIN v2 `data/fileNames.dat` selects `data/arii.dat` and `data/caiv.dat`. Its `source/hydro_mod.f90:499–672` reads those files and converts temperatures to `log10(T)`. `source/update_mod.f90:1268–1312` invokes `equilibrium` for available ions and sums its forbidden-line power at lines 1100–1109. `source/emission_mod.f90:2012–2050` calls `spline(logTemp,qq,1.e30,1.e30,qq2)` followed by `splint(logTemp,qq,qq2,log10Te,qomInt)`. The original `source/interpolation_mod.f90:194–328` defines these operations: a natural cubic spline of **linear Upsilon** against `log10(T)`, constant below and above the source nodes. Its arrays and interpolated Upsilon are default `real` (32-bit on this build). This audit compiles that original module in isolation and calls its routines on the parsed source nodes; no replacement cubic boundary condition was assumed.

MoCHII `src/nlevel_mod.f90:238–250,478–531,584–596` stores `log(T)` and `log(Upsilon)` for `TEMPERATURE_TABLE LOG_PCHIP_T`, uses weighted harmonic PCHIP slopes, and holds Upsilon constant below the first node. The diagnostic uses SciPy `PchipInterpolator` on those same transformed coordinates, with no extrapolation above the last Stout node. This matches the stated interpolation formula but is a diagnostic calculation, not a run of the full MoCHII cooling path.

| Input | MOCASSIN v2 Ar II | Stout Ar II | MOCASSIN v2 Ca IV | Stout Ca IV |
|---|---:|---:|---:|---:|
| Ground/upper statistical weights | 4 / 2 | 4 / 2 | 4 / 2 | 4 / 2 |
| Upper energy, cm^-1 | 1431.58 | 1431.583 | 3118.2 | 3118.200 |
| Radiative probability, s^-1 | 0.0530 | 0.0530 | 0.543 | 0.543 |
| First temperature node, K | 1585 | 1580 | 10000 | 10000 |
| Last temperature node, K | 158489 | 158000 | 1000000 | 1000000 |
| Number of Upsilon nodes | 11 | 11 | 11 | 11 |
| Upsilon ordinates | Same 11 values | Same 11 values | Same 11 values | Same 11 values |
| Below first node | Constant Upsilon | Constant Upsilon | Constant Upsilon | Constant Upsilon |

The largest paired node-temperature rounding difference is 489 K for both ions, or 0.3165% of the Stout Ar II temperature and 0.3101% of the Stout Ca IV temperature. Ar II's 0.003 cm^-1 energy difference is about 2.1 parts per million in the line energy; statistical weights and A values agree in the parsed files. MOCASSIN and MoCHII use different population and cooling implementations, so this small input difference must not be conflated with a total physical cooling difference. The local MOCASSIN sources are `data/arii.dat` and `data/caiv.dat`; the selected Stout files are `ar/ar_2/ar_2.{nrg,tp,coll}` and `ca/ca_4/ca_4.{nrg,tp,coll}` under `/home/kiseon/CLOUDY/c25.00/data/stout/`. Raw file hashes are in `results.json`.

## Complete sampled profiles and boundary checks

The [full-profile figure](upsilon_full_profiles.pdf) shows Upsilon and signed relative differences throughout 1,000–100,000 K. The complete numerical profiles are [Ar II](arii_full_profile.csv) and [Ca IV](caiv_full_profile.csv). Each includes an 801-point logarithmic temperature grid, every source node within this thermal interval, and points immediately below and above both cold boundaries. Additional direct node checks cover **all 11 nodes** of both formats, including nodes above 100,000 K. At the MOCASSIN nodes, the original single-precision spline reproduces its Upsilon values to at most 1.15e-7 absolute; log-PCHIP reproduces Stout ordinates to floating-point precision.

| Maximum absolute relative difference on sampled grid | Ar II | Ca IV |
|---|---:|---:|
| Actual MOCASSIN format versus actual Stout format | 0.2673% at 31,261 K | 0.5063% at 34,475 K |
| MOCASSIN spline versus log-PCHIP at the **same MOCASSIN nodes** | 0.2676% at 31,261 K | 0.4835% at 34,475 K |
| PCHIP using MOCASSIN nodes versus PCHIP using Stout nodes | 0.0293% at 14,791 K | 0.1039% at 17,278 K |
| Actual formats below the first MOCASSIN node | 0.0125% just below 1,585 K | 0% from 1,000 to 10,000 K |
| Actual formats at common MOCASSIN nodes | 0.0268% at 15,849 K | 0.0953% at 15,849 K |

The signed curves can cancel locally; the rows above are separate maxima and must not be added. The two cold boundaries are distinct for Ar II: Stout begins rising above 1,580 K while MOCASSIN stays constant through 1,585 K. For Ca IV the cold boundaries coincide. Both constant cold continuations are adopted numerical approximations. The published Pelan and Berrington Table 3 starts at log10(T/K)=3.2 for Ar II and 4.0 for Ca IV; neither calculation establishes collision strengths down to 1,000 K for both ions. The final journal scan is at `references/Pelan_Berrington_1995_AAS_110_209.pdf` in this repository.

This comparison holds the 11 Upsilon values fixed and isolates interpolation and tabulated-coordinate rounding. It does not solve level populations or vary electron density, ion fractions, recombination sources, abundance, or radiation field. Thus no nebular agreement claim follows from these curves. A later two-level cooling comparison would need matched energy, A, temperature, density, and population equations; it should report the unavoidable source-energy and coordinate mismatches explicitly.

## Reproduction

Run `MPLCONFIGDIR=/tmp/mochii_mpl_moc_stout python3 md/atomic_stage_expansion/validation_2026-10-10/moc_stout_upsilon/compare_upsilon.py` from the repository root. The script compiles the unmodified MOCASSIN `interpolation_mod.f90` together with `fortran_spline_values.f90` using `mpiifort -O0 -no-ipo`, then writes the profiles, `results.json`, and PDF. It takes about three seconds here. Compilation and intermediate tables stay in `/tmp/mochii_moc_stout_upsilon`; it does not change MOCASSIN or MoCHII production objects.
