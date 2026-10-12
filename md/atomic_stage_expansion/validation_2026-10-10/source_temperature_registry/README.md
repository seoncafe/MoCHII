# Source-limited cooling temperature grids

Date: October 10, 2026. The cooling registry previously built every atom on 401 temperatures through 1,000,000 K. A version-2 `TEMPERATURE_TABLE` collision row may end earlier, and `nlevel_mod::nlevel_upsilon` rejects temperatures above its final source node. Registration of the staged Ar II atom, whose last Upsilon node is 158,000 K, therefore reached unsupported temperatures during setup. This change makes the registered grid end at the smallest final temperature among its temperature-table collision rows. BT-source and legacy atoms retain all 401 temperature rows.

The registry retains every original 0.01-dex log-temperature row below the source ceiling and appends the last source temperature itself. Its last interval is usually shorter than 0.01 dex. Runtime interpolation uses the stored ion-specific vector for both absolute cooling and the low-density reference. The electron-density columns, interpolation formulas, cooling combination rules, and public signatures are unchanged. Temperatures above a source ceiling now produce an explicit error instead of extrapolating Upsilon. The final source temperature is reconstructed from the source reader's stored `ln(T)` coordinate; a one-representable-value downward adjustment is made only if `exp(ln T)` rounds beyond the node. The original physical source point still passes the runtime logarithmic bound check.

[`run_checks.py`](run_checks.py) builds copies of `define.f90`, `nlevel_mod.f90`, and `nlevel_cooling_mod.f90` in its private `corrected/build/` directory with GNU `mpifort -O1 -fcheck=all -ffpe-trap=invalid,zero,overflow`; it does not rewrite production objects. Its source hashes and complete 24-row profiles are in [`corrected/result.json`](corrected/result.json). The staged [Ar II source-node file](../fine_structure_native/ArII/nlevel_ar_2.txt) and a two-level synthetic temperature table ending at **158,123 K**, which falls between registry grid points, were compared with direct level solves at four densities from 1e-6 through 1e10 cm^-3. Exact grid/source-node errors were at most 3.89e-15 and 3.44e-15, respectively; the largest errors within the shortened final interval were 4.04e-5 for Ar II and 2.63e-5 for the synthetic table. An additional synthetic ceiling of 999,900 K tested a shortened final interval while the registry still has 401 rows; its maximum error was 3.35e-5. An unlimited legacy atom retained the full-grid behavior with maximum relative error 2.89e-15, including its existing upper-edge clamp. Ar II and both limited synthetic atoms explicitly rejected temperatures above their source ceilings. NaN, infinity, and zero were rejected explicitly for both temperature and electron density without a floating-point exception signal.

The existing four subnormal and exact-zero fixtures were linked to the same isolated objects and exercised through `tests/ion_stages/check_small_excitation_cooling.f90`. The `ratio`, `partial`, `zero_ratio`, and `zero_partial` fixtures each completed all 48 direct/table/combined-value rows with runtime finite, positivity, zero-preservation, and agreement checks. Their outputs and hashes are recorded in [`corrected/subnormal_result.json`](corrected/subnormal_result.json). The top-level `result.json` and `subnormal_result.json` are earlier measurements of the initial candidate; `attempt1/` preserves a diagnostic run interrupted by a stage-argument error in the test script. Neither is evidence for the corrected source. This is targeted verification of the changed registry code path; the full project build and wider production-profile review are owned by the main session.

## Main-session production validation and Ar II installation

The actual `make DEBUG=1 -j4` target rebuilt successfully. A driver linked
to the rebuilt Intel production objects passed the four physical/limited-grid
cases, all nine invalid-input or above-source-limit checks, and all four
48-row positive-subnormal/exact-zero fixtures. The measured Ar II, synthetic
158,123 K, synthetic 999,900 K, and unlimited-atom maximum relative errors
were 4.04255e-5, 2.63089e-5, 3.34688e-5, and 2.88658e-15.
[The production controller](production_review/controller/state.json) completed
both declared steps; [actual native results](production_review/native_result.json)
record the complete values and linked-object hashes. Focused tests were
linked with `-no-ftz`; production compiler defaults remain unchanged.

The preserved Ar II two-level source atom and its existing fitted cooling
were installed with identical bytes, then opened and registered through native
consumers. The source ceiling is 158,000 K; the fit interval is 1,000--100,000 K.
No numerical product was regenerated solely for installation. The fit's
801-point reader comparison was 5.53e-14 relative. The isolated limited-cell
checks do not certify arbitrary-temperature table errors or a spatial model.
[Installation evidence](argon_installation/installation.json) includes hashes.
Higher Ar II excitation is not represented by this ground fine-structure pair.

For BT/legacy atoms whose source rows do not define a physical-temperature
ceiling, the existing upper registry clamp at 1,000,000 K is retained. This
is a software-grid convention, not validation of cooling outside the accepted
fit or source interval. Finite temperature-table sources that extend beyond
that grid retain the same grid convention within their supported source
interval, while requests above their source ceiling are rejected.
