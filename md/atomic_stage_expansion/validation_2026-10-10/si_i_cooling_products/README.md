# Private Si I cooling-product successor

This workflow is prepared for the corrected, source-hash-frozen 35-level Si I
candidate in `../si_i_electric_dipole_source_audit/candidate_v2/`. It has **not
been launched**. The source/native/population/energy-ledger controller in
`../si_i_radiative_completion/` must first complete and produce its validated
manifest. Root review of the physical source model and native profiles is a
separate requirement from the numerical checks here. No production atom or
cooling file is written.

The synchronous entry point is `run_controller_sync.sh` with environment
`SII_PRODUCTS_APPROVED=yes`, `SII_CORRECTED_CANDIDATE=<absolute candidate>`, and
`SII_VALIDATION_MANIFEST=<absolute numerical manifest>`. It must run inside the
approved persistent outer controller after the source checks. Its final result
is `candidate/result.json` with `passed=true, ion="SiI"`; its own durable plan
and state are `plan.json` and `controller/state.json`. Any preexisting plan or
product causes a stop rather than an overwrite. A failed stage preserves its
files and logs.

The manifest must identify candidate SHA-256
`2c47d6e54cf4d1cfbc6340c5f01837889008a8fa8a0ac13afe0169f7f5d824df`,
35 levels, 250 radiation rows, 595 collision pairs, 8,925 source nodes, the
adopted 1003.02 K collision-strength floor, and 401 temperatures at four
densities. It must report finite relative cooling, population, collision
strength, and energy-ledger numerical checks. `source_accuracy_verified=false`
is retained explicitly. The exact required fields are checked in
`prepare_plan.py`, which freezes every source and linked-object hash before
product generation.

Three finite stages follow in order:

1. Copy the validated source atom byte for byte into a private input directory
   and verify its source identities and numerical records.
2. Fit the **corrected** full 35-state low-density source from
   `si_i_spontaneous_probabilities.py::source_atom()`, using process-local
   source overrides of the existing fitting functions. No shared file or
   public source dispatch is changed. Require a maximum 0.3% difference over
   801 temperatures from 1,000 to 100,000 K after coefficient serialization;
   save the complete source/fit/difference PDF. Check the actual compiled GNU
   coefficient reader and freeze its executable and linked objects. The reader
   calls `cooling_mod`, whose source and driver are checked against the build
   record; a later change in `nlevel_cooling_mod` does not enter this direct
   coefficient-reader check. The subsequent registry stage tests the current
   optimized n-level implementation separately.
3. Read the copied source atom and serialized fit through the actual optimized
   Intel n-level registry driver. Require 64,561 density-table cells, all 40
   native check samples, and maximum 0.3% interpolation/application errors.

The controller caps each stage at 8 GiB of virtual memory, one BLAS thread,
and finite CPU/wall-clock time (180 s, 3,600 s, and 1,800 s). Its shell runner
does not itself create a persistent session; the outer source controller
provides persistence. Fit accuracy, native agreement, and atom source
accuracy are distinct questions. The NIST 2002 radiative calculation and
original CHIANTI collision nodes are retained, while the final 2005 article's
transition-level methods and the atomic accuracy over this numerical
temperature range remain unverified.
