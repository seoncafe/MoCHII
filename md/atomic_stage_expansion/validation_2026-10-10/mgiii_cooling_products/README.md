# Prepared Mg III cooling product validation

This workflow is **prepared only**. It must not be launched until the corrected
level-167 source identity, NIST/CHIANTI ledger, and new 801-temperature native
comparison have been reviewed and the root session authorizes the product run.
The former Mg III candidate with SHA-256
`46e71f59c1c885a46ec4304f4a26bc1549e716a5490c6118d4cd242b16711873`
uses the physically wrong 628148.375 cm^-1 value for level 167 and is
explicitly rejected. The corrected atom must use the CHIANTI theoretical
634430.000 cm^-1 value for that state; no unrelated level energy is shifted.

The root session supplies a new candidate path and a JSON validation manifest
after the corrected physical and numerical checks. The manifest must contain:

```json
{
  "passed": true,
  "candidate_sha256": "SHA-256 of the corrected complete candidate",
  "temperature_points": 801,
  "densities_cm3": [1e-6, 100.0, 1e6, 1e10],
  "maximum_relative_cooling_errors": [0.0, 0.0, 0.0, 0.0],
  "level_167_identity_checked": true,
  "ledger_numerical_checks_passed": true
}
```

The four error entries must be measured native-versus-independent maxima
below `1e-12`; the example zeros are schema placeholders, not measured
values. The manifest and candidate are read and SHA-256 guarded when
`prepare_plan.py` runs, after approval and final shared source edits. The plan
also guards CHIANTI source files, the selected ionization source, fitting and
writer code, native source, executable files, and linked cooling reader object
files. It never freezes the currently changing source beforehand.

On authorization, launch from a persistent terminal with environment values
that name the reviewed files:

```sh
MGIII_PRODUCTS_APPROVED=yes \
MGIII_CORRECTED_CANDIDATE=/absolute/path/to/corrected_atom.txt \
MGIII_VALIDATION_MANIFEST=/absolute/path/to/validation.json \
bash md/atomic_stage_expansion/validation_2026-10-10/mgiii_cooling_products/launch_controller.sh
```

The script writes an immutable `plan.json` and starts
`calculation_workflow.py` in the `mochii_mgiii_cooling_products` tmux session.
No calculation is started by merely reading or importing this directory.
For a separately authorized outer controller that already runs inside tmux,
`run_controller_sync.sh` accepts the same environment values and runs the
product controller synchronously. It returns success only after its three
steps and final JSON checks complete. This avoids a detached nested tmux
process that would let an outer successor advance too early.

The finite sequence has three stages:

1. Serialize the full 167-level source atom with the final writer. Require
   3,364 radiative rows, 13,861 collision rows, 207,915 original nodes,
   `BT_TMIN_K 1800`, and exact equality of all noncomment numerical records
   with the newly validated candidate. Its staged header documents the
   unchanged first-27 A values, the higher E1 photon-frequency-cubed and
   upper-weight correction, the 112 unidentified `gf=0` effective-A values,
   and the source-accuracy limitation. Stop on any difference.
2. Run `fit_cooling.py --ion-spec mg:3` into private staging. Compare the
   serialized coefficients against the full source on 801 temperatures from
   1,000 to 100,000 K, require a maximum difference below 1%, and produce a
   full-profile PDF. Validate `cooling_load` and `cooling_eval` in the linked
   native reader on the same 801 temperatures. The cooling reader is the
   existing isolated `bt_floor_cooling_products/build_mpi/check_cooling_fit.x`;
   its build record, executable, source, and linked-object hashes are recorded.
3. Use the existing optimized isolated oxygen registry driver to build and
   test the private 401-by-161 table (64,561 cells). Require 40 sampled temperature/density
   conditions, exact-grid error below `1e-10`, interior and application
   differences below 0.5%, the correct 167/3364/13861 source counts, and the
   same atom SHA-256. Stop if the actual compiled application differs.

Each child has an 8 GiB address-space limit, one-thread BLAS environment,
CPU limit, wall timeout, and an expected JSON product checked by the durable
controller. Failed or interrupted products are preserved. The sequence never
copies an atom or fit to production. Samples below the published 1,800 K
collision lower bound test the **adopted constant-Upsilon continuation** only;
their numerical agreement is not a claim of source accuracy. The final source
and fit plots, coefficient file, and native differences require physical review
before any installation.
