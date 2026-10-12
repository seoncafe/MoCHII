# Prepared Ne I cooling product validation

This private workflow is prepared, not launched. The root session will review
the corrected Ne I source and native profiles before starting the product
controller. Its input is a complete 31-level version-2 atom with 144
spontaneous rows, 465 collision rows, 6,975 original BT nodes, and
`BT_TMIN_K 1000.36`. The two added spontaneous channels are
`A(3,15)=9.00e5 s^-1` and `A(5,15)=6.82e7 s^-1` from the NIST-dated
Stout Ne I source; all 142 original CHIANTI radiative rates remain unchanged.
The earlier 31/142 atom and its fit are physically stale and remain
preserved in their existing directories. This workflow rejects that atom's
SHA-256.

The root session supplies the corrected candidate and
`ne_i_radiative_completion/validation_manifest.json`. The manifest schema is:

```json
{
  "passed": true,
  "candidate_sha256": "SHA-256 of the complete corrected candidate",
  "temperature_points": 401,
  "densities_cm3": [1e-6, 100.0, 1e6, 1e10],
  "maximum_relative_cooling_errors": [0.0, 0.0, 0.0, 0.0],
  "maximum_relative_strength_error": 0.0,
  "stout_added_radiative_pairs_checked": true,
  "ledger_numerical_checks_passed": true
}
```

All five errors must be finite, nonnegative, and below `1e-12`. The example
zeros are schema placeholders, not measured results. The atom stage requires
exact equality of noncomment numerical records with the independently
validated candidate; it checks both added A values and source counts. Source
files, code, manifests, executables, and linked objects are hashed only when
the approved controller starts, after shared source edits and root review.
The workflow does not regenerate the earlier 31/142 products.

The synchronous successor interface for an already persistent outer
controller is:

```sh
NEI_PRODUCTS_APPROVED=yes \
NEI_CORRECTED_CANDIDATE=/absolute/path/to/corrected_atom.txt \
NEI_VALIDATION_MANIFEST=/absolute/path/to/validation_manifest.json \
bash md/atomic_stage_expansion/validation_2026-10-10/ne_i_cooling_products/run_controller_sync.sh
```

It returns success only after all three stages complete and
`candidate/result.json` has `passed=true`, `ion=NeI`. For a standalone
authorized launch, `launch_controller.sh` starts the same finite plan in a
tmux session. Neither script installs the atom or cooling coefficients.

The finite product sequence:

1. Write and compare the complete 31/144/465 atom with 6,975 nodes and its
   inferred 1,000.36 K BT floor. Stop on any numerical-record difference.
2. Fit `fit_cooling.py --ion-spec ne:1` privately. Compare the actual
   serialized coefficient file with the full bound source over 801
   temperatures from 1,000 to 100,000 K; require maximum relative error
   below 1%, a complete PDF profile, and the isolated GNU compiled cooling
   reader error below `1e-9`.
3. Build the private 401-by-161 density table (64,561 cells) with the
   existing isolated optimized registry executable. Require 40 sampled
   temperature/density conditions, exact-grid error below `1e-10`, interior
   and application errors below 0.5%, and matching atom hash and counts.

Each stage has an 8 GiB address-space limit, one-thread BLAS, CPU limit,
wall timeout, output checks, and durable state. Failed and interrupted
products are preserved. The source's published scattering method and atomic
accuracy in this interval remain unverified. The BT floor is inferred from
rounded source nodes, not a published physical lower bound. Numerical
agreement does not establish physical accuracy. Root will inspect the
profiles and products before any production installation.
