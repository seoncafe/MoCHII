# Native validation of the prepared 77-stage rate data

The actual linked species setup read eleven seven-stage chains and the
unchanged S II product-level dependency. The native entry point accepted a
100 eV photon band, including bare C VII. Three temperatures (5,000, 10,000,
and 30,000 K) covered 66 reactions each. Direct/cached populations, individual
and summed charge transfer, electron charge accounting, and the published or
adopted mean kinetic heats passed the driver's checks. Root verified source
and product hashes. See `result.json` and `root_review.json`.

The initial input directory contained only rate files and lacked S II, needed
by the sulfur product-term reaction. Its failure record is preserved in the
adjacent `staged77_charge_transfer` directory. No production file was replaced.
There was no matched plasma model, missing-line-atom validation, or claim of
published atomic accuracy. Extra reaction defaults still require the complete
v2 invocation audit. This is prepared rate coverage, not completion of the goal.
