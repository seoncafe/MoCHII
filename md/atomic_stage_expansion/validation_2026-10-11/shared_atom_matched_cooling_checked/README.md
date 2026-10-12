# Matching atoms and line cooling from the same source

The finite controller completed candidate preparation. Production was not
modified. Twelve matching atom/fit files are ready in `products/matched_data/`.
Original cooling files are preserved in `products/previous_cooling/`; original
atoms remain preserved by the earlier candidate preparation.

The 792 derived coefficients rescale measured native thermal coefficients by
the matching-to-existing low-density fit ratio. Matching fit files were already
checked by the native reader. These are compositional predictions, not a new
combined invocation. Maximum difference from direct source-atom emitted energy
is 0.457773% for Si II. All six full-profile pages were visually inspected;
no clipping or overlapping text was found. Input conditions and individual
errors are in [the review](../../SHARED_ATOM_EXPANSION_REVIEW_2026-10-11.md).
See [results](products/result.json), [profiles](products/matched_cooling_profiles.pdf),
and [root verification](root_review.json). The preceding preparation controller
failed because it selected an obsolete partial fit result; it remains failed.
This corrected preparation uses the completed 52-fit exact-endpoint record
without regenerating those fits. Explicit paired replacement approval remains
pending, the [actual combined-consumer check](../shared_atom_matched_native/README.md)
has subsequently passed.
