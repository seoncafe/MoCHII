# Approved excitation models: production selection and validation

The October 11, 2026 approval has been implemented. The default input is
`par%metal_excitation_model='mocassin_v2'`; `extended` explicitly selects the
larger previously installed models. The obsolete `fe_levels_full` input has
been removed. Rates remain in the common atomic directory. Excitation data
are selected from `excitation/<model>/`, including matching cooling and
optional hydrogen-impact data. The same selection reaches line emission,
registered thermal cooling, and charge-transfer product populations.
Missing model directories and unknown selections are rejected without
silently selecting other data. Ca VII remains absent pending physical correction.

Native readers loaded all 55 reference atoms (571 levels) and all 61 extended
atoms (1,946 levels), including the approved 142-level Fe II atom with 1,184
radiation rows and 10,011 collision rows. These are installed-file counts,
not active-ion counts. Production rate chains still contain 43 metal stages
and 32 adjacent transitions; activation of the 77-stage chains remains pending.

The new `UPSILON_CONSTANT_UPPER_K 100000` record after the collision rows
preserves original interpolation for C I, S III, and Ar VI. Above their last
source node it holds only Upsilon constant through 100,000 K; Maxwellian
factors use the actual temperature. Mg I retains the approved single-source-
node constant-strength representation. Atomic accuracy beyond the published
nodes is not established. Source curves were not refitted or regenerated.
The cold-only exponential fits for the three extended endpoint models were
not installed as matching hot cooling; these ions use absolute atom cooling.

The complete production build passed. Actual input selection, all installed
atom readers, and 3,216 direct population/cooling/emission samples passed.
The complete profiles exactly reproduce the previously inspected candidate
profiles at stored precision. Eighty registered cooling samples passed with
maximum relative direct difference 2.9171688e-5 (0.0029172%). Six invalid-
selection, missing-directory, and above-temperature-limit cases were rejected.
`root_review.json` verifies the source and executable hashes.

No spatial nebular comparison, full species thermal equilibrium validation,
or universal atomic accuracy is claimed. A preceding installation attempt
incorrectly checked an unused excited C II state against the hydrogen-impact
term mapping. The corrected check uses only the states actually addressed by
those hydrogen-impact records. Its failed staged products and program are
preserved in `../installed_excitation_models/`.

## Final common-hydrogen path correction

A caller audit confirmed that cooling_mod also registers H I through the
shared table API. The selector now applies only to metals and preserves the
common H I filename. The final production build passed, and the actual file
selection and reader loaded the existing 25-level H I model with no numerical
table regeneration. Five installed atoms and six malformed records were also
checked through the current parser. The sources executed in the earlier
full-profile validation are preserved in this directory; the subsequent diff
changes only the common H I selection and a defensive empty-row condition.
See [the final reader report](../approved_excitation_reader_guards_checked/README.md).
