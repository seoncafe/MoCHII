# Documentation review and update

Date: October 11, 2026 (KST).

## Scope and outcome

All Markdown and TeX files in the primary working tree were scanned for
atomic-stage coverage, excitation-model selection, option defaults, and the
status of historical calculations. The scan includes ignored documentation
and editor checkpoints; it excludes Git metadata, agent configuration, and
separate `.claude/worktrees` checkouts. The final inventory and file hashes
are in [document_scan.json](document_scan.json). This is a topic-specific
review, not independent verification of every unrelated scientific statement.

Current guidance now distinguishes the active 43-stage/32-transition registry
from the prepared, inactive 77-stage/66-transition registry. The reference
excitation selection contains 55 atoms, 571 levels, 49 primary fits, and
763 terms; the extended selection contains 62 atoms, 1951 levels, 38 primary
fits, and 377 terms. Installed excitation files are not active stage counts.
Ca VII remains unresolved, and the physics objective remains blocked.

The current optional physics defaults are disabled, except for the explicitly
requested `metal_freefree=.true.`. The default excitation model is
`'mocassin_v2'`. Historical run outputs were preserved and qualified wherever
their old defaults or source selection could be mistaken for current behavior.
No full spatial MOCASSIN v2 comparison is claimed.

## Updated documentation

- Project guidance: `README.md`, `CLAUDE.md`, `data/atomic/README.md`, and
  `examples/README.md`.
- Current result report: `md/ATOMIC_IONIZATION_STAGE_EXPANSION_RESULTS_2026-10-11.md`.
- Historical plans and reviews: the atomic expansion plan and October 9
  review, `COOLING_LOCAL_NE_PLAN.md`, `HEI_10833_PLAN.md`,
  `HEI_CASCADE_DIFFUSE_PLAN.md`, and `METAL_COUPLING_PLAN.md`.
- Atomic expansion status: `implementation_status.md`,
  `MOCASSIN_V2_MOCHII_BIDIRECTIONAL_COMPARISON.md`,
  `remaining_physical_requirements_2026-10-10.md`,
  `remaining_chlorine_calcium_excitation_sources.md`,
  `remaining_cooling_workflow.md`, `mocassin_v2_physics_scope.md`, and
  `automatic_continuation_status.md`.
- Fitting guidance: `tools/fitting/README.md`.
- TeX documents and corresponding PDFs: the seven documents listed below.

Three He I example input comments were corrected from “default on” to
“explicitly enabled below”; their explicit input values did not change.
The `src/nlevel_mod.f90` header was corrected to describe reading the
hydrogen-impact file beside the selected atom, rather than implying that
metal hydrogen-impact files always come from the historical root directory.
This was a comment-only correction.

## Implementation evidence

The current defaults were checked in `src/define.f90`. The selected atom and
cooling paths were checked against `metal_excitation_data_file` in
`src/nlevel_mod.f90` and its consumers. `hcoll_load` constructs the hydrogen-
impact filename beside the selected atom. The fitting documentation separates
the requested generator output directory from the installed runtime selection.
Stage counts and excitation inventory values were read from current files;
saved physics checks are identified as saved checks in the result report.

## PDF builds and visual inspection

Each changed TeX document was built at least twice for references. The changed
pages were rendered and inspected for overlaps, clipped text, and material
outside the margins. The primary session also inspected representative
rendered pages of every changed PDF, including the User Guide page with a
layout warning. The corresponding PDFs were installed under `docs/`.

| Document | Visually inspected PDF pages |
| --- | --- |
| `MoCHII_physics.tex` | 1, 6, 7, 8, 13, 14, 24, 25, 81 |
| `MoCHII_UserGuide.tex` | 1, 6, 7, 15, 16, 19, 20, 23, 31, 36 |
| `MoCHII_fitting.tex` | 1, 3, 4, 5, 8; generator-path change inspected after its final rebuild |
| `MoCHII_cooling_analysis.tex` | 1, 20, 26, 27 |
| `MoCHII_code_comparison.tex` | 4, 21, 28, 30, 54 |
| `PHOTON_ENERGY_SAMPLING_CODE_COMPARISON.tex` | 1, 4, 6, 7 |
| `Draine2026_HeI10833_tau.tex` | 1, 2 |

`MoCHII_QMC.tex` and `RTML_review.tex` were inspected for the affected topics
and did not require changes. `Draine2026_HeI10833_tau.tex` and its PDF are
ignored by an existing Git rule, so their updates do not appear in ordinary
Git status output.

There were no TeX errors in the successful builds. The User Guide has a
7.77757 pt overfull vertical-box warning on page 6; visual inspection found
no overlap, clipping, or content outside the page margins. An initially long
hydrogen-impact path in the fitting document was split and rechecked.
PDF hashes, page counts, and build-log locations are recorded in
[pdf_checks.json](pdf_checks.json).

Previous PDFs were preserved under `/tmp/mochii_tex_update/original/` and
`/tmp/mochii_tex_review/`. These are temporary build records, not a durable
archive. The updated PDFs and this review are the durable project products.

## Local links and validation limits

The Markdown scan found eight broken relative file links in the saved
`SHARED_ATOM_EXPANSION_REVIEW_2026-10-11` editor checkpoint. They were
corrected to resolve from its `.ipynb_checkpoints` directory, and a notice
identifies the checkpoint as historical. No numerical records were changed.
The final [local-link scan](link_scan.json) checks file targets only;
external URLs and section anchors were not validated by that scan.

`git diff --check` was run after the updates. No physics calculations,
reference outputs, or full regression suite were rerun: the changes in this
task affect documentation and comments only. Existing physics validation
and its limitations are documented in the
[implementation result report](../../../ATOMIC_IONIZATION_STAGE_EXPANSION_RESULTS_2026-10-11.md).
