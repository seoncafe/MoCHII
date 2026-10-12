> Current status: approved on October 11, 2026 and implemented with scoped
> validation. The proposal below is historical. See the
> [current implementation](validation_2026-10-11/approved_excitation_models/README.md)
> and [comparison](MOCASSIN_V2_MOCHII_BIDIRECTIONAL_COMPARISON.md).

# Shared-ion excitation model selection: measured gap and proposed design

## Physical judgment

Additional explicit excitation levels represent the same electron-impact and
radiative processes in a larger atomic model. They are not selected by the
additional-ion cooling option. Their thermal contribution cannot be assumed
small from the number of levels or their omission in the v2 files. The same
selected atomic model and collision source must determine thermal cooling,
level populations, and emitted lines. Mixing a reference-size atom with an
unrelated complete-source cooling fit is not a consistent reference model.

A direct inventory measured twelve installed shared-ion atoms larger than the
supplied v2 model; ten belong to current production rate chains. Production
still contains 43 metal stages and 32 adjacent transitions. Counts below were
read from actual installed files and the actual v2 file mapping, not copied
from the preceding coverage report. Counts do not prove matching physical
states or cooling accuracy.

| Element and stage | Supplied v2 levels | Installed levels | Production chain |
| --- | ---: | ---: | --- |
| C 1 | 6 | 42 | Active |
| C 2 | 10 | 53 | Active |
| C 4 | 15 | 24 | Active |
| N 1 | 13 | 50 | Active |
| N 2 | 15 | 58 | Active |
| O 2 | 13 | 196 | Active |
| O 3 | 15 | 177 | Active |
| O 4 | 10 | 147 | Inactive |
| O 5 | 10 | 166 | Inactive |
| S 2 | 5 | 70 | Active |
| S 3 | 5 | 53 | Active |
| CA 2 | 2 | 5 | Active |

See [measured inventory](validation_2026-10-11/shared_excitation_model_inventory/inventory.json)
and [root verification](validation_2026-10-11/shared_excitation_model_inventory/root_review.json).

## Existing invocation and controls

`species_setup` loads an exponential cooling file by element and stage, then
registers the finite-density atom with `nlevel_cooling_add`. That routine
selects the ordinary atomic filename, except for Fe II/III when
`fe_levels_full` selects `_full`. `lines_mod` independently repeats that
Fe filename rule. The hydrogen-excitation atom loaded in `species_setup`
also uses the ordinary filename. The current Fe exceptions have no hydrogen
collision table, so the last distinction does not establish an existing
H-impact defect; it would become unsafe in a general model selector.

`additional_line_cooling` selects extra ions, not the level count of shared
ions. O II and O III therefore still use their installed 196- and 177-level
atoms under default local-density cooling. The existing `fe_levels_full`
files declare only 16 levels for Fe II and 34 for Fe III. The name does not
mean the supplied 142-level v2 Fe II model, and is not a general complete-atom
control. Merely renaming or forwarding it would preserve an ambiguous input.

## Qualified charge-transfer model dependence

Actual native reverse-rate and mean-energy consumers passed 808 samples for
5-level and 70-level S II product atoms with identical sulfur reaction data.
The 4S, 2D, and 2P ranges and weighted term energies match both files, so these
particular candidates require no level-range remapping. Maximum source-model
reverse-rate difference is 44.64932194%; maximum mean-energy difference is
0.2184102926 eV. The experiment includes collision-data differences and does
not isolate additional levels or measure net cooling. The existing equations
correctly retain populations outside the three represented terms; an inaccurate
normalization comment was corrected. See [the actual comparison](validation_2026-10-11/sulfur_charge_exchange_models/README.md).
The public selection must consistently reach `cx_state_setup`; numerical
qualification does not approve or implement the new public input.

## Proposed input and data shape, not implemented

Use one explicit metal excitation-model input to choose a coherent atom/fit
pair. A proposed shape is `metal_excitation_model='mocassin_v2'` by default
and an explicitly requested `metal_excitation_model='extended'` for the
existing larger source models. The values describe source-model selections;
their numerical difference includes source-data differences, not only added
levels. This input is not currently present. Its exact public interface and
the disposition of `fe_levels_full` require approval before implementation.
The proposed replacement removes the Fe-specific control after confirming
its replacement semantics, rather than silently accepting conflicting inputs.

Use explicit source-model directories under `atomic_dir`, each containing
paired atom/cooling files and a source manifest. Rate-chain and charge-transfer
files remain common because this selection changes excitation physics, not
the ionic roster. A single selected model must be used by registration,
neutral-H excitation, and line emission. `cx_state_setup` also loads the
product-ion atom for form-9 charge-transfer populations, so it belongs to the
same selection. Its `CXSTATEK` level ranges must be rematched to physical terms
and checked before accepting a different ordering. The rate coefficients and
term definitions remain common; actual density-dependent reverse rates and
mean heat must be revalidated with the selected population model. Do not silently substitute the other
model when a requested file is missing. A validated atom can calculate absolute
cooling without a fitted representation, avoiding an unnecessary fit requirement
for Ne VI or Fe V whose exponential fits were rejected for conditioning.

Preserve the larger atoms and their matching existing cooling files before
moving or replacing paths. This document does not authorize those changes or
impose the proposed directory layout on the current twelve-file paired change.
The paired six-ion candidate is independently ready; the remaining general
selection is a separate incomplete part of the full objective.

## Temperature scope and source review

The original archived v2 `splint` holds endpoint strengths constant outside
multi-node tables; its one-node branch is also constant. Current MoCHII source
candidates instead reject unsupported upper temperatures. C I and S III end
at 20000 K, Ar VI at 40000 K, and Mg I has only 10000 K. This distinction is
confirmed in [the original source review](validation_2026-10-11/v2_collision_temperature_endpoints/README.md).
The default source-model choice must state whether it accepts a qualified hot
endpoint continuation or retains these bounds. Do not silently broaden the
accepted cold approximation. The proposed model selector must not claim v2
usable-temperature equivalence while rejecting an unacknowledged part of the
v2 calculation. No production upper-temperature policy changed. The
[private endpoint-model candidates](validation_2026-10-11/upper_temperature_models_checked/README.md)
now have numerical direct and registered-consumer verification. They retain
original in-range polynomials; a synthetic appended node would change prior
interpolation and is not the selected candidate representation.

## Implementation acceptance after the public choice

1. Match configurations, terms, J values, statistical weights, and transition
   endpoints for the reference coverage. Level-count truncation alone is not
   sufficient. Record corrected external-source errors separately.
2. Resolve the atom and cooling representation from one source selection in
   all invoked consumers, including charge-transfer product-term populations.
   Check omitted input and each explicit model value.
3. Validate full temperature/density cooling and population profiles with
   matching spectra, abundances, and fixed physical inputs where comparisons
   include ionization. Isolate source-model differences from new processes.
4. Check thermal cooling against the sum of emitted line energies, option
   selection in the actual line writer, cold/hot bounds, missing-file errors,
   and saved-state interpretation. Keep numerical and source-accuracy claims
   separate. Run a spatial comparison only with fully documented input matches.
5. Update implementation status only after production invocation is verified.
   Prepared files, a compiled library, or a callable solver alone do not finish
   this requirement. Existing figures and tables need no path-only regeneration.
