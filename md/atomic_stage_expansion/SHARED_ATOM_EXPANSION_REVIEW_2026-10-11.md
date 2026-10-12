> Current status: approved on October 11, 2026 and implemented with scoped
> validation. The proposal below is historical. See the
> [current implementation](validation_2026-10-11/approved_excitation_models/README.md)
> and [comparison](MOCASSIN_V2_MOCHII_BIDIRECTIONAL_COMPARISON.md).

# Reviewable expansion of six shared-ion atoms and matching cooling fits

## Physical judgment and proposed change

Replace each of six smaller installed atoms together with its matching
low-density cooling fit. Replacing only the atom would mix MOCASSIN collision
strengths with existing CHIANTI-based thermal coefficients. That combination
would not conserve the same line energy calculated by the new finite-level
model. The original six-atom-only proposal is superseded by this twelve-file
proposal. This is a prospective inconsistency in the candidate combination;
the six production atoms have not been replaced. These legacy source data are
not established to be more accurate than modern CHIANTI data.

| Ion | Atomic file and cooling file in `data/atomic/` | Levels before -> proposed |
| --- | --- | --- |
| C III | `nlevel_c_3.txt`, `cooling_c_3.txt` | 5 -> 10 |
| N III | `nlevel_n_3.txt`, `cooling_n_3.txt` | 5 -> 10 |
| O I | `nlevel_o_1.txt`, `cooling_o_1.txt` | 5 -> 7 |
| Ne III | `nlevel_ne_3.txt`, `cooling_ne_3.txt` | 5 -> 9 |
| S IV | `nlevel_s_4.txt`, `cooling_s_4.txt` | 5 -> 12 |
| Si II | `nlevel_si_2.txt`, `cooling_si_2.txt` | 5 -> 12 |

Existing input syntax and rate chains remain unchanged. Preserve all twelve
original files and record hashes before any authorized replacement. Fe II is
excluded because its existing forbidden transitions need the separate pending
physical choice. `metal_freefree=.true.` remains the user-requested exception.
The reference model interval for these six pairs is 1000-100000 K; installation
must verify that the invoked calculation respects the source bounds.

## Completed numerical prerequisite

The production cooling table now uses shape-preserving cubic interpolation of
ln Lambda in -1/T, with the original 401 temperature nodes. The full build,
797 compiled interpolation samples, 1584 actual thermal samples, 1584 emitted
energy samples, and 24 zero-density checks passed. Root verified recorded
source and artifact hashes and inspected all six full-profile figure pages.
See [the interpolation report](validation_2026-10-11/cooling_temperature_pchip/README.md)
and [root review](validation_2026-10-11/cooling_temperature_pchip/root_review.json).
This corrects the earlier N III interpolation failure without increasing the
number of temperature solves. It does not certify atomic accuracy.

The candidate preparation verifies unchanged retained weights and positive
radiative transitions. Maximum retained energy differences are 0, 0.002,
0.038, 0.002, 0.030, and 0.040 cm^-1 in the ion order above. The original
six atoms and six fits remain unchanged, verified against preserved hashes.

## Source-consistent candidate comparison

The matching fit files already passed the actual native fit reader checks in
[the 52-fit result](validation_2026-10-11/mocassin_cooling_native_exact_endpoints/native/result.json).
The candidate prediction rescales measured thermal coefficients by the ratio
of matching to existing low-density fits. This follows the invoked routine's
linear dependence on the fit. Those first profiles are derived values. A subsequent actual combined native
invocation completed and independently confirmed the predicted differences.
The 792 predicted samples are compared with directly computed full-atom emitted
energy; all six satisfy the 0.5% representation criterion.

| Element and stage | Proposed levels | Maximum mixed-source energy discrepancy | Maximum matched-source energy discrepancy |
| --- | ---: | ---: | ---: |
| C 3 | 10 | 31.424768% | 0.021045% |
| N 3 | 10 | 35.476949% | 0.021010% |
| O 1 | 7 | 2.153019% | 0.049893% |
| NE 3 | 9 | 9.827943% | 0.143944% |
| S 4 | 12 | 34.960954% | 0.022069% |
| SI 2 | 12 | 78.420181% | 0.457773% |

At 10000 K and ne=1e-6 cm^-3, the C III mixed-source thermal coefficient is
0.6993497370437626 times its new atom's emitted-energy coefficient. This is
an isolated coefficient comparison, not a measured change in a nebular model.
The largest matched-source discrepancy is 0.457773%, for Si II. This is a
sampled numerical comparison, not a continuous bound or atomic accuracy claim.
See [the candidate result](validation_2026-10-11/shared_atom_matched_cooling_checked/products/result.json),
[full profiles](validation_2026-10-11/shared_atom_matched_cooling_checked/products/matched_cooling_profiles.pdf),
and [root review](validation_2026-10-11/shared_atom_matched_cooling_checked/root_review.json).

| Condition | Actual fixture / derived prediction | Independent reference |
| --- | --- | --- |
| Atom | Supplied full MOCASSIN candidate | Identical states, A values, collision strengths |
| Temperature | 33 values, 1000-100000 K | Identical values |
| Electron density | 1e-6, 100, 1e6, 1e10 cm^-3 | Identical values |
| Cooling fit | Matching source fit divided by existing fit rescales measured coefficient | Direct emitted-energy matrix sum |
| Hydrogen density and abundance | 100 cm^-3 and 1e-4 | Divided out to obtain coefficients |
| Stage fractions | Current production chains with imposed photo rates | Not independently validated; divided out |
| S II dependency | Existing 70-level atom, fixture-only zero S II cooling to isolate S IV | No claim of altered production S II |
| Spatial domain, velocity, gravity, spectrum, boundaries | No spatial model | No matched nebular comparison |

## Approval and implementation acceptance

Automatic approval review rejected the earlier installation before execution:
it found insufficient explicit approval to overwrite six public atomic-data
files. Candidate-only preparation did not execute that rejected action.
The governing instruction in
[`~/.codex/AGENTS.md`](/nfs/mocafe/kiseon/.codex/AGENTS.md) states:
"Ask before changes to public interfaces, data layout, user input parameters,
or other changes costly to reverse."

The [paired installation preparation](validation_2026-10-11/shared_atom_paired_installation_preparation/README.md)
preserved all twelve originals without changing production files. The updated
tool rejects atom-only production replacement.

Approval is now requested for the concrete twelve-file paired replacement
above. The actual paired thermal/emission consumer and affected O I hydrogen-option
checks have now passed in isolated fixtures. See [the actual paired result](validation_2026-10-11/shared_atom_matched_native/README.md).
After approval, require matching source hashes, preserve originals, verify
installed temperature bounds and file readers, and check affected consumers
using the installed paths. Do not regenerate
unrelated products. The global 77-stage objective remains incomplete.

## Preserved unsuccessful diagnostics

The original linear 401-node candidate failed N III's 0.5% interpolation
criterion. A private 801-node diagnostic reduced that error but first exceeded
a finite S II registration limit. Its later continuation produced numerical
outputs, then failed its frozen-source check when production interpolation was
edited. Those controllers remain failed and their products are preserved.
The completed PCHIP controller uses its own verified source and outputs;
no earlier failed workflow is relabeled as successful.
