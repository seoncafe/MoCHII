# Stationary Low-Density Cooling Validation

Date: October 9, 2026 (KST).

## Physical judgment and implementation

Assigning every ion to its ground state is not the stationary low-density
limit of the supplied C III and Si III atoms. Each has an excited state with
no outgoing spontaneous transition. Its population remains finite as electron
density tends to zero and must be obtained from collision balance, including
radiative cascades following excitation. The generator now solves that
balance. This statement describes the supplied electron-collision and
spontaneous-radiation model; it does not include hyperfine decay, additional
two-photon transitions, photoexcitation, or ionization of excited states.

For a collision from stable state i to level u, let B(u,j) be the probability
that spontaneous cascades end in stable state j. The effective rate into j
is q(i,u) B(u,j). Stationary populations satisfy the resulting rate matrix
and normalization. The emitted energy is E(u)-sum_j B(u,j) E(j). Summing that
energy over collisions and their source populations yields Lambda. A unique
radiative end state reduces this expression to ground-state excitation power.

An independent test assembles the full finite-density level rate matrix,
including spontaneous transitions and both collision directions, and sums
its emitted power. Actual C III and Si III agree with the limiting calculation
at 1,000, 10,000, and 100,000 K for n_e=1e-12 cm^-3 within 1e-8 relative
accuracy. A synthetic three-level atom tests transfer between two stable
states; C II tests reduction to the single-state expression.

## Matched inputs and limits

This is a rate-coefficient comparison, not a nebular-model comparison.

| Condition | Stationary source | Independent balance / fit comparison |
|---|---|---|
| Atomic source | Same CHIANTI 11.0.2 `.elvlc`, `.scups`, `.wgfa` | Identical files; source SHA-256 records retained |
| Level set | All collision-connected states reachable from ground | Full supplied collision/radiation state set for finite-density tests |
| Electron distribution | Maxwellian | Same distribution and collision-strength descaling |
| Temperature | 1,000-100,000 K | Three balance points; 801 temperatures for serialized fits |
| Electron density | Analytic n_e -> 0 | 1e-12 cm^-3 in full balance; this deliberate mismatch measures approach to the limit |
| Radiation | Spontaneous emission; optically thin | Same; no external photoexcitation |
| Spectral integral, gravity, spatial domain, boundaries | Not applicable to a local cooling coefficient | No nebular transport model was run |
| Ionization and composition | One specified ion; coefficient normalized to n_e n_ion | No ionization-equilibrium or abundance comparison is claimed |

## Measured results and changed products

| Ion | Maximum serialized fit error on 801 temperatures | Maximum Fortran/Python difference | Physical source change at 10,000 K relative to imposed ground population |
|---|---:|---:|---:|
| C III | 0.127901% | 1.08913e-13 | -0.593587% |
| Si III | 0.438816% | 3.04423e-13 | -1.146694% |
| S V (staged) | 0.278725% | 1.84408e-13 | Single stable-state source limit |

The physical source correction spans -0.847669% to -0.422925% for C III and
-3.337890% to -0.461394% for Si III across the full temperature range.
Differences from old production fits also include their previous fitting
error: corrected source / old fit - 1 spans -5.9953% to +2.6094% for C III
and -24.3069% to +14.1644% for Si III. These are coefficient differences,
not measured nebular temperature or spectrum changes.

Both production cooling files and their two aggregate-table blocks are
corrected. Old coefficients are preserved as `previous_production.txt` beside
the staged files. Existing model outputs are preserved and stale until
recalculated. Unrelated cooling files and n-level products were not regenerated.

The existing exponential file layout is preserved. Positive-only expansions
failed the 1% accuracy requirement for C III. Signed coefficients provide an
accurate mathematical expansion of total cooling; individual terms are not
emission rates. The total is positive across the fit interval. Additional
checks over 1-10^9 K found no negative sum, but do not establish validity or
accuracy outside the fit range.

Full temperature curves and residuals were visually inspected; labels and
curves fit within the figure boundaries:

- [C III source and serialized fit](validation_2026-10-09/stationary_cooling_exponentials/products/CIII/CIII.pdf)
- [Si III source and serialized fit](validation_2026-10-09/stationary_cooling_exponentials/products/SiIII/SiIII.pdf)
- [S V source and serialized fit](validation_2026-10-09/stationary_cooling_exponentials/products/SV/SV.pdf)
- [Compiled reader/evaluator results](validation_2026-10-09/sv_conditioned_exponentials/products/fortran/result.json)
- [Stationary population correction](validation_2026-10-09/stationary_cooling_exponentials/products/stationary_population_comparison.pdf)
- [Population-comparison numerical records](validation_2026-10-09/stationary_cooling_exponentials/products/stationary_population_comparison.json)
- [Production-file update and hashes](validation_2026-10-09/stationary_cooling_exponentials/products/production_cooling_update.json)

## Source defects and remaining work

The official CHIANTI `pro/read/read_elvlc.pro` reads
`(i7,a30,a5,i5,a5,f5.1,f15.3,f15.3)` from the first 87 characters.
The old Python parser read the last three whitespace fields and misidentified
N IV's extra g, n, and L metadata as J and two energies. The corrected reader
obtains level 3 as g=3, E_obs=67,272 and E_th=67,458 cm^-1. Tests cover main
records, comma-separated extra columns, and whitespace-separated extra columns.

Even after correct parsing, N IV has 14 positive radiative rows whose adopted
energies reverse the supplied downward transition. For example, source levels
107 and 116 have observed energies 595,296 and 593,690 cm^-1 although the
source gives a positive 116 -> 107 spontaneous rate. The generator rejects
this multi-stable-state atom rather than emit a physically inconsistent fit.
External source files remain unchanged. Review of a consistent energy set is
required before N IV activation. Some existing single-end-state atoms also
contain inconsistent intermediate energies; their aggregate excitation-power
identity does not use those intermediate gaps, but detailed emission needs a
separate energy consistency check.

S V generation, the 801-point fit check, and the optimized Fortran evaluation
are complete. Its 14-term staged fit has 0.278725% maximum error. A previous
11-term candidate had large cancelling coefficients and failed the unchanged
1e-10 Fortran/Python check with a 1.34094e-9 difference. Removing nearly
linearly dependent columns with the SVD cutoff and checking the sum condition
produced a stable fit that passes at 1.84408e-13. The failed candidate is
preserved in `sv_fourteen_terms/products/failed_sv_fit/`.

S V source includes an extremely small positive spontaneous rate; its strict
zero-density limit does not establish a useful finite-density reference.
Density and selected-level convergence remain required before activation.
An additional source audit reads all 20 added CHIANTI atoms: nine have
radiative rows with reversed adopted energy ordering, and 14 have some
nonpositive collision gaps. Zero gaps can represent degenerate levels;
these counts require transition-specific review, not automatic rejection of
every listed atom. [Detailed energy audit](validation_2026-10-09/added_energy_order_audit.json)
records the rows. Intermediate gaps do not enter the single-end-state
aggregate cooling identity, but detailed line balance needs consistent levels.

No additional production ion stage is active. General density dependence,
reduced-atom convergence, and matched photoionization-model profiles remain
unverified.

## Calculation records

Generation is controlled by `~/.codex/scripts/calculation_workflow.py` in tmux
with recorded commands, source hashes, process IDs, exit status, and numerical
JSON acceptance checks. Failed attempts are preserved in `stationary_cooling/`
and `stationary_cooling_checked/`. The C III and Si III steps in `stationary_cooling_exponentials/` passed;
its S V step failed the original eight-term accuracy requirement. The waiting
`stationary_cooling_final/` correctly refused further checks after that
failure. A 14-term search in `sv_fourteen_terms/` passed fitting but failed
Fortran precision because of cancellation. The accepted S V recovery and
all three compiled evaluations are recorded in
`sv_conditioned_exponentials/`. An earlier direct Fortran check started
before S V existed; its C III/Si III samples are preserved in
`stationary_cooling_fortran/`, and the missing-file failure is not counted as
a passed sequence. All controllers have stopped. No calculation or nebular
model is running; no model-event delivery or automatic scientific judgment
is claimed. Generated dates now use the actual KST date. Production energy
headers were corrected without changing coefficient rows; update records
retain both the validated numerical-file hashes and final production hashes.
