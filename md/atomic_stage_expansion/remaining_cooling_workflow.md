# Remaining coherent-source cooling workflow

This is a historical October 10 controller plan. Its present-tense execution
description below records the plan at launch; it does not mean the controller
is still running. The current production selection and incomplete 77-stage
goal are recorded in the [current result](../ATOMIC_IONIZATION_STAGE_EXPANSION_RESULTS_2026-10-11.md).

The user authorized continuation through the atomic-stage plan and subsequently
adopted the ten high-charge Langevin limits. The persistent tmux session
`mochii-remaining-cooling` executes
`~/.codex/scripts/calculation_workflow.py` with the saved
[plan](validation_2026-10-10/remaining_coherent_cooling/plan.json).
The actual progress is recorded in
[controller state](validation_2026-10-10/remaining_coherent_cooling/controller/state.json).

The sequence is N V, O IV, Ne V, S VI, Ar V, Ar VI, Ar VII, Mg IV,
Fe V, Fe VI, Fe VII and Ca VI. Each fit must pass its complete 801-point
source comparison with the existing 1% error target. Its linked native
cooling-reader/evaluator check runs immediately afterward and must produce
`passed: true` before the next ion begins. Each fitting step is limited to
300 seconds, each native check to 120 seconds; OpenMP and BLAS use one thread.
Source-code and linked-object hashes are guarded by the controller. The
[input-source record](validation_2026-10-10/remaining_coherent_cooling/input_sources.json)
checks all 36 CHIANTI files against the inventory created before this sequence;
each fitted-ion JSON independently records its input hashes.

Failure or timeout stops the sequence and preserves outputs; no uncertain
partial product is automatically restarted. This controller produces staged
coefficients, vector figures, dense source checks, and native checks only.
It does not activate production data or make a physical source judgment.
No model-resume event watcher is configured. A completion event does not
resume this chat or perform unregistered successor work. Figure inspection,
source accuracy and remaining ion/energy completeness require subsequent
assessment against the actual products.

N IV and S V were already validated separately; C V, N VI, S VII and Ne IV
are complete in [high-excitation validation](high_excitation_cooling_validation.md).
O V, Cl VI, and Cl VII were excluded from this original sequence while their
energy consistency and low-temperature source decisions were unresolved.
Subsequent checks established coherent O V/Cl VI energy vectors and adopted
constant lowest-temperature Upsilon for Cl VII; see
[ICFT validation](icft_energy_atom_validation.md) and
[implementation status](implementation_status.md). The original products
remain preserved; no zero-cooling substitution was used.

## Measured completion status

The controller completed all 24 registered fitting/native steps on October 10,
2026, at 09:09:57 KST. No calculation remains running in this sequence.
All 12 fits pass the unchanged 1% dense-grid source requirement; all linked
native checks pass. Each ion uses 801 temperatures from 1,000 to 100,000 K.
These are staged coefficient checks, not accepted finite-density atoms or
a nebular model comparison. The complete source curves and residuals of all 12 generated vector figures
were visually inspected in the subsequent continuation. No clipped labels,
overlapping plot elements or unexplained curve discrepancies were found.
These checks do not replace finite-density atom convergence or source accuracy
assessment.

| Ion | Maximum relative fit error (%) | Native evaluator |
|---|---:|---|
| NV | 0.9800655 | Passed |
| OIV | 0.5069555 | Passed |
| NeV | 0.1598017 | Passed |
| SVI | 0.4440276 | Passed |
| ArV | 0.1556865 | Passed |
| ArVI | 0.2215969 | Passed |
| ArVII | 0.6653188 | Passed |
| MgIV | 0.8119473 | Passed |
| FeV | 0.6312646 | Passed |
| FeVI | 0.5644621 | Passed |
| FeVII | 0.2738698 | Passed |
| CaVI | 0.0212163 | Passed |
