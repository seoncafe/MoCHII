# Registered cooling validation within collision-table bounds

The old diagnostic required 100-1000000 K even for atoms whose collision
strength table rejects higher temperatures. It also reported a fixed 401*161
cell count without measuring the shortened registered table. Both issues are
corrected in the diagnostic, without changing production physics or inputs.
Finite-temperature atoms now select actual regular cooling nodes and interior
cells within the minimum collision-table ceiling. The wrapper reports tested
temperature bounds and marks the cell count unmeasured instead of inventing it.

The persistent controller completed compilation and both numerical checks.
The six-level supplied C I candidate passed 40 samples from 100 to
19952.62314968883 K, within its 20000 K source ceiling; maximum interior
cooling discrepancy was 0.0010735306425899704. The existing five-level C III
passed 40 samples from 100 to 1000000 K, with maximum interior discrepancy
0.0014182223424196172. Exact-node errors were at most 4.45e-15. These are
numerical registration comparisons with the same atom's direct solver;
low/high extrapolated legacy fit accuracy is not established. The cold
constant-strength approximation remains explicit and this test does not
extend its physical validity.

Root verified the controller, source hashes, and completed artifacts. Neither
production atomic data nor unrelated reference products were regenerated.
See [C I](finite_ci/result.json), [C III](legacy_ciii/result.json), and
[root review](root_review.json). The actual tested temperature/density records
are preserved in the corresponding `registered.txt` files.
