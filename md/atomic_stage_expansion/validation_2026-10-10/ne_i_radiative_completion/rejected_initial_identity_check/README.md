# Preserved initial identity-check failure

The first Ne I controller stopped in its pilot step before writing any
native or direct population profile. Its `run_profiles.py` SHA-256 is
`f8a80639815056774fb5f4e97bec4a943887df1f74a83b085fd7406a583a5ee2`.
The script demanded exact equality between CHIANTI-selected energies and the
Stout `.nrg` entries for all 31 levels. Direct inspection showed the same
level indices and statistical weights, but Stout rounds three energy entries
by 0.001 cm^-1: levels 20, 25, and 28. The three levels used by the added
A values, 3, 5, and 15, agree exactly. The exact-float check was therefore
too strict for the published source precision; it did not expose a different
level identity.

The original plan, controller state and log, and exact rejected script are
preserved here. `../run_profiles_v2.py` requires matching indices and exact
statistical weights, permits no more than 0.002 cm^-1 energy difference,
and records the measured maximum difference. The separate
`../recovery_v2/plan.json` does not overwrite this failed controller state.
