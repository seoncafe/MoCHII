# Paired O I hydrogen-excitation option check

The actual neutral-hydrogen thermal increment and the emitted-energy
increment agreed in 18 samples with the paired seven-level O I atom and
matching cooling fit. The input option was absent and explicitly enabled;
absence returned a zero increment, and enabling it changed emitted energy.
The diagnostic uses 1000, 3000, and 10000 K, ne=1e-4, 1, and 1e4 cm^-3,
and nHI=1000 cm^-3 with the neutral stage selected. These temperatures remain
within the candidate cooling interval. It exercises the actual production
routines and input reader in isolated fixtures; production was not modified.
No spatial calculation, independent atomic accuracy verification, or effect
size in a nebular model is claimed. The controller completed and root verified
recorded source hashes. See [result](native/result.json) and
[root review](root_review.json).
