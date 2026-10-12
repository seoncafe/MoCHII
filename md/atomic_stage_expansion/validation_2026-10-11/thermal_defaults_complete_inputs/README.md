# Optional thermal terms: default selection

`charge_exchange_heating` now defaults to `.false.`. `metal_freefree` remains
`.true.` at the user's explicit request on October 11. Both input settings
can be changed independently. Their thermal conditions already guard the
production `net_rate` calls. Charge transfer in ionization balance remains
independent of the kinetic-energy option.

The inspected MOCASSIN v2 `update_mod.f90::thermBalance` includes H and He
free-free and recombination terms, metal forbidden-line cooling,
photoionization heating, and optional dust terms. It does not include metal
free-free cooling or a charge-transfer kinetic-energy term. This statement
concerns the invoked thermal balance, not all radiative emission routines.

The controller completed a rebuilt debug executable and 108 native thermal
cases: three temperatures, three ionizing rates, metal-only/helium-only/combined
reaction consumers, and four input choices. At fixed temperature, the flags
leave the solved electron density and nonselected energy terms unchanged.
The resulting net-energy differences match the selected transfer heat and
metal free-free loss. The diagnostic compilation adds only `public ::
net_rate` to a temporary copy of the complete production thermal source;
it does not change a production API or copy a simplified physics routine.
No spatial calculation or new atomic-accuracy claim is made.

The first diagnostic attempt omitted `data_dir` and failed before numerical
checks. Its products remain in `../thermal_defaults/`. The corrected fixture
uses the installed Gaunt-factor data. Root verified all frozen source/product
and controller artifact hashes; see [root_review.json](root_review.json)
and [native/checks_result.json](native/checks_result.json).
