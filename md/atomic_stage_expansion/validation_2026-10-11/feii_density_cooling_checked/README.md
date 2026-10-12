# Fe II matching finite-density cooling continuation

The first controller terminated before calculating a table because its copied
four-argument driver was called with five arguments. The original failed
products and the exact generator are preserved in
[the first attempt](../feii_density_cooling/controller/state.json).
The corrected generator changes only this private diagnostic argument count;
the production driver still has its original four-argument interface.

The first compiler also reported that runtime debug checks disabled its
requested optimization. This continuation removes those flags from the
isolated compilation and requests `-O2 -no-ipo -fp-model precise`, retaining
floating-point exception trapping. The candidate atom equations are unchanged.
Earlier direct-profile diagnostics used the checked compilation.

The tmux controller `mochii_feii_density_checked_20261011`, on socket
`mochii-workflows`, runs one finite-density validation step. Its child timeout
is 1200 seconds; the controller step timeout is 1500 seconds. It verifies a
complete result, 40 samples, a measured 301-by-161 table, and no production
changes. It stops on failure. No model-completion event watcher is configured.
Current authoritative state is in [controller/state.json](controller/state.json).
Do not infer completion from the existence of this document.

The diagnostic limits the original temperature-grid prefix to the fit's
1000-100000 K interval. It reads the matching Fe II fit and invokes actual
`nlevel_cooling_apply`, comparing the result with direct atom cooling at
regular and intermediate temperature and density points. It does not invoke
species thermal balance, the line-output writer, or a spatial simulation.

## Completed verification

The corrected controller completed, and root verified the controller source
hashes and all native-result artifact hashes before any later comment change.
The measured table has 301 temperature nodes and 161 density nodes, totaling
48461 cells. Its registration wall time was 33.096096 seconds.
All 40 application samples passed. Maximum relative differences are
0.1173885499% for actual matching-fit application versus direct cooling and
0.0941962222% for intermediate-point table interpolation versus direct cooling.
The largest application error occurs at a high intermediate electron density;
this is an interpolation-and-fit error, not an atomic uncertainty estimate.

Twenty regular samples also matched the earlier checked, unoptimized direct
calculations to a maximum relative difference of 2.88658e-15. Changing compiler
optimization therefore did not measurably change those direct solutions.
All four pages of [the full-profile comparison](matching_density_profiles.pdf)
were visually inspected. The continuous curves are reused 401-point direct
profiles; the optimized native application is sampled at five regular points
in each density panel. The error panels also show intermediate points at the
explicitly labeled nearby densities. They are not dense native-application
profiles or spatial nebular profiles.

| Input condition | Direct comparison and matching-fit application |
| --- | --- |
| Ion and level energies | Identical 142-level Fe II candidate |
| Electron collision strengths | Identical 10011 collision pairs |
| Radiation probabilities | Identical 1144 allowed and 40 preserved forbidden pairs |
| Temperature interval | 1000-100000 K |
| Regular electron densities | 1e-6, 100, 1e6, 1e10 cm^-3 |
| Intermediate electron densities | 10^-5.95, 10^2.05, 10^6.05, 10^9.95 cm^-3 |
| Ion abundance | Coefficient normalized by electron and Fe II densities |
| Radiation spectrum and transport | No imposed photoionizing field; excitation calculation only |
| Domain, gravity, and boundaries | No spatial domain or hydrodynamic calculation |
| Other physical terms | Electron excitation and radiative de-excitation only |

The exact tested generator is retained as
[frozen_preparation_tool.py](frozen_preparation_tool.py). Its output column
header inherited the label `direct_floor`; this column actually contains the
matching fit from `cooling_eval`, as demonstrated by the retained driver.
The current generator corrects that future header to `matching_fit` without
regenerating or altering the hashed completed products.

[Native result](native/result.json) and [root review](root_review.json) record
completion and limits. No production data changed, no physical default was
accepted, and the 77-stage production objective remains incomplete.
