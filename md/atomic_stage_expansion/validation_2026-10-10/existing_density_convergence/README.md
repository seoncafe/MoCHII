# Existing low-ion density convergence, 2026-10-10

The installed five-level C I, C II, N I, N II, and S II atoms do not converge
to their complete collision-connected bound sources over the tested interval.
The omitted excited levels carry substantial cooling, especially toward
100,000 K and high electron density. Applying the complete bound atom with
`DENSITY_RATIO 1` removes the spurious positive residual that the partial-atom
rule leaves after the low-density fit and five-state source nearly cancel.
This conclusion concerns optically thin, electron-impact line cooling, not a
spatial photoionization benchmark or the accuracy of the original scattering
calculations.

| Quantity | Complete-source and selected-source calculations | Installed atom calculation |
| --- | --- | --- |
| Electron temperature | 801 logarithmic samples, 1,000-100,000 K | Identical grid |
| Electron density | 1e-6, 1e2, 1e6, 1e10 cm^-3 | Identical grid |
| Radiation field | No photoexcitation or trapping | None in the isolated atom driver |
| Population balance | Source collision strengths, spontaneous emission, Maxwellian detailed balance | Frozen optimized native n-level solver, installed five-state Chebyshev atoms |
| Energy selection | `chianti_cooling.bound_excitation_atom`, including the ionization threshold | Values serialized in the installed atom |
| Cooling | Sum of spontaneous photon energy, divided by electron density | Same native emitted-power quantity |
| Atomic hydrogen | Absent from electron-density profiles | Separate full C I/C II sidecar check; the original first-level identities and energies are unchanged |

The complete-source calculation and the native comparison both use four
densities at every temperature. Numerical interpolation of the registered
401-by-161 cooling table was checked separately at its nodes and interior
points. Native validation includes source collision strengths, populations,
and line-emission power. The complete profiles and differences were visually
inspected; the figure links below contain every sampled temperature and
density. Raw source, code, driver, and input hashes are in `result.json`,
`n_s_result.json`, the native result files, and the workflow state files.

| Ion | Complete levels | Radiative rows | Collision rows | Largest five-state source deficit at 1e-6 cm^-3 | Largest five-state source deficit over all four densities | Largest native cooling error with complete atom | Largest registered interior interpolation error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| C I | 42 | 243 | 861 | 92.678% | 99.9968% | 3.53e-13 | 0.1689% |
| C II | 53 | 314 | 1,378 | 93.175% | 97.8201% | 3.55e-15 | 0.0965% |
| N I | 50 | 388 | 559 | 83.926% | 99.9974% | 4.66e-15 | 0.0837% |
| N II | 58 | 386 | 327 | 93.000% | 99.9939% | 4.44e-15 | 0.1372% |
| S II | 70 | 619 | 2,415 | 91.865% | 99.9972% | 5.44e-15 | 0.1402% |

Every complete atom passed the compiled source-node comparison and actual
cooling registration with `DENSITY_RATIO 1`. At the four sampled densities,
the remaining difference between the complete-source density ratio multiplied
by the installed low-density fit and the source cooling is exactly the fit
error: 0.9485% for C I, 0.2959% for C II, 0.5894% for N I, 0.6625% for N II,
and 0.3793% for S II. The existing partial-atom expression instead has
unbounded amplification when the fit and five-state low-density reference
nearly cancel; the measured largest relative errors at 1e10 cm^-3 are
3.31e6, 1.31e6, 4.10e3, 5.88e5, and 5.61e3, respectively, in the same ion
order. These numbers are relative errors, not percent values. The comparison
used the native installed-atom profile in the physical complete/partial
equations; the registered-table interpolation has its own measured error.

The original C I and C II atomic-hydrogen collision sidecars use source level
IDs 1-3 and 1-2, respectively. Those level IDs, weights, and selected energies
are identical in the full atoms. A private build of the existing H-impact gate
with the full atoms and unchanged sidecars passed. The main-session independent
debug-object gate records its stronger rate, Boltzmann, and emitted-power checks
in `hcoll_main_review/result.json`. The full-atom sidecars must remain beside
their atoms when installed; the isolated electron-density profile directories
did not include hydrogen impacts.

At 1,000 and 100,000 K, none of the source Burgess-Tully coordinates for these
five ions lies outside its tabulated scaled-node interval. This establishes
mathematical interpolation coverage, not published scattering accuracy at all
temperatures. C I has one exact-degenerate collision pair (source IDs 37 and
39), and C II has one (29 and 30); these collisions redistribute population
at zero energy and are retained. No complete-source radiative row has a zero
or reversed selected-energy gap for these five ions.

O III and O II require a distinct source-energy decision before claiming that
their complete bound atoms are physical. The selected observed/theoretical
energy mixture reverses 48 and 31 positive radiative rows, respectively, while
the original theoretical gaps of all those rows are positive. For the 48 O III
and 30 O II reversed rows with positive oscillator strength, the source
`gf`, `A`, and theoretical wavelength imply the *original* upper-level
statistical weight within 1%; this holds for 48/48 and 30/30 rows. Simply
reversing the upper identity while retaining `A` does not preserve the
source's transition physics. The remaining O II row has zero `gf` and cannot
be assessed by that relation. The row-by-row source and selected energies,
wavelengths, rates, and weights are preserved in
[`oxygen_radiation_audit.json`](oxygen_radiation_audit.json). No O II or O III
full-atom cooling was accepted or installed on this evidence.

Complete-profile figures: [C I](CI/cooling_profiles_layout.pdf),
[C II](CII/cooling_profiles_layout.pdf), [N I](NI/cooling_profiles_layout.pdf),
[N II](NII/cooling_profiles_layout.pdf), and
[S II](SII/cooling_profiles_layout.pdf). Source-versus-compiled atom plots
are in each ion's `native_full` directory. Earlier diagnostic figures remain
preserved under their original filenames.
