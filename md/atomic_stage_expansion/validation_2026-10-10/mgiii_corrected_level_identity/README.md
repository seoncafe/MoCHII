# Mg III corrected-level private validation

The raw CHIANTI Mg III level 167 is labeled `2s2 2p5 7p 1S`, `J = 0`.
Its raw observed-energy field is 628148.375 cm^-1, but that energy belongs
to a NIST `7d`, `J = 1` level. The selected energy for the same CHIANTI
level identity is therefore its original theoretical value, 634430 cm^-1.
The [source-identity audit](../mgiii_radiative_candidate/source_identity/README.md)
records the official NIST response and the other 67 higher observed-energy
matches. This is an explicit theoretical-energy approximation, not a claim
that all high-lying wavefunctions or collision strengths are accurate.

`prepare_candidate.py` used the production source-node writer to generate
`atom/nlevel_mg_3.txt` (SHA-256
`ca27317e9a76235bcd45243c61ab4659fcfb562928c917f25606c4dd4197764f`).
It checked all 167 selected levels, 3,364 radiative rows, 13,861 collision
rows, 207,915 BT nodes, and the 1,800 K cold Upsilon bound. All selected
energies except level 167, all collision records and nodes, and the first-27
atom are unchanged relative to the preserved rejected candidate. Twenty-four
radiative rows involving level 167 changed; a direct reconstruction from the
raw `.wgfa` rates agrees exactly with the serialized corrected rates. The
focused Mg III electric-dipole test also passed.

The finite controller plan `plan.json` freezes 40 source, code, object, and
input hashes at start. It runs a small full-atom pilot, optimized compiled
full-atom profiles, independent Python full-atom profiles, comparison with
preserved first-27 profiles, an independent stationarity and energy ledger,
and a validation manifest. Only after these numerical checks pass does it
run the private cooling-fit and registered-table controller. The outer
controller's state is in `controller/state.json`; the nested product state
is under `../mgiii_cooling_products/controller/`. The sequence is private
and does not install production atomic data.

The numerical interval is 1,000--100,000 K at electron densities of
`1e-6`, `1e2`, `1e6`, and `1e10 cm^-3`. Upsilon is held constant below its
published lower temperature of 1,800 K, while Maxwellian rates and detailed
balance use actual gas temperature. The later products fit this declared
cold continuation. Numerical agreement and energy conservation do not
establish the accuracy of the underlying source radiation and collision
data. The final publisher article for the high-level Mg III radiation has
not been verified here.

The full-atom native and independent Python comparison passed at all 801
temperatures and four densities. The maximum relative cooling difference was
`2.58e-14`, and the maximum absolute population difference was `2.00e-15`.
The native collision strengths agreed with 69,305 direct source samples to
`2.66e-15` relative. From the saved native populations, the maximum level
stationarity residual relative to total gross collisional flow was
`3.49e-15`; collisional energy input minus radiative output relative to gross
collisional energy flow reached `4.14e-14`. The most weakly populated upper
states can have a large **local** relative residual at 1,000 K when their
saved population is zero, so the ledger also reports the local residual for
levels carrying more than `1e-12` of the gross flow. These conditions and
hashes are recorded in `independent_balance_result.json`.

The first 27 levels differ from the corrected full atom by up to `17.19%`
in cooling, so that truncation fails the `0.3%` criterion within this
candidate model. The [complete cooling profiles](full_and_first27_cooling_profiles.pdf)
and [relative-difference profiles](first27_truncation_differences.pdf) were
visually inspected. The earlier candidate with the misassigned level-167
energy is physically rejected, regardless of its maximum `0.0687%` cooling
difference from the corrected candidate on this grid.

Both persistent controllers completed all guarded steps. The private
low-density cooling fit has a maximum relative error of `0.2107%` over the
801-temperature grid. The actual serialized-coefficient reader agreed to
`1.14e-11` relative at all 801 resolved temperatures. The 64,561-cell
registered table validation reached `7.11e-15` relative error at grid
points, `4.15e-5` between grid points, and `5.35e-6` for the applied cooling
value. The [fit figure](../mgiii_cooling_products/candidate/fit/MgIII_fit.pdf)
was visually inspected. The full private product record is
[`candidate/result.json`](../mgiii_cooling_products/candidate/result.json).
These numerical checks do not authorize production installation or certify
the physical accuracy of the source calculation.
