# Ne I radiative completion and source-node validation

The former 31-level Ne I atom omitted two identified spontaneous decays into
level 15, making that excited level falsely closed. The corrected atom retains
all 31 selected energies, 142 original CHIANTI radiative rates, 465 collision
pairs, and 6,975 collision nodes. It adds only the Stout rates
`A(3,15) = 9.0e5 s^-1` and `A(5,15) = 6.82e7 s^-1`. Its 144 radiative rows use
the inferred `BT_TMIN_K = 1000.36 K` continuation: collision strengths are
constant below that source temperature, while Maxwellian rates use the actual
gas temperature. This floor is inferred from finite, rounded source nodes;
published accuracy over 1,000–100,000 K has not been established.

The [candidate result](candidate_result.json) records structural checks and
the [candidate atom](atom/nlevel_ne_1.txt), SHA-256
`95b241d6a27b9f42fda9a39c895a9145f28825da9d9a68908196c8826d6e15d4`.
The successful [recovery controller](recovery_v2/controller/state.json) froze
its source files and completed all seven numerical and product stages. The
initial, overly strict Stout/CHIANTI energy equality check stopped before any
population calculation; its unchanged files and explanation remain in
[rejected_initial_identity_check](rejected_initial_identity_check/README.md).
The corrected identity check matches indices and statistical weights and
allows 0.002 cm^-1 energy rounding. The measured maximum is
0.0010000000184 cm^-1 at levels 20, 25, and 28; levels 3, 5, and 15 match
exactly.

For 401 temperatures from 1,000 to 100,000 K and electron densities
`1e-6`, `100`, `1e6`, and `1e10 cm^-3`, the optimized native source reader and
an independently reconstructed Python solution from raw CHIANTI and the two
Stout rates agree in cooling to at most `2.22e-15` relative and in populations
to at most `5.56e-16` absolute. The raw-source collision-strength comparison
is within `2.22e-16` relative. [Complete cooling profiles](cooling_profiles.pdf),
[native differences](native_differences.pdf), and
[all level populations](all_level_populations.pdf) were visually inspected.
The independent [balance ledger](independent_balance_result.json), evaluated
from saved native populations without rerunning either solver, has maximum
stationarity residual `1.09e-15` of gross collisional flow and maximum
collision-to-radiation energy residual `4.93e-15` of gross collisional power.
See the [validation manifest](validation_manifest.json) for source hashes,
condition-specific maxima, and the explicit physical accuracy limitation.

The [private cooling product result](../ne_i_cooling_products/candidate/result.json)
passes its 801-temperature fit (maximum relative error `3.18e-4`), actual
compiled coefficient-reader check (`7.78e-12`), and 64,561-cell registered
table check (grid `1.15e-14`, midpoint `7.13e-5`, application `4.15e-5` maximum
relative errors). The [fit plot](../ne_i_cooling_products/candidate/fit/NeI_fit.pdf)
was visually inspected: the source and fitted cooling curves overlap throughout
the range, with a smooth residual of approximately `0.032%` or less.

The main-session review [accepted and installed](main_review.json) the exact
validated [atom and fit bytes](../ne_i_cooling_products/installation_manifest.json).
An installed debug reader at 10,000 K and four densities reproduced every
population and cooling record exactly. This is numerical and source-identity
validation of the adopted model, not a verification of the underlying
collision-strength accuracy or a spatial plasma calculation. The historical
31-level, 142-rate atom is physically stale because level 15 lacked two
identified radiative losses.
