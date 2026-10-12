# Absolute atomic line cooling without an exponential fit

The production `species_setup` now registers an available level atom when its
cooling fit is absent, in either existing density mode. `metal_cooling` uses
`nlevel_cooling_lam_se` for this case. Existing fitted calculations retain their
combination rule. Zero hydrogen or electron density returns zero electron-impact
line cooling before ion fractions are evaluated. Neutral-hydrogen excitation
has its own calculation and is not removed by that electron-density guard.

## Actual linked validation

The debug production objects were linked into the focused test. Fixtures use
staged seven-stage neon and iron rates, but supply only the Ne VI or Fe V atom
and no exponential cooling files. Both the default `local_ne` and explicit
`low_density` namelist choices were exercised on 33 temperatures from 1,000 to
100,000 K and four imposed electron densities: 1e-6, 100, 1e6, and 1e10 cm^-3.
All 528 samples and eight zero-density checks passed. The maximum relative
coefficient difference from independent full-level matrix solutions was
0.0012425075345752212 (0.124251%), within the 0.5% test criterion.

| Input condition | Production calculation | Independent comparison |
| --- | --- | --- |
| Level energies, A values, collision nodes | Identical candidate source atom | Same supplied finite-level atom |
| Temperature | Same 33 points | Same 33 points |
| `local_ne` density | Imposed density at each point | Same imposed density |
| `low_density` density | Table reference column at 1e-6 cm^-3 | Same reference density |
| Ion abundance | Calculated fraction and fixed abundance | Removed when comparing the cooling coefficient |
| Geometry, spectrum, transport | No spatial calculation | No spatial calculation |

The density floor is a finite-density approximation to the low-density limit,
not an exact zero-density population solution. Electron density is imposed for
this diagnostic; it is not a self-consistent charge-neutral plasma model.

Both pages of [the complete curves](native/absolute_cooling_profiles.pdf) were
visually inspected. [Root review](root_review.json) records output hashes and
numerical and visual acceptance. The controller completed its recorded steps;
no background model-resume service is claimed. These checks cover the changed
cooling consumer, not an equilibrium-temperature solver or spatial transport.

Ne VI and Fe V no longer require an exponential fit to contribute physical
line cooling. Their earlier rejected fits remain preserved. Candidate atoms
are not installed; production remains at 43 metal stages and 32 transitions.
Atomic source accuracy and expanded production integration remain unfinished.
