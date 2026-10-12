# Installation of 23 missing supplied finite-level atoms

Twenty-three previously absent native atom files were installed under
`data/atomic`, preserving exact validated candidate bytes and every existing
file. [Installation record](installation/result.json) lists the ions, levels,
source bounds, and hashes. These additions do not install exponential fits or
replace any shared-ion atom, and do not activate the expanded rate chains.
Mg I and Ca VII remain uninstalled pending their source limitations.

## Invoked physics and file checks

The actual compiled reader opened all 23 installed locations and evaluated an
original collision node. Maximum relative source-node difference was
8.881784197001252e-16. The existing complete native/source excitation profiles
were preserved; unchanged numerical physics was not regenerated merely because
its file was installed at a new location.

Cl V is the only newly added model within an already active production
chain. The actual species setup and no-fit cooling consumer were exercised
using the existing five-stage chlorine rate file and the newly installed
two-level atom. The default line-scope option accepts this common v2 ion.
Two existing density modes passed 264 samples from 1,000 to 100,000 K and
four zero-density cases. Maximum relative independent-matrix coefficient
difference was 6.00158299411957e-5 (0.00600158%).

| Input condition | Linked production consumer | Independent excitation calculation |
| --- | --- | --- |
| Levels, radiation, strengths | Installed Cl V supplied model | Identical model |
| Temperature | 33 points, 1,000--100,000 K | Same points |
| Electron density (cm^-3), `local_ne` | 1e-6, 100, 1e6, 1e10 | Same densities |
| Electron density, `low_density` | Native reference at 1e-6 cm^-3 | Same reference density |
| Chlorine ionization | Existing five-stage rates, imposed Gamma | Removed when comparing the coefficient |
| Other chlorine line models | Omitted from isolated fixture | Omitted |
| Spectrum, geometry, spatial solution | Not calculated | Not calculated |

This fixture isolates the new Cl V term; it does not measure a complete
chlorine thermal solution or certify a photon band for the terminal Cl V
stage. Electron density and photoionization rates are imposed independently.
The [full cooling profiles](native/clv_installed_cooling_profiles.pdf) were
visually inspected, and [root review](root_review.json) verifies source,
installed-byte, and output hashes. The persistent controller completed both
steps; no calculation from this sequence remains live.

## Current coverage and limitations

A [new direct inventory](../../phase1_77_mocassin_v2/mochii_coverage_after_missing_installation.json)
now finds two absent atom files among 56 available MOCASSIN source models,
25 absent exponential fits, and seven smaller installed explicit atoms.
An absent fit does not suppress a supplied atom's cooling: the implemented
absolute-cooling consumer supplies that calculation in either density mode.
The active production chains still contain 43 metal stages and 32 adjacent
transitions. Upper-stage file presence is preparation, not invocation.

Published absolute atomic accuracy remains unverified as stated in each
source header. Corrected positive interpolation and adopted cold constant
strengths remain approximations. Shared-ion reference selection, complete
Mg I and corrected Ca VII data, full 77-stage activation, remaining rate and
energy controls, and matched spatial validation remain unfinished.
